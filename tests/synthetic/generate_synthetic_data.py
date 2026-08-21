#!/usr/bin/env python3
from pathlib import Path
import gzip, hashlib, random

ROOT = Path(__file__).resolve().parents[2]
TEST = ROOT / 'tests' / 'synthetic'
DATA = TEST / 'data'
EXPECTED = TEST / 'expected'
REF = TEST / 'reference'
for p in (DATA, EXPECTED, REF):
    p.mkdir(parents=True, exist_ok=True)

SEED = 20250124
REF_LEN = 8000
QUAL = 'I'

# ---------------------------------------------------------------------------
# Deterministic synthetic reference. Avoid homopolymers >= 6, making it less
# likely that the QuantSeq poly(A)-trim test is confounded by reference sequence.
# ---------------------------------------------------------------------------
rng = random.Random(SEED)
bases = []
for i in range(REF_LEN):
    choices = list('ACGT')
    if len(bases) >= 5 and len(set(bases[-5:])) == 1:
        choices.remove(bases[-1])
    bases.append(rng.choice(choices))
refseq = ''.join(bases)

# gene_id -> strand, exons (1-based inclusive, ascending genomic order), name
GENES = {
    'SYN_GENE_A': {'strand': '+', 'exons': [(1001,1300),(1601,1900)], 'name': 'SYN_A'},
    'SYN_GENE_B': {'strand': '-', 'exons': [(3501,3800),(4101,4400)], 'name': 'SYN_B'},
    'SYN_GENE_C': {'strand': '+', 'exons': [(6001,6500)],             'name': 'SYN_C'},
}

def rc(s):
    return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]

def gslice(start, end):
    return refseq[start-1:end]

def transcript(gid):
    g = GENES[gid]
    exons = g['exons'] if g['strand'] == '+' else list(reversed(g['exons']))
    seqs = [gslice(s,e) for s,e in exons]
    if g['strand'] == '-':
        seqs = [rc(x) for x in seqs]
    return ''.join(seqs)

TX = {g: transcript(g) for g in GENES}
assert len(TX['SYN_GENE_A']) == 600
assert len(TX['SYN_GENE_B']) == 600
assert len(TX['SYN_GENE_C']) == 500

# Check uniqueness of 50-mers in reference (+ reverse complement). This makes
# exact short-read placement unambiguous for the unspliced synthetic reads.
kmers = {}
for i in range(len(refseq)-49):
    k = refseq[i:i+50]
    kmers[k] = kmers.get(k,0)+1
assert max(kmers.values()) == 1

# Write reference FASTA.
with open(REF/'synthetic.fa', 'w') as fh:
    fh.write('>chrSynthetic pRCC-TREAT deterministic miniature RNA-seq test reference\n')
    for i in range(0, len(refseq), 80):
        fh.write(refseq[i:i+80]+'\n')

# Write a compact but STAR/HTSeq-compatible GTF.
with open(REF/'synthetic.gtf', 'w') as fh:
    fh.write('##description: pRCC-TREAT deterministic miniature RNA-seq test annotation\n')
    fh.write('##seed: %d\n' % SEED)
    for gid, g in GENES.items():
        gene_start = min(s for s,e in g['exons'])
        gene_end = max(e for s,e in g['exons'])
        tid = gid.replace('GENE','TX') + '.1'
        attrs_gene = f'gene_id "{gid}"; gene_type "protein_coding"; gene_name "{g["name"]}";'
        attrs_tx = attrs_gene + f' transcript_id "{tid}"; transcript_type "protein_coding"; transcript_name "{g["name"]}-001";'
        fh.write(f'chrSynthetic\tpRCC-TREAT\tgene\t{gene_start}\t{gene_end}\t.\t{g["strand"]}\t.\t{attrs_gene}\n')
        fh.write(f'chrSynthetic\tpRCC-TREAT\ttranscript\t{gene_start}\t{gene_end}\t.\t{g["strand"]}\t.\t{attrs_tx}\n')
        exon_order = g['exons'] if g['strand'] == '+' else list(reversed(g['exons']))
        for n,(s,e) in enumerate(exon_order,1):
            attrs_ex = attrs_tx + f' exon_number "{n}"; exon_id "{tid}.E{n}";'
            fh.write(f'chrSynthetic\tpRCC-TREAT\texon\t{s}\t{e}\t.\t{g["strand"]}\t.\t{attrs_ex}\n')

# FASTQ helpers. Gzip mtime=0 gives stable bytes/checksums on regeneration.
def write_fastq(path, records):
    path = Path(path)
    with open(path, 'wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as gz:
            for name, seq in records:
                text = f'@{name}\n{seq}\n+\n{QUAL*len(seq)}\n'
                gz.write(text.encode())

def pe_fragment(gid, start, frag_len=220, r1_len=100, r2_len=100):
    tx = TX[gid]
    frag = tx[start:start+frag_len]
    assert len(frag) == frag_len
    return frag[:r1_len], rc(frag[-r2_len:])

def pe_fragment_umi(gid, start, umi, frag_len=180):
    tx = TX[gid]
    frag = tx[start:start+frag_len]
    assert len(frag) == frag_len
    # Raw R1 is 100 nt: 6-nt UMI followed by 94 biological bases.
    return umi + frag[:94], rc(frag[-100:])

def qs_read(gid, start, bio_len=60, polya=15):
    seq = TX[gid][start:start+bio_len]
    assert len(seq) == bio_len
    return seq + 'A'*polya

def qs_read_umi(gid, start, umi, bio_len=60, polya=9):
    seq = TX[gid][start:start+bio_len]
    assert len(seq) == bio_len
    return umi + seq + 'A'*polya

read_manifest = []
def add_manifest(sample, rid, gid, umi, group, exp_aligned=True, exp_trimmed=True):
    read_manifest.append((sample,rid,gid,umi,group,int(exp_aligned),int(exp_trimmed)))

# ---------------------------------------------------------------------------
# 1) Full-length paired-end, no UMI: 5/4/3 fragments over A/B/C, plus one
# deliberately unmappable pair. Several A/B fragments cross splice junctions.
# ---------------------------------------------------------------------------
fl_r1, fl_r2 = [], []
fl_starts = {
    'SYN_GENE_A': [40,120,220,320,360],
    'SYN_GENE_B': [30,150,220,350],
    'SYN_GENE_C': [20,180,260],
}
for gid, starts in fl_starts.items():
    for j, st in enumerate(starts,1):
        rid = f'FL_noUMI|{gid}|fragment{j}|start{st}'
        r1,r2 = pe_fragment(gid, st)
        fl_r1.append((rid+'/1',r1)); fl_r2.append((rid+'/2',r2))
        add_manifest('FL_noUMI', rid, gid, '-', f'{gid}:start{st}')
rid='FL_noUMI|UNMAPPED|fragment1'
fl_r1.append((rid+'/1','N'*100)); fl_r2.append((rid+'/2','N'*100))
add_manifest('FL_noUMI', rid, 'UNMAPPED', '-', 'unmapped', False, True)
write_fastq(DATA/'FL_noUMI_R1.fastq.gz', fl_r1)
write_fastq(DATA/'FL_noUMI_R2.fastq.gz', fl_r2)

# ---------------------------------------------------------------------------
# 2) Full-length paired-end with 6-nt UMI at R1 start. Deliberate PCR duplicate
# groups test that UMI handling is assay-independent after the refactor.
# Raw gene counts: A/B/C = 6/5/4; UMI molecules = 4/4/3.
# ---------------------------------------------------------------------------
fl_u_r1, fl_u_r2 = [], []
fl_umi_blueprint = {
 'SYN_GENE_A': [(100,'ACGTAC',3),(100,'TGCATG',1),(200,'ACGTAC',1),(300,'GATCGA',1)],
 'SYN_GENE_B': [(80,'CCGGTT',2),(80,'TTAACT',1),(180,'CCGGTT',1),(280,'AGCTGA',1)],
 'SYN_GENE_C': [(40,'CACGTC',2),(140,'CACGTC',1),(140,'GATCGA',1)],
}
for gid, groups in fl_umi_blueprint.items():
    for st, umi, copies in groups:
        for c in range(1,copies+1):
            rid=f'FL_UMI|{gid}|start{st}|umi{umi}|copy{c}'
            r1,r2=pe_fragment_umi(gid, st, umi)
            fl_u_r1.append((rid+'/1',r1)); fl_u_r2.append((rid+'/2',r2))
            add_manifest('FL_UMI',rid,gid,umi,f'{gid}:start{st}:umi{umi}')
rid='FL_UMI|UNMAPPED|umiTGCATG|copy1'
fl_u_r1.append((rid+'/1','TGCATG'+'N'*94)); fl_u_r2.append((rid+'/2','N'*100))
add_manifest('FL_UMI',rid,'UNMAPPED','TGCATG','unmapped',False,True)
write_fastq(DATA/'FL_UMI_R1.fastq.gz', fl_u_r1)
write_fastq(DATA/'FL_UMI_R2.fastq.gz', fl_u_r2)

# ---------------------------------------------------------------------------
# 3) QuantSeq SE, no UMI. Each normal read has 15 nt poly(A) to exercise BBDuk.
# Includes one unmapped read and one read that becomes <20 nt after poly(A) trim.
# Expected aligned gene counts A/B/C = 5/4/3.
# ---------------------------------------------------------------------------
qs = []
qs_starts = {'SYN_GENE_A':[500,505,510,515,520], 'SYN_GENE_B':[500,510,520,530], 'SYN_GENE_C':[400,410,420]}
for gid, starts in qs_starts.items():
    for j,st in enumerate(starts,1):
        rid=f'QS_noUMI|{gid}|read{j}|start{st}'
        qs.append((rid,qs_read(gid,st)))
        add_manifest('QS_noUMI',rid,gid,'-',f'{gid}:start{st}')
rid='QS_noUMI|UNMAPPED|read1'
qs.append((rid,'N'*60+'A'*15)); add_manifest('QS_noUMI',rid,'UNMAPPED','-','unmapped',False,True)
rid='QS_noUMI|TRIM_DROP|read1'
qs.append((rid,TX['SYN_GENE_A'][450:460]+'A'*65)); add_manifest('QS_noUMI',rid,'TRIM_DROP','-','trim_drop',False,False)
write_fastq(DATA/'QS_noUMI_R1.fastq.gz',qs)

# ---------------------------------------------------------------------------
# 4) QuantSeq SE + Lexogen-like 6-nt UMI at R1 start. Normal reads are 75 nt:
# 6 UMI + 60 biological + 9 poly(A). Same-coordinate/same-UMI PCR duplicates,
# same-coordinate/different-UMI molecules, and same-UMI/different-coordinate
# molecules are all represented.
# Expected raw A/B/C = 6/5/4; dedup molecules = 4/4/3.
# ---------------------------------------------------------------------------
qsu=[]
qs_umi_blueprint = {
 'SYN_GENE_A': [(500,'ACGTAC',3),(500,'TGCATG',1),(510,'ACGTAC',1),(520,'GATCGA',1)],
 'SYN_GENE_B': [(500,'CCGGTT',2),(500,'TTAACT',1),(510,'CCGGTT',1),(520,'AGCTGA',1)],
 'SYN_GENE_C': [(400,'CACGTC',2),(410,'CACGTC',1),(410,'GATCGA',1)],
}
for gid, groups in qs_umi_blueprint.items():
    for st,umi,copies in groups:
        for c in range(1,copies+1):
            rid=f'QS_UMI|{gid}|start{st}|umi{umi}|copy{c}'
            qsu.append((rid,qs_read_umi(gid,st,umi)))
            add_manifest('QS_UMI',rid,gid,umi,f'{gid}:start{st}:umi{umi}')
rid='QS_UMI|UNMAPPED|umiTGCATG|copy1'
qsu.append((rid,'TGCATG'+'N'*60+'A'*9)); add_manifest('QS_UMI',rid,'UNMAPPED','TGCATG','unmapped',False,True)
rid='QS_UMI|TRIM_DROP|umiACGTAC|copy1'
qsu.append((rid,'ACGTAC'+TX['SYN_GENE_A'][450:460]+'A'*59)); add_manifest('QS_UMI',rid,'TRIM_DROP','ACGTAC','trim_drop',False,False)
write_fastq(DATA/'QS_UMI_R1.fastq.gz',qsu)

# Expected gene lengths from union-of-exons.
with open(EXPECTED/'gene_lengths.tsv','w') as fh:
    fh.write('gene_id\tgene_length\tgene_type\tgene_name\n')
    for gid in GENES:
        length = sum(e-s+1 for s,e in GENES[gid]['exons'])
        fh.write(f'{gid}\t{length}\tprotein_coding\t{GENES[gid]["name"]}\n')

# Expected raw STAR unstranded counts, after preprocessing/alignment.
raw_counts = {
 'SYN_GENE_A': {'FL_noUMI':5,'FL_UMI':6,'QS_noUMI':5,'QS_UMI':6},
 'SYN_GENE_B': {'FL_noUMI':4,'FL_UMI':5,'QS_noUMI':4,'QS_UMI':5},
 'SYN_GENE_C': {'FL_noUMI':3,'FL_UMI':4,'QS_noUMI':3,'QS_UMI':4},
}
with open(EXPECTED/'raw_gene_counts.tsv','w') as fh:
    samples=['FL_noUMI','FL_UMI','QS_noUMI','QS_UMI']
    fh.write('gene_id\tgene_name\t'+'\t'.join(samples)+'\n')
    for gid in GENES:
        fh.write(gid+'\t'+GENES[gid]['name']+'\t'+'\t'.join(str(raw_counts[gid][s]) for s in samples)+'\n')

# Expected molecule counts for UMI-bearing samples.
dedup_counts = {
 'SYN_GENE_A': {'FL_UMI':4,'QS_UMI':4},
 'SYN_GENE_B': {'FL_UMI':4,'QS_UMI':4},
 'SYN_GENE_C': {'FL_UMI':3,'QS_UMI':3},
}
with open(EXPECTED/'umi_dedup_gene_counts.tsv','w') as fh:
    fh.write('gene_id\tgene_name\tFL_UMI\tQS_UMI\n')
    for gid in GENES:
        fh.write(f'{gid}\t{GENES[gid]["name"]}\t{dedup_counts[gid]["FL_UMI"]}\t{dedup_counts[gid]["QS_UMI"]}\n')

# Read-level blueprint for debugging failures.
with open(EXPECTED/'read_manifest.tsv','w') as fh:
    fh.write('sample\tread_id\tgene_id\tumi\tmolecule_group\texpected_to_align\texpected_to_survive_trim\n')
    for row in read_manifest:
        fh.write('\t'.join(map(str,row))+'\n')

# Basic expected input/preprocessing summary.
summary = [
 ('FL_noUMI', len(fl_r1), len(fl_r1), 12, '-', 1, 0),
 ('FL_UMI', len(fl_u_r1), len(fl_u_r1), 15, 11, 1, 0),
 ('QS_noUMI', len(qs), len(qs)-1, 12, '-', 1, 1),
 ('QS_UMI', len(qsu), len(qsu)-1, 15, 11, 1, 1),
]
with open(EXPECTED/'expected_summary.tsv','w') as fh:
    fh.write('sample\tinput_fragments_or_reads\texpected_after_trim\texpected_raw_gene_assigned\texpected_umi_molecules\texpected_unmapped\texpected_trim_dropped\n')
    for row in summary:
        fh.write('\t'.join(map(str,row))+'\n')

# Future-facing sample sheet: only adds UMI metadata to the present schema.
with open(TEST/'samples.tsv','w') as fh:
    fh.write('sample\tassay\tfq1\tfq2\tpatient\tbatch\tstrandedness\thas_umi\tumi_pattern\tumi_location\n')
    rows=[
      ('FL_noUMI','full_length_pe','tests/synthetic/data/FL_noUMI_R1.fastq.gz','tests/synthetic/data/FL_noUMI_R2.fastq.gz','SYN_FL1','synthetic','forward','false','-','-'),
      ('FL_UMI','full_length_pe','tests/synthetic/data/FL_UMI_R1.fastq.gz','tests/synthetic/data/FL_UMI_R2.fastq.gz','SYN_FL2','synthetic','forward','true','NNNNNN','read1_start'),
      ('QS_noUMI','quantseq_3prime_se','tests/synthetic/data/QS_noUMI_R1.fastq.gz','-','SYN_QS1','synthetic','forward','false','-','-'),
      ('QS_UMI','quantseq_3prime_se','tests/synthetic/data/QS_UMI_R1.fastq.gz','-','SYN_QS2','synthetic','forward','true','NNNNNN','read1_start'),
    ]
    for row in rows: fh.write('\t'.join(row)+'\n')

# A target config for the refactored smoke-test mode. It intentionally points to
# a local miniature reference rather than GDC downloads.
config = '''\
# pRCC-TREAT synthetic smoke test (target configuration)
samples: tests/synthetic/samples.tsv
results: tests/synthetic/results
tmpdir: tests/synthetic/results/tmp

reference:
  mode: local
  dir: tests/synthetic/reference
  download_references: false
  genome_fasta: synthetic.fa
  gtf: synthetic.gtf
  star_index: star_index
  build_star_index: true
  sjdb_overhang: 100
  # Tiny genomes need a smaller STAR suffix-array seed than the default.
  genome_sa_index_nbases: 4

star:
  threads: 2
  gdc_params: >-
    --alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000
    --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignSoftClipAtReferenceEnds Yes
    --chimJunctionOverhangMin 15 --chimMainSegmentMultNmax 1
    --chimOutType Junctions SeparateSAMold WithinBAM SoftClip --chimSegmentMin 15
    --genomeLoad NoSharedMemory --limitSjdbInsertNsj 1200000
    --outFilterIntronMotifs None --outFilterMatchNminOverLread 0.33
    --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1
    --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33
    --outFilterType BySJout --outSAMattributes NH HI AS nM NM ch
    --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --outSAMunmapped Within
    --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic

full_length:
  trim_adapters: false
  count_column: unstranded
  compute_fpkm_tpm: true

quantseq:
  bbduk_polyA: true

modules:
  fusion: false
  te: false
  ase: false
  rseqc: false
'''
(TEST/'config.yaml').write_text(config)


# ---------------------------------------------------------------------------
# Checksums
# ---------------------------------------------------------------------------
# The checksum manifest fingerprints only the canonical synthetic test fixture:
# config, sample sheet, FASTQs, miniature reference, and golden expected TSVs.
# Documentation and test-running code are deliberately excluded.
CHECKSUM_FILES = [
    TEST / 'config.yaml',
    TEST / 'samples.tsv',

    DATA / 'FL_noUMI_R1.fastq.gz',
    DATA / 'FL_noUMI_R2.fastq.gz',
    DATA / 'FL_UMI_R1.fastq.gz',
    DATA / 'FL_UMI_R2.fastq.gz',
    DATA / 'QS_noUMI_R1.fastq.gz',
    DATA / 'QS_UMI_R1.fastq.gz',

    REF / 'synthetic.fa',
    REF / 'synthetic.gtf',

    EXPECTED / 'expected_summary.tsv',
    EXPECTED / 'gene_lengths.tsv',
    EXPECTED / 'raw_gene_counts.tsv',
    EXPECTED / 'read_manifest.tsv',
    EXPECTED / 'umi_dedup_gene_counts.tsv',
]

checksum_path = TEST / 'checksums.sha256'
with checksum_path.open('w') as fh:
    for p in CHECKSUM_FILES:
        if not p.is_file():
            raise FileNotFoundError(f'Cannot checksum missing fixture file: {p}')
        digest = hashlib.sha256(p.read_bytes()).hexdigest()
        fh.write(f'{digest}  {p.relative_to(ROOT)}\n')

# Validate gzip FASTQs and lengths/counts.
def read_fastq(path):
    out=[]
    with gzip.open(path,'rt') as fh:
        while True:
            h=fh.readline().rstrip('\n')
            if not h: break
            seq=fh.readline().rstrip('\n'); plus=fh.readline().rstrip('\n'); qual=fh.readline().rstrip('\n')
            assert h.startswith('@') and plus.startswith('+') and len(seq)==len(qual)
            out.append((h[1:],seq))
    return out

assert len(read_fastq(DATA/'FL_noUMI_R1.fastq.gz')) == 13
assert len(read_fastq(DATA/'FL_noUMI_R2.fastq.gz')) == 13
assert len(read_fastq(DATA/'FL_UMI_R1.fastq.gz')) == 16
assert len(read_fastq(DATA/'FL_UMI_R2.fastq.gz')) == 16
assert len(read_fastq(DATA/'QS_noUMI_R1.fastq.gz')) == 14
assert len(read_fastq(DATA/'QS_UMI_R1.fastq.gz')) == 17
# All QS raw reads are exactly 75 nt; FL no-UMI = 100; FL UMI raw R1/R2 = 100.
assert {len(s) for _,s in read_fastq(DATA/'QS_noUMI_R1.fastq.gz')} == {75}
assert {len(s) for _,s in read_fastq(DATA/'QS_UMI_R1.fastq.gz')} == {75}
assert {len(s) for _,s in read_fastq(DATA/'FL_noUMI_R1.fastq.gz')} == {100}
assert {len(s) for _,s in read_fastq(DATA/'FL_UMI_R1.fastq.gz')} == {100}

print('Synthetic test fixture built successfully at', TEST)
print(f'Checksum manifest written for {len(CHECKSUM_FILES)} canonical fixture files: {checksum_path}')
print('FASTQ sizes:')
for p in sorted(DATA.glob('*.fastq.gz')):
    print(f'  {p.name}: {p.stat().st_size} bytes')
