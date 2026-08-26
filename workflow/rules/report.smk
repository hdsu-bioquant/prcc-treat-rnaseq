# report.smk — portable run metadata, restricted originals, manifests/checksums

import csv
import datetime
import os
import shutil
import yaml

from build_software_versions import build_software_versions


def _portable_config():
    ref = config.get("reference", {})
    modules = config.get("modules", {})
    portable = {
        "pipeline": {
            "name": "pRCC-RNA-Seq",
            "release": config.get("pipeline_release", "unversioned"),
            "output_contract": 1,
        },
        "reference": {
            "mode": ref.get("mode", "gdc"),
            "genome_fasta": os.path.basename(str(ref.get("genome_fasta", ""))),
            "gtf": os.path.basename(str(ref.get("gtf", ""))),
            "star_index": os.path.basename(str(ref.get("star_index", "")).rstrip("/")),
            "sjdb_overhang": ref.get("sjdb_overhang", 100),
        },
        "star": {
            "gdc_params": config.get("star", {}).get("gdc_params", ""),
            "primary_gene_count": "unstranded",
        },
        "full_length": {
            "trim_adapters": config.get("full_length", {}).get("trim_adapters", False),
            "normalized_expression": ["fpkm", "fpkm_uq", "tpm"],
        },
        "quantseq": {
            "bbduk_polyA": config.get("quantseq", {}).get("bbduk_polyA", True),
            "gene_length_normalization": False,
        },
        "modules": {
            "fusion": bool(modules.get("fusion", False)),
            "te": bool(modules.get("te", False)),
            "ase": bool(modules.get("ase", False)),
            "rseqc": bool(modules.get("rseqc", False)),
        },
    }
    if modules.get("fusion", False) and config.get("ctat_genome_lib"):
        portable["modules"]["ctat_genome_lib"] = os.path.basename(config["ctat_genome_lib"].rstrip("/"))
    if modules.get("te", False) and config.get("te_gtf"):
        portable["modules"]["te_gtf"] = os.path.basename(config["te_gtf"])
    if modules.get("ase", False) and config.get("ase_germline_vcf"):
        portable["modules"]["ase_germline_vcf"] = os.path.basename(config["ase_germline_vcf"])
    return portable


rule run_metadata:
    input:
        samples = config["samples"],
        fasta = FASTA,
        gtf = GTF,
        star_index = STAR_IDX_DONE
    output:
        libraries = join(RESULTS, "run/libraries.tsv"),
        portable_config = join(RESULTS, "run/config.yaml"),
        provenance = join(RESULTS, "run/provenance.yaml"),
        references = join(RESULTS, "run/references.tsv"),
        original_libraries = join(RESTRICTED, "run/libraries.original.tsv"),
        effective_config = join(RESTRICTED, "run/config.effective.yaml")
    run:
        os.makedirs(os.path.dirname(output.libraries), exist_ok=True)
        os.makedirs(os.path.dirname(output.original_libraries), exist_ok=True)

        # Restricted copies retain the original run paths and all extra metadata.
        shutil.copyfile(input.samples, output.original_libraries)
        with open(output.effective_config, "w") as fh:
            yaml.safe_dump(dict(config), fh, sort_keys=False, default_flow_style=False)

        # Portable technical library manifest: no absolute input paths and no
        # arbitrary extra metadata columns are copied implicitly.
        cols = [
            "library_id", "sample_id", "assay", "layout", "strandedness",
            "has_umi", "umi_pattern", "umi_location", "umi_discard_bases",
        ]
        samples[cols].to_csv(output.libraries, sep="\t", index=False)

        with open(output.portable_config, "w") as fh:
            yaml.safe_dump(_portable_config(), fh, sort_keys=False, default_flow_style=False)

        provenance = {
            "pipeline": "pRCC-RNA-Seq",
            "pipeline_release": config.get("pipeline_release", "unversioned"),
            "output_contract": 1,
            "generated_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "library_count": len(LIBRARIES),
            "biological_sample_count": int(samples["sample_id"].nunique()),
            "primary_expression_measure": "STAR GeneCounts unstranded raw counts",
        }
        with open(output.provenance, "w") as fh:
            yaml.safe_dump(provenance, fh, sort_keys=False, default_flow_style=False)

        ref = config.get("reference", {})
        mode = ref.get("mode", "gdc")
        with open(output.references, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(["role", "file", "reference_mode", "source"])
            if mode == "gdc":
                genome_source = "NCI GDC GRCh38.d1.vd1"
                gtf_source = "NCI GDC GENCODE v36"
                index_source = "NCI GDC STAR 2.7.5c pre-built index"
            else:
                genome_source = "local"
                gtf_source = "local"
                index_source = "local_build"
            writer.writerow(["genome_fasta", os.path.basename(FASTA), mode, genome_source])
            writer.writerow(["annotation_gtf", os.path.basename(GTF), mode, gtf_source])
            writer.writerow(["star_index", os.path.basename(STAR_IDX), mode, index_source])


rule software_versions:
    input:
        manifest = SOFTWARE_MANIFEST
    output:
        tsv = join(RESULTS, "run/software_versions.tsv"),
        multiqc = temp(join(INTERMEDIATE, "qc/prcc_pipeline_mqc_versions.yml"))
    run:
        build_software_versions(
            manifest_path=input.manifest,
            used_tools=used_tool_names(),
            resolved_containers={name: IMG[name] for name in used_tool_names()},
            output_tsv=output.tsv,
            output_multiqc_yml=output.multiqc,
        )


def manifest_inputs(wc):
    out = []
    out += expand(join(RESULTS, "libraries/{lib}/gene_expression.tsv"), lib=LIBRARIES)
    out += expand(join(RESULTS, "libraries/{lib}/qc_metrics.tsv"), lib=LIBRARIES)
    out += [
        join(RESULTS, "matrices/raw_gene_counts.tsv"),
        join(RESULTS, "qc/qc_metrics.tsv"),
        join(RESULTS, "qc/multiqc_report.html"),
        join(RESULTS, "run/libraries.tsv"),
        join(RESULTS, "run/config.yaml"),
        join(RESULTS, "run/provenance.yaml"),
        join(RESULTS, "run/software_versions.tsv"),
        join(RESULTS, "run/references.tsv"),
    ]
    if UMI_LIBRARIES:
        out.append(join(RESULTS, "matrices/umi_molecule_counts.tsv"))
    if config.get("modules", {}).get("rseqc", False):
        out += expand(join(RESULTS, "libraries/{lib}/{lib}.rseqc_read_distribution.txt"), lib=FL_LIBRARIES)
    return out


rule result_manifests:
    input:
        manifest_inputs
    output:
        manifest = join(RESULTS, "run/manifest.tsv"),
        checksums = join(RESULTS, "run/checksums.sha256"),
        validation = join(RESULTS, "run/validation_checksums.sha256")
    params:
        results = RESULTS,
        script = join(SCRIPT_DIR, "build_manifest.py")
    container: IMG["py"]
    resources:
        mem_mb = 4000, runtime = 60
    shell:
        "python {params.script:q} {params.results:q} {output.manifest:q} {output.checksums:q} {output.validation:q}"
