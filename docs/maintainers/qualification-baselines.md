# Qualification baseline maintenance

This is a maintainer procedure for reviewing and updating the versioned deterministic validation
baselines used by the synthetic and realistic qualification tests. Ordinary qualification compares
new results with the baselines shipped in the repository and does not regenerate them.

The repository maintains four distinct integrity/validation manifests:

- `tests/synthetic/checksums.sha256` — integrity of the version-controlled synthetic **input
  fixture**;
- `tests/synthetic/expected/validation_checksums.sha256` — expected deterministic synthetic
  **outputs**;
- `tests/real/expected/validation_checksums.sha256` — expected deterministic realistic
  **outputs**; and
- `resources/gdc_installed_reference.sha256` — expected bytes of the maintained installed GDC
  reference bundle.

Do not refresh these together merely because one qualification baseline changes. In particular,
changing expected pipeline outputs does not imply a change to the synthetic input-fixture checksum
manifest or to the canonical GDC reference manifest.

See `release-policy.md` for the validation level appropriate to a repository change.

## Default qualification logic

The maintained expected checksums are the first reproducibility check. Start with one ordinary
qualification run against the existing baseline.

- If the generated validation manifest matches the maintained baseline, qualification is complete;
  no second run is required merely to reconfirm an unchanged expected result.
- If the manifest differs, investigate the changed deterministic files before updating anything.
- If the difference is unintended or unexplained, fix or investigate the implementation and retain
  the existing baseline.
- If the difference is intentional and limited to deterministic metadata/provenance, review the
  exact diff and update the affected baseline when appropriate. A second computational run is not
  required solely to prove that a directly generated metadata value changed as intended.
- If deterministic computational outputs change, finalize the intended implementation first, then
  require two fresh clean runs of the affected qualification under that finalized implementation.
  Their generated `validation_checksums.sha256` files must agree byte-for-byte before the maintained
  expected baseline is replaced.

This makes the maintained baseline the normal guard against drift while reserving repeated runs for
cases where a new computational result actually needs a new reproducibility baseline.

## Baseline update helper

Use the maintainer helper to inspect or accept deterministic output baselines. It previews the diff
by default and changes a maintained baseline only with `--apply`.

It intentionally does **not** update `tests/synthetic/checksums.sha256` or
`resources/gdc_installed_reference.sha256`.

### Synthetic

After an ordinary synthetic run reports an understood expected-baseline mismatch, the completed
workflow output can be validated again while skipping only the frozen-baseline comparison:

```bash
python tests/synthetic/validate_results.py --skip-frozen-baseline
```

This avoids a second workflow execution when the original run completed successfully and failed only
on comparison with the old expected manifest. If the workflow or any other validation failed,
investigate that failure instead of entering baseline-update mode.

Preview the candidate baseline:

```bash
bash tests/maintainers/update_validation_baseline.sh synthetic
```

If only reviewed deterministic metadata/provenance changed, the candidate may be accepted after
that review. If deterministic computational outputs changed, first produce two fresh clean
synthetic runs with the finalized implementation and preserve/compare their generated manifests,
for example:

```bash
bash tests/synthetic/run_test.sh --skip-frozen-baseline
cp tests/synthetic/output/results/run/validation_checksums.sha256 /tmp/prcc-synthetic-run1.sha256

bash tests/synthetic/run_test.sh --skip-frozen-baseline
cmp /tmp/prcc-synthetic-run1.sha256 \
    tests/synthetic/output/results/run/validation_checksums.sha256
```

Only after the required review/reproducibility check succeeds, update the expected baseline:

```bash
bash tests/maintainers/update_validation_baseline.sh synthetic --apply
```

Then rerun ordinary qualification against the newly maintained baseline:

```bash
bash tests/synthetic/run_test.sh
```

The helper updates only `tests/synthetic/expected/validation_checksums.sha256`.

## Realistic qualification baseline

Run the realistic qualification normally first. If its generated deterministic manifest differs
from the maintained baseline, investigate the difference before entering baseline-update mode.

For a deliberate candidate baseline, run the maintained validator with only the frozen-baseline
comparison skipped:

```bash
python tests/real/validate_results.py \
  --run-dir "$RUN_DIR" \
  --skip-frozen-baseline
```

Preview the candidate against the maintained realistic baseline:

```bash
bash tests/maintainers/update_validation_baseline.sh realistic "$RUN_DIR"
```

If the reviewed difference is limited to deterministic metadata/provenance, the candidate may be
accepted without a second workflow run. If deterministic computational outputs changed, first
produce two new clean realistic runs from the finalized implementation using the same pinned public
FASTQs, a qualified GDC installation, and the site's production-style persistent execution profile.
For example, preserve the first candidate manifest:

```bash
cp "$RUN_DIR/output/results/run/validation_checksums.sha256" \
   "$RUN_DIR/validation_checksums.run1.sha256"
```

Remove the first output tree, repeat the workflow and `validate_results.py --skip-frozen-baseline`
from the same run definition, and compare the second candidate:

```bash
cmp "$RUN_DIR/validation_checksums.run1.sha256" \
    "$RUN_DIR/output/results/run/validation_checksums.sha256" \
  && echo "PASS: realistic deterministic validation manifest reproduced"
```

If the two manifests differ, investigate nondeterminism or unintended behavior and do not update
the baseline.

After the required review/reproducibility check succeeds, update the maintained realistic baseline:

```bash
bash tests/maintainers/update_validation_baseline.sh realistic "$RUN_DIR" --apply
```

Then rerun the ordinary validator on the completed run without `--skip-frozen-baseline`:

```bash
python tests/real/validate_results.py --run-dir "$RUN_DIR"
```

A new workflow execution is not required merely to perform this final comparison; the validator
checks the completed run against the newly maintained expected manifest.

## Canonical GDC installed-reference manifest

`resources/gdc_installed_reference.sha256` defines the expected installed bytes for the maintained
GDC reference bundle. It is separate from pipeline output baselines and should change only when the
maintained reference bundle intentionally changes.

Generate a maintainer candidate with:

```bash
bash resources/verify_gdc_references.sh \
  --write-canonical-manifest "$GDC_DIR"
```

Review the candidate before committing it. Then qualify the site installation against the reviewed
manifest:

```bash
bash resources/verify_gdc_references.sh --qualify "$GDC_DIR"
```

Sites qualify their installation against the repository manifest rather than creating a separate
site-specific canonical baseline.
