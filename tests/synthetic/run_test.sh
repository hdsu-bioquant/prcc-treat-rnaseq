#run_test.sh
rm -rf tests/synthetic/results

snakemake \
  --snakefile workflow/Snakefile \
  --configfile tests/synthetic/config.yaml \
  --cores 2 \
  --software-deployment-method apptainer

python tests/synthetic/validate_results.py