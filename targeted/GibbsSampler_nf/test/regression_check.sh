#!/usr/bin/env bash
#
# regression_check.sh
#
# Confirms that the patched engine (src/Deep_seq_tree_GS.jl) produces output
# identical to the baseline engine (test/engine_baseline.jl) on a real sample.
# This is the acceptance check for the performance fixes: same seed, same
# parameters, only the engine differs, so identical output proves the changes
# did not alter the sampler.
#
# Requirements: a Julia with Phylo, RCall, DataFrames, Distributions, CodecZlib
# and ArgParse available, plus R with ape, in the environment you run this from.
#
# Usage: test/regression_check.sh <path_to_gibbs_info.RDS> [sample_id]

set -euo pipefail

RDS="${1:?provide the path to a *_gibbs_info.RDS file}"
SID="${2:-regtest}"

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
OUT="$HERE/out"
rm -rf "$OUT"
mkdir -p "$OUT/patched" "$OUT/baseline"

# A short run is sufficient to exercise both blocks and the recording step.
COMMON=(--input "$RDS" --sample-id "$SID" --seed 42 \
        --iter 2000 --burn-in 1000 --thin 100 --scale-pm 50)

echo "Running patched engine..."
julia "$ROOT/bin/run_gibbs.jl" "${COMMON[@]}" \
    --outdir "$OUT/patched" --engine "$ROOT/src/Deep_seq_tree_GS.jl"

echo "Running baseline engine..."
julia "$ROOT/bin/run_gibbs.jl" "${COMMON[@]}" \
    --outdir "$OUT/baseline" --engine "$HERE/engine_baseline.jl"

echo "Comparing outputs..."
ok=1
for f in "${SID}_posterior_VAFs.txt.gz" "${SID}_branch_VAFs.txt.gz"; do
    if diff <(zcat "$OUT/patched/$f") <(zcat "$OUT/baseline/$f") > /dev/null; then
        echo "  IDENTICAL: $f"
    else
        echo "  DIFFERENT: $f"
        ok=0
    fi
done

if [ "$ok" -eq 1 ]; then
    echo "PASS: patched and baseline outputs are identical on this input."
else
    echo "FAIL: outputs differ; the performance changes are not output-neutral here."
    exit 1
fi
