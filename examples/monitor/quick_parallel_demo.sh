#!/bin/bash
# Quick parallel + qxub monitor demo with new features
# Usage: ./quick_parallel_demo.sh

# Submit 5 simple jobs in parallel, limit to 2 concurrent submissions
echo "Submitting 5 test jobs with parallel (max 2 concurrent)..." >&2

seq 1 5 | parallel -j 2 "
    echo 'Submitting job {}' >&2
    qxub --terse --name 'test_{}' --default -- sleep {}0
" | qxub monitor --suffix .gadi-pbs --interval 10

echo "All jobs completed!" >&2
