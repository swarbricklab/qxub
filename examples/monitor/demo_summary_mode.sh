#!/bin/bash
# Demo script showing the difference between regular and summary mode
# when monitoring many jobs

echo "qxub monitor --summary demonstration" >&2
echo "=====================================" >&2
echo "" >&2

# Use known job IDs from previous tests
job_ids=(
    "152366620" "152366621" "152366622"  # 3 jobs from quick test
    "152387053" "152387054" "152387055"  # 3 more jobs
    "152387166" "152387167" "152387168"  # Latest 3 jobs
)

echo "Monitoring ${#job_ids[@]} jobs..." >&2
echo "" >&2

echo "1. Regular mode (full table):" >&2
printf '%s\n' "${job_ids[@]}" | qxub monitor --suffix .gadi-pbs --quiet

echo "" >&2
echo "2. Summary mode (compact counts):" >&2
printf '%s\n' "${job_ids[@]}" | qxub monitor --suffix .gadi-pbs --summary

echo "" >&2
echo "Benefits of --summary mode:" >&2
echo "• Compact display for many jobs" >&2
echo "• Easy to see overall progress at a glance" >&2
echo "• No terminal scrolling issues" >&2
echo "• Quick state overview with counts" >&2
