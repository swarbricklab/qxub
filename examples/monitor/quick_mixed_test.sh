#!/bin/bash
# Quick test of qxub monitor with success/failure scenarios
# Uses shorter wait times for rapid testing

echo "Quick demo: Mixed job outcomes with qxub monitor" >&2

# Submit a few quick jobs with known outcomes
job_ids=()

echo "Submitting jobs..." >&2

# Job 1: Will succeed (exit 0)
echo "  Job 1: Success scenario" >&2
job1=$(qxub --terse --name "success_job" --default -- bash -c "echo 'Success job running on node \$HOSTNAME'; sleep 5; echo 'Success job completed'; exit 0")
job_ids+=($job1)

# Job 2: Will fail (exit 1)
echo "  Job 2: Failure scenario" >&2
job2=$(qxub --terse --name "failure_job" --default -- bash -c "echo 'Failure job running on node \$HOSTNAME'; sleep 3; echo 'Failure job encountered an error!' >&2; exit 1")
job_ids+=($job2)

# Job 3: Will succeed but take longer
echo "  Job 3: Slow success scenario" >&2
job3=$(qxub --terse --name "slow_success" --default -- bash -c "echo 'Slow job running on node \$HOSTNAME'; sleep 8; echo 'Slow job completed successfully'; exit 0")
job_ids+=($job3)

echo "" >&2
echo "Monitoring ${#job_ids[@]} jobs with mixed outcomes..." >&2
echo "Expected results:" >&2
echo "  success_job: ✅ (exit code 0)" >&2
echo "  failure_job: ❌ (exit code 1)" >&2
echo "  slow_success: ✅ (exit code 0)" >&2
echo "" >&2

# Monitor the jobs
printf '%s\n' "${job_ids[@]}" | qxub monitor --suffix .gadi-pbs --interval 5

# Show the exit code
exit_code=$?
echo "" >&2
echo "Monitor exit code: $exit_code" >&2
if [ $exit_code -eq 0 ]; then
    echo "Result: All jobs succeeded ✅" >&2
else
    echo "Result: Some jobs failed ❌" >&2
fi

echo "" >&2
echo "This demonstrates:" >&2
echo "• Mixed success/failure job outcomes" >&2
echo "• Proper exit code tracking (✅ vs ❌)" >&2
echo "• Monitor exit code reflects overall success" >&2
echo "• Suffix removal (.gadi-pbs hidden)" >&2
echo "• Job name display for context" >&2
