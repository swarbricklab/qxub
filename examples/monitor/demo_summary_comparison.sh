#!/bin/bash
# Demo script comparing regular vs summary monitor modes
# Shows how summary mode is better for large numbers of jobs

echo "=== qxub monitor: Regular vs Summary Mode Demo ===" >&2
echo "" >&2

# Activate environment
source venv/bin/activate

# Create multiple jobs to demonstrate summary mode benefits
> batch_job_ids.txt

echo "Creating 6 test jobs with various outcomes..." >&2

for i in {1..6}; do
    case $i in
        1|2|5)
            # Success jobs
            echo "  Creating success job $i..." >&2
            job=$(qxub --terse --name "success_job_$i" --default --cmd "echo 'Success job $i on $(hostname)'; sleep $((10 + RANDOM % 15)); echo 'Job $i completed successfully'")
            ;;
        3|4)
            # Failure jobs
            echo "  Creating failure job $i..." >&2
            job=$(qxub --terse --name "failure_job_$i" --default --cmd "echo 'Failure job $i on $(hostname)'; sleep $((8 + RANDOM % 10)); echo 'Job $i failed!' >&2; exit 1")
            ;;
        6)
            # Long running job
            echo "  Creating long job $i..." >&2
            job=$(qxub --terse --name "long_job_$i" --default --cmd "echo 'Long job $i on $(hostname)'; sleep 35; echo 'Long job $i finally completed'")
            ;;
    esac
    echo $job >> batch_job_ids.txt
done

echo "" >&2
echo "Submitted 6 jobs (3 success, 2 failure, 1 long-running)" >&2
echo "" >&2

# Wait a moment for jobs to start
sleep 5

echo "1. Regular mode (full table):" >&2
echo "   Shows individual job details - good for few jobs" >&2
echo "" >&2

# Show regular mode briefly
timeout 25s cat batch_job_ids.txt | qxub monitor --interval 8 || true

echo "" >&2
echo "2. Summary mode (compact counts):" >&2
echo "   Shows job state summary - perfect for many jobs" >&2
echo "" >&2

# Show summary mode
cat batch_job_ids.txt | qxub monitor --summary --interval 8

echo "" >&2
echo "Demo completed!" >&2
echo "" >&2
echo "Summary mode benefits:" >&2
echo "• Compact display regardless of job count" >&2
echo "• Easy to see overall progress at a glance" >&2
echo "• No terminal scrolling with many jobs" >&2
echo "• Perfect for pipeline monitoring" >&2
echo "• Shows success/failure breakdown" >&2
