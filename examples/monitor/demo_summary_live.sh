#!/bin/bash
# Quick demo showing summary mode with live updates
# Perfect for demonstrating the --summary option

echo "=== qxub monitor --summary Live Demo ===" >&2
echo "" >&2

# Activate environment
source venv/bin/activate

# Create some quick test jobs
> quick_summary_jobs.txt

echo "Creating 5 quick jobs for summary demo..." >&2

for i in {1..5}; do
    duration=$((15 + i * 5))  # 20, 25, 30, 35, 40 seconds
    if [ $i -eq 3 ]; then
        # Make job 3 fail
        echo "  Job $i: ${duration}s (will fail)..." >&2
        job=$(qxub --terse --name "summary_job_$i" --default --cmd "echo 'Summary job $i starting'; sleep $duration; echo 'Job $i failed!' >&2; exit 1")
    else
        echo "  Job $i: ${duration}s..." >&2
        job=$(qxub --terse --name "summary_job_$i" --default --cmd "echo 'Summary job $i starting'; sleep $duration; echo 'Job $i completed'")
    fi
    echo $job >> quick_summary_jobs.txt
done

echo "" >&2
echo "Monitoring with --summary mode (10s intervals)..." >&2
echo "Summary mode shows compact state counts instead of individual rows:" >&2
echo "" >&2

# Monitor in summary mode
cat quick_summary_jobs.txt | qxub monitor --summary --interval 10

echo "" >&2
echo "Summary mode demo completed!" >&2
echo "" >&2
echo "Benefits of --summary mode:" >&2
echo "• Compact display shows job state counts" >&2
echo "• Perfect for monitoring many jobs (10s, 100s, 1000s)" >&2
echo "• No screen scrolling issues" >&2
echo "• Clear success/failure breakdown" >&2
echo "• Updates every interval with state transitions" >&2
