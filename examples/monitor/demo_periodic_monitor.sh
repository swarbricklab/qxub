#!/bin/bash
# Demo script showing qxub monitor periodic updates
# Creates jobs with different run times to demonstrate real-time monitoring

echo "=== qxub monitor Periodic Updates Demo ===" >&2
echo "" >&2

# Activate environment
source venv/bin/activate

# Clear any existing job list
> demo_job_ids.txt

echo "Creating jobs with different durations..." >&2

# Job 1: Quick job (10 seconds)
echo "  Creating quick job (10s)..." >&2
job1=$(qxub --terse --name "quick_job" --default --cmd 'echo "Quick job starting on $(hostname)"; sleep 10; echo "Quick job completed"')
echo $job1 >> demo_job_ids.txt

# Job 2: Medium job (25 seconds)
echo "  Creating medium job (25s)..." >&2
job2=$(qxub --terse --name "medium_job" --default --cmd 'echo "Medium job starting on $(hostname)"; sleep 25; echo "Medium job completed"')
echo $job2 >> demo_job_ids.txt

# Job 3: Slow job (40 seconds)
echo "  Creating slow job (40s)..." >&2
job3=$(qxub --terse --name "slow_job" --default --cmd 'echo "Slow job starting on $(hostname)"; sleep 40; echo "Slow job completed"')
echo $job3 >> demo_job_ids.txt

# Job 4: Failure job (15 seconds, then fails)
echo "  Creating failure job (15s, then fails)..." >&2
job4=$(qxub --terse --name "failure_job" --default --cmd 'echo "Failure job starting on $(hostname)"; sleep 15; echo "Failure job encountered an error!" >&2; exit 1')
echo $job4 >> demo_job_ids.txt

echo "" >&2
echo "Submitted 4 jobs:" >&2
echo "  • quick_job: 10 seconds → success" >&2
echo "  • medium_job: 25 seconds → success" >&2
echo "  • slow_job: 40 seconds → success" >&2
echo "  • failure_job: 15 seconds → failure" >&2
echo "" >&2

echo "Jobs submitted:" >&2
cat demo_job_ids.txt | sed 's/^/  /' >&2
echo "" >&2

echo "Monitoring with 8-second intervals..." >&2
echo "You should see different states as jobs complete at different times:" >&2
echo "  ⏳ Queued → 🔄 Running → ✅ Success / ❌ Failed" >&2
echo "" >&2

# Monitor the jobs with 8-second intervals to see periodic updates
cat demo_job_ids.txt | qxub monitor --interval 8

echo "" >&2
echo "Demo completed!" >&2
echo "" >&2
echo "This demonstrated:" >&2
echo "• Jobs transitioning from Queued → Running → Completed states" >&2
echo "• Periodic updates every 8 seconds with countdown timers" >&2
echo "• Mixed success (✅) and failure (❌) outcomes" >&2
echo "• Real-time job name display from PBS metadata" >&2
echo "• Automatic suffix removal from configuration" >&2
