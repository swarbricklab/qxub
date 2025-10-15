#!/bin/bash
# Demo script showing qxub monitor with longer jobs to see periodic updates
# Creates jobs that will definitely show state transitions

echo "=== qxub monitor Live Updates Demo ===" >&2
echo "" >&2

# Activate environment
source venv/bin/activate

# Clear any existing job list
> live_demo_jobs.txt

echo "Creating jobs with staggered completion times..." >&2

# Job 1: 30 seconds
echo "  Job 1: 30 seconds..." >&2
job1=$(qxub --terse --name "job_30s" --default --cmd 'echo "30s job starting at $(date)"; for i in {1..6}; do echo "30s job progress: $i/6"; sleep 5; done; echo "30s job completed at $(date)"')
echo $job1 >> live_demo_jobs.txt

# Job 2: 50 seconds
echo "  Job 2: 50 seconds..." >&2
job2=$(qxub --terse --name "job_50s" --default --cmd 'echo "50s job starting at $(date)"; for i in {1..10}; do echo "50s job progress: $i/10"; sleep 5; done; echo "50s job completed at $(date)"')
echo $job2 >> live_demo_jobs.txt

# Job 3: 70 seconds
echo "  Job 3: 70 seconds..." >&2
job3=$(qxub --terse --name "job_70s" --default --cmd 'echo "70s job starting at $(date)"; for i in {1..14}; do echo "70s job progress: $i/14"; sleep 5; done; echo "70s job completed at $(date)"')
echo $job3 >> live_demo_jobs.txt

# Job 4: 40 seconds, then fails
echo "  Job 4: 40 seconds, then fails..." >&2
job4=$(qxub --terse --name "job_40s_fail" --default --cmd 'echo "40s fail job starting at $(date)"; for i in {1..8}; do echo "40s fail job progress: $i/8"; sleep 5; done; echo "40s fail job encountered error!" >&2; exit 1')
echo $job4 >> live_demo_jobs.txt

echo "" >&2
echo "Jobs created with staggered completion:" >&2
echo "  job_30s: completes at ~30s" >&2
echo "  job_50s: completes at ~50s" >&2
echo "  job_70s: completes at ~70s" >&2
echo "  job_40s_fail: fails at ~40s" >&2
echo "" >&2

echo "Job IDs:" >&2
cat live_demo_jobs.txt | sed 's/^/  /' >&2
echo "" >&2

echo "Starting monitor with 10-second intervals..." >&2
echo "Watch for jobs transitioning from Queued → Running → Completed!" >&2
echo "(Press Ctrl+C to interrupt if needed)" >&2
echo "" >&2

# Monitor with 10-second intervals to see the transitions
cat live_demo_jobs.txt | qxub monitor --interval 10

echo "" >&2
echo "Live demo completed!" >&2
echo "" >&2
echo "You should have seen:" >&2
echo "• Initial state: All jobs Queued (⏳)" >&2
echo "• Transition: Jobs move to Running (🔄)" >&2
echo "• Staggered completion: Jobs finish at different times" >&2
echo "• Final states: Success (✅) and Failure (❌)" >&2
echo "• Countdown timers between updates" >&2
