#!/bin/bash
# Comprehensiv    case $job_num in
        1|2)
            # These jobs will succeed (exit code 0)
            echo "  Submitting successful job $job_num" >&2
            qxub --terse \
                --name "$job_name" \
                --env "$CONDA_ENV" \
                --joblog "${job_name}.log" \
                -- \
                bash -c "echo 'Job $job_num: Starting successful task on node \$HOSTNAME'; sleep \$((RANDOM % 20 + 10)); echo 'Job $job_num: Task completed successfully'; exit 0"
            ;;
        3|4)
            # These jobs will fail (exit code 1)
            echo "  Submitting failing job $job_num" >&2
            qxub --terse \
                --name "$job_name" \
                --env "$CONDA_ENV" \
                --joblog "${job_name}.log" \
                -- \
                bash -c "echo 'Job $job_num: Starting task that will fail on node \$HOSTNAME'; sleep \$((RANDOM % 15 + 5)); echo 'Job $job_num: ERROR - Task failed!' >&2; exit 1"
            ;;
        5)
            # This job will have random outcome
            echo "  Submitting random outcome job $job_num" >&2
            qxub --terse \
                --name "$job_name" \
                --env "$CONDA_ENV" \
                --joblog "${job_name}.log" \
                -- \
                bash -c "echo 'Job $job_num: Starting random outcome task on node \$HOSTNAME'; sleep 15; outcome=\$((RANDOM % 2)); if [ \$outcome -eq 0 ]; then echo 'Job $job_num: Random outcome - SUCCESS'; exit 0; else echo 'Job $job_num: Random outcome - FAILURE' >&2; exit 1; fi"
            ;;
        6)
            # This job will be long-running but successful
            echo "  Submitting long-running job $job_num" >&2
            qxub --terse \
                --name "$job_name" \
                --env "$CONDA_ENV" \
                --joblog "${job_name}.log" \
                -- \
                bash -c "echo 'Job $job_num: Starting long task on node \$HOSTNAME'; for i in \$(seq 1 9); do echo \"Job $job_num: Progress \$i/9\"; sleep 5; done; echo 'Job $job_num: Long task completed successfully'; exit 0"
            ;;
    esac2.3 parallel execution with success/failure scenarios
# This demonstrates the new --terse and monitor functionality with realistic job outcomes

set -e

echo "=== qxub v2.3 Parallel Execution Demo ===" >&2
echo "This demo shows:" >&2
echo "• Parallel job submission with --terse output" >&2
echo "• Mixed success/failure scenarios" >&2
echo "• Enhanced monitor display with exit code tracking" >&2
echo "• Suffix removal and emoji status indicators" >&2
echo "" >&2

# Configuration
MAX_JOBS=6
CONDA_ENV="dvc3"  # Change this to an available environment

echo "Creating $MAX_JOBS test jobs with mixed outcomes..." >&2

# Function to submit a single test job
submit_test_job() {
    local job_num=$1
    local job_name="demo_job_${job_num}"

    # Create different scenarios
    case $job_num in
        1|2)
            # These jobs will succeed (exit code 0)
            echo "  Submitting successful job $job_num" >&2
            qxub --terse \
                --name "$job_name" \
                --env "$CONDA_ENV" \
                --joblog "${job_name}.log" \
                -- \
                bash -c "
                    echo 'Job $job_num: Starting successful task on node \$HOSTNAME'
                    sleep \$((RANDOM % 20 + 10))  # 10-30 seconds
                    echo 'Job $job_num: Task completed successfully'
                    exit 0
                "
            ;;
        3|4)
            # These jobs will fail (exit code 1)
            echo "  Submitting failing job $job_num" >&2
            qxub --terse \
                --name "$job_name" \
                --env "$CONDA_ENV" \
                --joblog "${job_name}.log" \
                -- \
                bash -c "
                    echo 'Job $job_num: Starting task that will fail on node \$HOSTNAME'
                    sleep \$((RANDOM % 15 + 5))   # 5-20 seconds
                    echo 'Job $job_num: ERROR - Simulated failure occurred!' >&2
                    exit 1
                "
            ;;
        5)
            # This job will have a random outcome
            echo "  Submitting random outcome job $job_num" >&2
            qxub --terse \
                --name "$job_name" \
                --env "$CONDA_ENV" \
                --joblog "${job_name}.log" \
                -- \
                bash -c "
                    echo 'Job $job_num: Starting random outcome task on node \$HOSTNAME'
                    sleep \$((RANDOM % 15 + 8))   # 8-23 seconds
                    outcome=\$((RANDOM % 2))      # 0 or 1
                    if [ \$outcome -eq 0 ]; then
                        echo 'Job $job_num: Random success!'
                        exit 0
                    else
                        echo 'Job $job_num: Random failure!' >&2
                        exit 1
                    fi
                "
            ;;
        6)
            # This job will succeed but take longer
            echo "  Submitting long-running successful job $job_num" >&2
            qxub --terse \
                --name "$job_name" \
                --env "$CONDA_ENV" \
                --joblog "${job_name}.log" \
                -- \
                bash -c "
                    echo 'Job $job_num: Starting long task on node \$HOSTNAME'
                    for i in {1..5}; do
                        echo \"Job $job_num: Progress step \$i/5\"
                        sleep \$((RANDOM % 8 + 5))  # 5-13 seconds per step
                    done
                    echo 'Job $job_num: Long task completed successfully'
                    exit 0
                "
            ;;
    esac
}

# Export function for parallel
export -f submit_test_job
export CONDA_ENV

echo "" >&2
echo "Submitting jobs in parallel (max 3 concurrent submissions)..." >&2

# Submit all jobs in parallel and pipe to monitor
seq 1 $MAX_JOBS | parallel -j 3 submit_test_job | qxub monitor --suffix .gadi-pbs --interval 8

# Capture the exit code from monitor
monitor_exit_code=$?

echo "" >&2
echo "=== Demo Results ===" >&2

if [ $monitor_exit_code -eq 0 ]; then
    echo "🎉 All jobs completed successfully!" >&2
    echo "   Monitor exit code: $monitor_exit_code (success)" >&2
else
    echo "⚠️  Some jobs failed or encountered errors." >&2
    echo "   Monitor exit code: $monitor_exit_code (indicating failures)" >&2
fi

echo "" >&2
echo "Job logs created:" >&2
ls -la demo_job_*.log 2>/dev/null || echo "No job logs found yet" >&2

echo "" >&2
echo "=== Features Demonstrated ===" >&2
echo "✅ --terse output: Only job IDs emitted during submission" >&2
echo "✅ Parallel submission: GNU parallel with concurrency control" >&2
echo "✅ Enhanced monitor: Job names, exit code tracking, emoji status" >&2
echo "✅ Suffix removal: .gadi-pbs suffix hidden from display" >&2
echo "✅ Exit code propagation: Monitor returns 0 only if ALL jobs succeed" >&2
echo "✅ Mixed outcomes: Success (✅), failure (❌), and random scenarios" >&2
echo "✅ Real-time updates: Live status display with countdown timers" >&2

echo "" >&2
echo "Demo completed! This showcases qxub v2.3's pipeline capabilities." >&2
