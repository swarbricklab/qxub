#!/bin/bash
# Simple example combining GNU parallel with qxub --terse and qxub monitor
# This demonstrates parallel job submission and monitoring

set -e

# Configuration
DATA_DIR="sample_data"
OUTPUT_DIR="results"
CONDA_ENV="myenv"
MAX_CONCURRENT_JOBS=4

echo "=== Parallel Job Submission with qxub Example ===" >&2

# Create sample data if it doesn't exist
if [[ ! -d "$DATA_DIR" ]]; then
    echo "Creating sample data directory..." >&2
    mkdir -p "$DATA_DIR"
    for i in {01..08}; do
        echo "Sample data file $i" > "$DATA_DIR/file_${i}.txt"
        echo "Processing this will take some time..." >> "$DATA_DIR/file_${i}.txt"
    done
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Function to submit a single processing job
submit_job() {
    local input_file="$1"
    local basename=$(basename "$input_file" .txt)
    local output_file="$OUTPUT_DIR/${basename}_processed.txt"

    echo "  Submitting job for $input_file" >&2

    # Submit job with qxub --terse (outputs only job ID)
    qxub --terse \
        --name "process_${basename}" \
        --env "$CONDA_ENV" \
        --joblog "${basename}.log" \
        -- \
        bash -c "
            echo 'Processing $input_file on node \$HOSTNAME'
            sleep \$((RANDOM % 30 + 10))  # Simulate 10-40 seconds of work
            wc -l '$input_file' > '$output_file'
            echo 'Completed processing $input_file'
        "
}

# Export function and variables for parallel
export -f submit_job
export OUTPUT_DIR CONDA_ENV

echo "Submitting jobs in parallel (max $MAX_CONCURRENT_JOBS concurrent)..." >&2

# Method 1: Using find with parallel
echo "Method 1: Using find + parallel + qxub monitor" >&2
find "$DATA_DIR" -name "*.txt" | \
    parallel -j "$MAX_CONCURRENT_JOBS" submit_job | \
    qxub monitor

echo "" >&2
echo "All jobs completed!" >&2

# Alternative Method 2: Collect job IDs first, then monitor
echo "" >&2
echo "Method 2: Collect job IDs, then monitor" >&2

job_ids=()
for file in "$DATA_DIR"/*.txt; do
    if [[ -f "$file" ]]; then
        basename=$(basename "$file" .txt)
        echo "  Submitting job for $file" >&2

        job_id=$(qxub --terse \
            --name "alt_${basename}" \
            --env "$CONDA_ENV" \
            --joblog "alt_${basename}.log" \
            -- \
            bash -c "
                echo 'Alternative processing $file on node \$HOSTNAME'
                sleep \$((RANDOM % 20 + 5))  # Simulate 5-25 seconds of work
                echo 'Alternative completed $file'
            ")

        job_ids+=("$job_id")
    fi
done

echo "Monitoring ${#job_ids[@]} jobs..." >&2
printf '%s\n' "${job_ids[@]}" | qxub monitor

echo "" >&2
echo "All alternative jobs completed!" >&2

# Show results
echo "" >&2
echo "=== Results ===" >&2
if [[ -d "$OUTPUT_DIR" ]]; then
    echo "Processed files:" >&2
    ls -la "$OUTPUT_DIR/" >&2
else
    echo "No output directory found" >&2
fi

echo "" >&2
echo "Job logs:" >&2
ls -la *.log 2>/dev/null || echo "No job logs found" >&2

echo "" >&2
echo "Example completed successfully!" >&2
