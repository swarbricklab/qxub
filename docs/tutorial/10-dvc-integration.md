# DVC Integration: Using qxub in Data Science Pipelines

DVC (Data Version Control) is perfect for managing data science workflows, and qxub integrates seamlessly as the execution engine for DVC stages. This section shows how to use qxub in `dvc.yaml` pipelines, handle environment switching, and coordinate parallel execution.

## Why qxub + DVC?

DVC provides:
- **Reproducible pipelines** with dependency tracking
- **Data versioning** and artifact management
- **Pipeline orchestration** with automatic stage execution

qxub adds:
- **HPC execution** with proper resource allocation
- **Real-time monitoring** of long-running stages
- **Environment management** for different pipeline stages
- **Automatic cleanup** and error handling

Together, they create powerful, reproducible HPC data science workflows.

## Basic DVC Pipeline with qxub

### Simple Pipeline Structure

Create a basic DVC pipeline that uses qxub for execution:

```yaml
# dvc.yaml
stages:

  data_preprocessing:
    cmd: qxub py -- python3 src/preprocess.py
    deps:
      - src/preprocess.py
      - data/raw/input.csv
    outs:
      - data/processed/clean_data.csv

  feature_engineering:
    cmd: qxub py --mem 8GB -- python3 src/features.py
    deps:
      - src/features.py
      - data/processed/clean_data.csv
    outs:
      - data/features/feature_matrix.csv

  model_training:
    cmd: qxub py --mem 16GB --ncpus 4 --walltime 2:00:00 -- python3 src/train.py
    deps:
      - src/train.py
      - data/features/feature_matrix.csv
    outs:
      - models/trained_model.pkl
    metrics:
      - metrics/training_metrics.json

  model_evaluation:
    cmd: qxub py -- python3 src/evaluate.py
    deps:
      - src/evaluate.py
      - models/trained_model.pkl
      - data/features/feature_matrix.csv
    outs:
      - results/predictions.csv
    metrics:
      - metrics/evaluation_metrics.json
```

### Creating the Pipeline Scripts

Let's create the supporting scripts:

```bash
# Create directory structure
mkdir -p src data/{raw,processed,features} models results metrics

# Sample preprocessing script
cat > src/preprocess.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys

print("Starting data preprocessing...")

# Simulate reading raw data
print("Reading raw data...")
np.random.seed(42)
data = pd.DataFrame({
    'feature1': np.random.randn(1000),
    'feature2': np.random.randn(1000),
    'target': np.random.randint(0, 2, 1000)
})

print(f"Raw data shape: {data.shape}")

# Simulate preprocessing
print("Cleaning data...")
data_clean = data.dropna()
data_clean = data_clean[data_clean['feature1'] > -3]
data_clean = data_clean[data_clean['feature2'] < 3]

print(f"Clean data shape: {data_clean.shape}")

# Save processed data
data_clean.to_csv('data/processed/clean_data.csv', index=False)
print("Preprocessing completed successfully!")
EOF

# Sample feature engineering script
cat > src/features.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import numpy as np

print("Starting feature engineering...")

# Read processed data
data = pd.read_csv('data/processed/clean_data.csv')
print(f"Input data shape: {data.shape}")

# Create new features
print("Engineering features...")
features = data.copy()
features['feature1_squared'] = features['feature1'] ** 2
features['feature2_log'] = np.log(np.abs(features['feature2']) + 1)
features['interaction'] = features['feature1'] * features['feature2']

print(f"Feature matrix shape: {features.shape}")

# Save feature matrix
features.to_csv('data/features/feature_matrix.csv', index=False)
print("Feature engineering completed successfully!")
EOF

# Sample training script
cat > src/train.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import json
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score

print("Starting model training...")

# Read features
features = pd.read_csv('data/features/feature_matrix.csv')
print(f"Feature matrix shape: {features.shape}")

# Prepare data
X = features.drop('target', axis=1)
y = features['target']

X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
print(f"Training set: {X_train.shape}, Validation set: {X_val.shape}")

# Train model
print("Training Random Forest...")
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Evaluate on validation set
y_pred = model.predict(X_val)
accuracy = accuracy_score(y_val, y_pred)
f1 = f1_score(y_val, y_pred)

print(f"Validation Accuracy: {accuracy:.4f}")
print(f"Validation F1-Score: {f1:.4f}")

# Save model
with open('models/trained_model.pkl', 'wb') as f:
    pickle.dump(model, f)

# Save metrics
metrics = {
    'accuracy': float(accuracy),
    'f1_score': float(f1),
    'n_estimators': 100,
    'training_samples': len(X_train)
}

with open('metrics/training_metrics.json', 'w') as f:
    json.dump(metrics, f, indent=2)

print("Model training completed successfully!")
EOF

# Sample evaluation script
cat > src/evaluate.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import pickle
import json
from sklearn.metrics import classification_report, confusion_matrix

print("Starting model evaluation...")

# Load model
with open('models/trained_model.pkl', 'rb') as f:
    model = pickle.load(f)

# Load features
features = pd.read_csv('data/features/feature_matrix.csv')
X = features.drop('target', axis=1)
y = features['target']

print(f"Evaluating on {len(X)} samples...")

# Make predictions
predictions = model.predict(X)
probabilities = model.predict_proba(X)

# Create predictions dataframe
pred_df = pd.DataFrame({
    'true_label': y,
    'predicted_label': predictions,
    'probability_class_0': probabilities[:, 0],
    'probability_class_1': probabilities[:, 1]
})

# Save predictions
pred_df.to_csv('results/predictions.csv', index=False)

# Calculate evaluation metrics
report = classification_report(y, predictions, output_dict=True)
conf_matrix = confusion_matrix(y, predictions)

# Save evaluation metrics
eval_metrics = {
    'accuracy': float(report['accuracy']),
    'precision_class_0': float(report['0']['precision']),
    'recall_class_0': float(report['0']['recall']),
    'f1_class_0': float(report['0']['f1-score']),
    'precision_class_1': float(report['1']['precision']),
    'recall_class_1': float(report['1']['recall']),
    'f1_class_1': float(report['1']['f1-score']),
    'confusion_matrix': conf_matrix.tolist()
}

with open('metrics/evaluation_metrics.json', 'w') as f:
    json.dump(eval_metrics, f, indent=2)

print(f"Final Accuracy: {eval_metrics['accuracy']:.4f}")
print("Model evaluation completed successfully!")
EOF

chmod +x src/*.py
```

### Running the DVC Pipeline

```bash
# Initialize DVC in your project
dvc init

# Run the entire pipeline
dvc repro
```

**What happens:**
- DVC executes each stage in dependency order
- Each `qxub` command runs as an HPC job
- You see real-time output from each stage
- DVC tracks all outputs and metrics
- Pipeline stops if any stage fails

## Environment Switching Between Stages

Different pipeline stages often need different software environments:

```yaml
# dvc.yaml with environment switching
stages:

  data_download:
    cmd: qxub exec --env aws -- python3 src/download_data.py
    outs:
      - data/raw/dataset.csv

  quality_control:
    cmd: qxub exec --env dvc3 --mem 8GB -- python3 src/quality_check.py
    deps:
      - src/quality_check.py
      - data/raw/dataset.csv
    outs:
      - data/qc/qc_report.html
      - data/processed/clean_dataset.csv

  genomic_analysis:
    cmd: qxub exec --env pysam --mem 32GB --ncpus 8 -- python3 src/genomic_analysis.py
    deps:
      - src/genomic_analysis.py
      - data/processed/clean_dataset.csv
    outs:
      - results/genomic_results.csv

  visualization:
    cmd: qxub exec --env tidyverse --mem 4GB -- Rscript src/create_plots.R
    deps:
      - src/create_plots.R
      - results/genomic_results.csv
    outs:
      - plots/analysis_plots.pdf

  report_generation:
    cmd: qxub exec --env pandoc -- pandoc report.md -o report.pdf
    deps:
      - report.md
      - plots/analysis_plots.pdf
    outs:
      - report.pdf
```

Each stage automatically uses the appropriate software environment for its task.

## Resource Scaling Based on Data Size

Use different resource allocations for different stages:

```yaml
# dvc.yaml with resource scaling
stages:

  small_data_prep:
    cmd: qxub test -- python3 src/prep_small.py
    deps:
      - src/prep_small.py
      - data/small_input.csv
    outs:
      - data/small_processed.csv

  medium_analysis:
    cmd: qxub py --mem 16GB --ncpus 4 -- python3 src/analyze_medium.py
    deps:
      - src/analyze_medium.py
      - data/medium_input.csv
    outs:
      - results/medium_results.csv

  large_scale_processing:
    cmd: qxub bigmem -- python3 src/process_large.py
    deps:
      - src/process_large.py
      - data/large_input.csv
    outs:
      - results/large_results.csv

  memory_intensive_ml:
    cmd: qxub exec --mem 128GB --ncpus 16 --walltime 8:00:00 --queue hugemem -- python3 src/train_large_model.py
    deps:
      - src/train_large_model.py
      - results/large_results.csv
    outs:
      - models/large_model.pkl
```

## Parallel Execution Within DVC Stages

### Simple Parallel Stage with Blocking

```yaml
# dvc.yaml with parallel execution
stages:

  parallel_analysis:
    cmd: >
      bash -c '
      # Submit parallel jobs
      job_ids=()
      for sample in sample1 sample2 sample3 sample4; do
        job_id=$(qxub exec --terse --name "analyze-$sample" sc --mem 16GB -- python3 src/analyze_sample.py --sample $sample)
        job_ids+=($job_id)
        echo "Submitted $sample: $job_id"
      done

      # Wait for all jobs to complete
      echo "Waiting for ${#job_ids[@]} parallel jobs to complete..."
      qxub monitor --wait-for-completion "${job_ids[@]}"
      echo "All parallel analysis completed"
      '
    deps:
      - src/analyze_sample.py
      - data/samples/
    outs:
      - results/sample1_results.csv
      - results/sample2_results.csv
      - results/sample3_results.csv
      - results/sample4_results.csv
```

### Advanced Parallel Stage with Error Handling

```yaml
stages:

  robust_parallel_processing:
    cmd: >
      bash -c '
      set -e  # Exit on any error

      echo "Starting robust parallel processing..."

      # Create job submission function
      submit_analysis_job() {
        local sample=$1
        local job_id=$(qxub exec --terse --name "process-$sample" \
          --env dvc3 --mem 8GB --ncpus 2 --walltime 1:00:00 \
          -- python3 src/process_sample.py --sample "$sample" --output "results/${sample}_output.csv")
        echo "$job_id"
      }

      # Submit all jobs
      job_ids=()
      for sample in $(ls data/samples/ | sed "s/.csv//"); do
        job_id=$(submit_analysis_job "$sample")
        job_ids+=($job_id)
        echo "Submitted processing for $sample: $job_id"
      done

      echo "Submitted ${#job_ids[@]} parallel processing jobs"

      # Monitor with failure detection
      if qxub monitor --wait-for-completion --fail-fast "${job_ids[@]}"; then
        echo "✅ All parallel processing completed successfully"
      else
        echo "❌ Some parallel jobs failed"
        exit 1
      fi

      # Combine results
      echo "Combining results..."
      python3 src/combine_results.py
      echo "Results combination completed"
      '
    deps:
      - src/process_sample.py
      - src/combine_results.py
      - data/samples/
    outs:
      - results/combined_results.csv
```

## Parameter Sweeps in DVC

### Parameterized Pipeline

```yaml
# dvc.yaml with parameters
stages:

  parameter_sweep:
    cmd: >
      bash -c '
      # Read parameters from params.yaml
      learning_rates=(0.001 0.01 0.1)
      batch_sizes=(32 64 128)

      job_ids=()
      for lr in "${learning_rates[@]}"; do
        for bs in "${batch_sizes[@]}"; do
          job_name="train-lr${lr}-bs${bs}"
          job_id=$(qxub exec --terse --name "$job_name" \
            --env dvc3 --mem 16GB --ncpus 4 --walltime 2:00:00 \
            -- python3 src/train_model.py --lr $lr --batch_size $bs --output_dir "models/lr${lr}_bs${bs}")
          job_ids+=($job_id)
          echo "Submitted training: lr=$lr, batch_size=$bs ($job_id)"
        done
      done

      echo "Submitted ${#job_ids[@]} parameter sweep jobs"
      qxub monitor --wait-for-completion "${job_ids[@]}"

      # Select best model
      python3 src/select_best_model.py
      '
    deps:
      - src/train_model.py
      - src/select_best_model.py
      - data/features/
    outs:
      - models/best_model.pkl
    metrics:
      - metrics/parameter_sweep_results.json
```

## Cross-Validation Pipeline

```yaml
stages:

  cross_validation:
    cmd: >
      bash -c '
      echo "Starting 5-fold cross-validation..."

      # Submit cross-validation jobs
      job_ids=()
      for fold in {1..5}; do
        job_id=$(qxub exec --terse --name "cv-fold-$fold" \
          --env dvc3 --mem 12GB --ncpus 3 --walltime 1:30:00 \
          -- python3 src/cross_validate.py --fold $fold --output_dir "cv_results/fold_$fold")
        job_ids+=($job_id)
        echo "Submitted fold $fold: $job_id"
      done

      # Wait for all folds
      echo "Waiting for cross-validation to complete..."
      qxub monitor --wait-for-completion "${job_ids[@]}"

      # Aggregate results
      echo "Aggregating cross-validation results..."
      python3 src/aggregate_cv_results.py
      echo "Cross-validation completed"
      '
    deps:
      - src/cross_validate.py
      - src/aggregate_cv_results.py
      - data/features/feature_matrix.csv
    outs:
      - cv_results/
    metrics:
      - metrics/cross_validation_metrics.json
```

## Handling Large-Scale Data Processing

### Chunked Processing Pipeline

```yaml
stages:

  chunk_and_process:
    cmd: >
      bash -c '
      echo "Processing large dataset in chunks..."

      # Split data into chunks
      python3 src/create_chunks.py --input data/large_dataset.csv --chunk_size 10000 --output_dir data/chunks/

      # Process each chunk in parallel
      job_ids=()
      for chunk_file in data/chunks/chunk_*.csv; do
        chunk_name=$(basename "$chunk_file" .csv)
        job_id=$(qxub exec --terse --name "process-$chunk_name" \
          --env dvc3 --mem 8GB --ncpus 2 --walltime 30:00 \
          -- python3 src/process_chunk.py --input "$chunk_file" --output "results/processed_$chunk_name.csv")
        job_ids+=($job_id)
        echo "Submitted $chunk_name: $job_id"
      done

      echo "Processing ${#job_ids[@]} chunks in parallel..."
      qxub monitor --wait-for-completion "${job_ids[@]}"

      # Merge processed chunks
      echo "Merging processed chunks..."
      python3 src/merge_chunks.py --input_dir results/ --output results/final_processed_data.csv
      echo "Large-scale processing completed"
      '
    deps:
      - src/create_chunks.py
      - src/process_chunk.py
      - src/merge_chunks.py
      - data/large_dataset.csv
    outs:
      - results/final_processed_data.csv
```

## Integration with DVC Experiments

### Experiment Tracking with qxub

```yaml
# dvc.yaml for experiments
stages:

  experiment:
    cmd: >
      bash -c '
      # Get experiment parameters
      model_type=${MODEL_TYPE:-"random_forest"}
      n_estimators=${N_ESTIMATORS:-100}
      max_depth=${MAX_DEPTH:-10}

      echo "Running experiment: $model_type with $n_estimators estimators, max_depth=$max_depth"

      # Choose resources based on model type
      if [ "$model_type" = "neural_network" ]; then
        alias_name="ml"  # More resources for neural networks
      else
        alias_name="py"  # Standard resources for tree models
      fi

      # Run experiment
      job_id=$(qxub exec --terse --name "exp-$model_type" $alias_name \
        -- python3 src/run_experiment.py \
          --model_type "$model_type" \
          --n_estimators "$n_estimators" \
          --max_depth "$max_depth" \
          --output_dir "experiments/exp_$(date +%Y%m%d_%H%M%S)")

      echo "Submitted experiment: $job_id"
      qxub monitor --wait-for-completion "$job_id"
      '
    deps:
      - src/run_experiment.py
      - data/features/
    outs:
      - experiments/
    metrics:
      - metrics/experiment_results.json
```

Run experiments with different parameters:

```bash
# Run different experiments
MODEL_TYPE=random_forest N_ESTIMATORS=100 dvc repro
MODEL_TYPE=random_forest N_ESTIMATORS=200 dvc repro
MODEL_TYPE=neural_network N_ESTIMATORS=50 dvc repro

# Compare experiments
dvc exp show
```

## Best Practices for DVC + qxub

### 1. Resource-Aware Stage Design

```yaml
# Good: Match resources to stage requirements
stages:
  data_exploration:
    cmd: qxub test -- python3 src/explore.py  # Quick, low resources

  feature_engineering:
    cmd: qxub py -- python3 src/features.py   # Standard resources

  model_training:
    cmd: qxub ml -- python3 src/train.py      # High resources for ML

  visualization:
    cmd: qxub r -- Rscript src/plots.R        # R environment
```

### 2. Error Handling and Debugging

```yaml
stages:
  robust_stage:
    cmd: >
      bash -c '
      set -e  # Exit on error
      set -x  # Show commands (for debugging)

      echo "Starting robust stage..."

      # Test command with dry run first (in development)
      # qxub exec --dry py -- python3 src/script.py

      job_id=$(qxub exec --terse py -- python3 src/script.py)
      echo "Job submitted: $job_id"

      if qxub monitor --wait-for-completion "$job_id"; then
        echo "Stage completed successfully"
      else
        echo "Stage failed - check job logs"
        qxub history show "$job_id" --error
        exit 1
      fi
      '
```

### 3. Reusable Stage Templates

Create reusable templates for common patterns:

```bash
# Create template script
cat > scripts/run_parallel_analysis.sh << 'EOF'
#!/bin/bash
# Reusable parallel analysis template

set -e

SCRIPT_NAME=$1
INPUT_PATTERN=$2
OUTPUT_DIR=$3
ALIAS=${4:-py}

echo "Running parallel analysis: $SCRIPT_NAME"
echo "Input pattern: $INPUT_PATTERN"
echo "Output directory: $OUTPUT_DIR"

job_ids=()
for input_file in $INPUT_PATTERN; do
    basename=$(basename "$input_file")
    job_id=$(qxub exec --terse --name "analyze-$basename" $ALIAS -- python3 "$SCRIPT_NAME" --input "$input_file" --output "$OUTPUT_DIR")
    job_ids+=($job_id)
    echo "Submitted $basename: $job_id"
done

echo "Waiting for ${#job_ids[@]} jobs to complete..."
qxub monitor --wait-for-completion "${job_ids[@]}"
echo "Parallel analysis completed"
EOF

chmod +x scripts/run_parallel_analysis.sh
```

Use in DVC stages:

```yaml
stages:
  sample_analysis:
    cmd: ./scripts/run_parallel_analysis.sh src/analyze_sample.py "data/samples/*.csv" results/ sc
    deps:
      - scripts/run_parallel_analysis.sh
      - src/analyze_sample.py
      - data/samples/
    outs:
      - results/
```

## Key Takeaways

1. **Seamless integration**: qxub works naturally as DVC stage commands
2. **Environment flexibility**: Different stages can use different software environments
3. **Resource optimization**: Match resources to stage requirements
4. **Parallel processing**: Use `--terse` + `qxub monitor` for parallel stages
5. **Error handling**: Robust pipelines with proper error checking
6. **Reproducibility**: DVC + qxub ensures consistent, trackable workflows

## Next Steps

Now that you understand DVC integration:
- **[Remote Execution](11-remote-execution.md)** - Run DVC pipelines from your laptop
- Apply these patterns to your own data science workflows

DVC + qxub creates powerful, reproducible data science pipelines that scale from laptop development to HPC production.

---

**💡 Pro Tips:**
- Use `qxub exec --dry` during pipeline development to test resource allocations
- Create project-specific aliases for common stage patterns
- Use `qxub monitor --wait-for-completion` to make DVC stages block properly
- Include error handling in complex parallel stages
- Version control your `.qxub/config.yaml` with DVC project settings
