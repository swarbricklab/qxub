# Platform Configuration

Platforms define HPC system capabilities and queue limits for intelligent queue selection.

## Auto Queue Selection

Use `--queue auto` to automatically select the best queue:

````markdown
# Platform Configuration

Platforms define HPC system capabilities and queue limits for intelligent, cost-optimized queue selection.

## Auto Queue Selection

Use `--queue auto` to automatically select the most cost-effective queue for your job:

```bash
# Small jobs - uses normal queue
qxub --queue auto -l mem=8GB --env myenv -- python small_job.py
# → Selects 'normal' queue (2.0 SU/CPU·hr)

# High-memory jobs - cost-optimized selection
qxub --queue auto -l mem=300GB --env myenv -- python memory_job.py
# → Selects 'hugemembw' queue (1.25 SU/CPU·hr) - most cost-effective

qxub --queue auto -l mem=1200GB --env myenv -- python huge_job.py
# → Selects 'megamem' queue (1.25 SU/CPU·hr) - 58% cheaper than hugemem!

# GPU jobs
qxub --queue auto -l ngpus=1 -l ncpus=12 --env pytorch -- python train.py
# → Selects 'gpuvolta' queue

# Large-scale jobs
qxub --queue auto -l ncpus=5000 --env myenv -- python parallel_job.py
# → Selects 'normalsr' queue (high-core queue)
```

## Cost-Optimized Selection Rules

qxub automatically chooses the most cost-effective queue based on your requirements:

### Memory-Based Selection (Cost-Optimized)
- **≥1000GB memory**: `megamem` (1.25 SU/CPU·hr) - Best value for very high memory
- **≥500GB memory**: `hugemem` (3.0 SU/CPU·hr) - When megamem unavailable
- **≥192GB memory**: `hugemembw` (1.25 SU/CPU·hr) - Most cost-effective for high memory
- **<192GB memory**: Uses CPU-based selection

### CPU-Based Selection
- **>3200 CPUs**: `normalsr` (2.0 SU/CPU·hr) - Large-scale parallel
- **>1000 CPUs**: `normalbw` (1.25 SU/CPU·hr) - Cost-effective parallel
- **≤48 CPUs**: `normal` (2.0 SU/CPU·hr) - Default queue
- **≤3200 CPUs**: `normalsl` (1.5 SU/CPU·hr) - Balanced option

### Special Cases
- **GPU required**: `gpuvolta` (3.0 SU/CPU·hr) - GPU compute

## Cost Savings Example

For a 1200GB memory job:
- **Before**: Would select `hugemem` at 3.0 SU/CPU·hr
- **Now**: Selects `megamem` at 1.25 SU/CPU·hr
- **Savings**: 58% cost reduction!

## Queue Requirements (Validated)

Some queues have minimum requirements enforced by PBS:

- **hugemem**: Minimum 192GB memory per node
- **hugemembw**: Minimum 192GB memory per node
- **megamem**: Minimum 1000GB memory per node
- **gpuvolta**: Minimum 1 GPU

## All Available Queues

| Queue | Type | Billing Rate | Max Memory | Max CPUs | Use Case |
|-------|------|-------------|------------|----------|----------|
| `normal` | Standard | 2.0 | 190GB | 20,736 | General compute |
| `express` | Standard | 6.0 | 190GB | 3,168 | High priority |
| `normalsl` | Standard | 1.5 | 192GB | 3,200 | Balanced option |
| `normalbw` | Standard | 1.25 | 256GB | 10,080 | Cost-effective parallel |
| `expressbw` | Standard | 3.75 | 256GB | 1,848 | High priority parallel |
| `normalsr` | Standard | 2.0 | 500GB | 10,400 | Large-scale jobs |
| `expresssr` | Standard | 6.0 | 500GB | 2,080 | High priority large-scale |
| `hugemembw` | High-memory | 1.25 | 1,020GB | 140 | Cost-effective high memory |
| `hugemem` | High-memory | 3.0 | 1,470GB | 192 | High memory |
| `megamem` | Mega-memory | 1.25 | 2,990GB | 96 | Extreme memory (best value) |
| `gpuvolta` | GPU | 3.0 | 382GB | 960 | GPU compute |
| `copyq` | Data | 2.0 | 190GB | 1 | Data operations |

````

## Selection Rules

- **Memory > 100GB**: High-memory queues
- **GPU required**: GPU queues
- **Small jobs (≤2 CPUs, ≤2GB)**: Express/small queues
- **Default**: Normal queue

## Platform Files

Platform definitions are YAML files in:
- `/etc/xdg/qxub/platforms/`
- `~/.config/qxub/platforms/`
- `$QXUB_PLATFORM_PATHS`

## Custom Platforms

Set `QXUB_PLATFORM_PATHS` to use custom platform definitions:

```bash
export QXUB_PLATFORM_PATHS="/path/to/platforms:/other/path"
```

## Commands

```bash
qxub platform list              # List available platforms
qxub platform info              # Show current platform details
qxub select-queue --help        # Queue selection help
```
