# Platform Configuration

Platforms define HPC system capabilities and queue limits for intelligent queue selection.

## Auto Queue Selection

Use `--queue auto` to automatically select the best queue:

```bash
qxub --queue auto -l mem=2GB --env myenv -- python small_job.py
# → Selects 'express' or 'small' queue

qxub --queue auto -l mem=500GB --env myenv -- python big_job.py
# → Selects 'hugemem' queue

qxub --queue auto -l ngpus=1 -l ncpus=12 --env pytorch -- python train.py
# → Selects 'gpu' queue
```

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
