# --config Option for Testing

The `--config` option allows you to override the default configuration hierarchy with a custom config file. This is particularly useful for:
- Testing from laptops/workstations with remote execution
- CI/CD pipelines with specific platform configurations
- Development and testing without modifying user configs

## Usage

```bash
qxub exec --config <path-to-config-file> [options] -- command
```

## Configuration Precedence

With `--config`, the precedence order becomes:
```
system < user < project < local < test < override
```

The override config file merges with (rather than replaces) lower precedence configs, allowing you to:
- Add new platforms while keeping existing ones
- Override specific settings (defaults, platforms, etc.)
- Test different configurations without affecting user config

## Example Config File

```yaml
# tests/configs/test_remote.yaml
defaults:
  platform: ci_test_remote
  project: a56

platforms:
  # Remote platform for CI/laptop testing
  ci_test_remote:
    name: ci_test_remote
    definition: file:///g/data/a56/config/xdg/qxub/platforms/nci_gadi.yaml
    remote:
      host: gadi  # SSH hostname from ~/.ssh/config
      working_dir: /scratch/a56/{user}
      conda_init: 'eval "$(conda shell.bash hook)"'

  # Local platform for direct PBS submission
  ci_test_local:
    name: ci_test_local
    definition: file:///g/data/a56/config/xdg/qxub/platforms/nci_gadi.yaml
```

## Usage Examples

### Local Execution with Override Config
```bash
# Use local platform (direct PBS submission)
qxub exec --config tests/configs/test_remote.yaml \
  --platform ci_test_local \
  --env pytorch \
  -- python train.py
```

### Remote Execution with Override Config
```bash
# Execute remotely via SSH (requires SSH access to platform)
qxub exec --config tests/configs/test_remote.yaml \
  --platform ci_test_remote \
  --env pytorch \
  -- python train.py
```

### Using Default Platform from Override Config
```bash
# Uses defaults.platform from override config
qxub exec --config tests/configs/test_remote.yaml \
  --env pytorch \
  -- python train.py
```

## CI/CD Integration

### GitHub Actions Example
```yaml
- name: Submit job to Gadi via qxub
  run: |
    qxub exec \
      --config .github/configs/ci_gadi.yaml \
      --platform gadi_remote \
      --env test_env \
      -- pytest tests/
```

### GitLab CI Example
```yaml
test_on_gadi:
  script:
    - qxub exec --config .gitlab/configs/gadi.yaml --platform gadi_remote -- pytest
```

## Testing

Run the test script to verify the `--config` option works:
```bash
./tests/test_config_option.sh
```

## Notes

- The override config file path is relative to the current directory
- Override config has highest precedence - overrides everything
- Platform definitions in override config merge with existing platforms
- Variable expansion works in override configs (`{user}`, `{project}`)
- SSH hostnames reference `~/.ssh/config` (not managed by qxub)
