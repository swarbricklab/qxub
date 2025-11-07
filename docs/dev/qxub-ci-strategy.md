# qxub CI runners: agent runbook (concise)

This is an internal runbook for operating qxub’s self-hosted GitHub Actions runners. It replaces container-per-job workflows with pre-configured runners to cut startup time, remove OpenSSL mismatch issues, and simplify auth.

## Objectives

- Primary: run CI on self-hosted runners labeled `qxub-runner` with a stable, pre-installed qxub.
- Secondary: allow opt-in dev overrides (`pip install -e .`) when testing features.
- Constraints: keep auth simple (WIF + SSH secrets), keep NCI environments in a known state, avoid multiple runner image variants.

## Runner model

- Single stable image: `ghcr.io/swarbricklab/qxub-runner:stable` (Ubuntu 20.04 for OpenSSL 1.1.1).
- Pre-installs: qxub (stable), gcloud SDK; ships SSH defaults and RNG setting.
- Usage modes:
  - Stable: use preinstalled qxub.
  - Dev override: checkout repo + `pip install -e .` as-needed.
- Label: `runs-on: qxub-runner`.

## NCI execution environments

- Stable (unchanged):
  - Source: `/g/data/a56/software/qsub_tools/`
  - Conda env: `qxub`
- Dev (managed by CI):
  - Source: `/scratch/a56/qxub-work/qxub/`
  - Conda env: `qxub-dev`
- Keep dev env current via a periodic/commit-driven job (core steps):
  - `git fetch && git checkout <sha>`
  - `conda activate qxub-dev`
  - `pip uninstall qxub -y && pip install -e .`

## Remote config sets (platform-driven)

Use two different config files that both reference the same platform name from the platform definition; vary only the remote activation/settings:

```yaml
# ~/.config/qxub/stable-testing.yaml
defaults: { platform: nci_gadi }
platforms:
  gadi:
    definition: file:///g/data/a56/config/xdg/qxub/platforms/nci_gadi.yaml
    remote:
      host: gadi
      working_dir: /scratch/a56/{user}/ci-work
      conda_init: |
        eval "$(conda shell.bash hook)"
        conda activate qxub

# ~/.config/qxub/dev-testing.yaml
defaults: { platform: nci_gadi }
platforms:
  gadi:
    definition: file:///g/data/a56/config/xdg/qxub/platforms/nci_gadi.yaml
    remote:
      host: gadi
      working_dir: /scratch/a56/{user}/ci-work
      conda_init: |
        eval "$(conda shell.bash hook)"
        conda activate qxub-dev
```

Note: Configs are flexible and reference platform definitions; platform names come from the definition and must not be altered per scenario.

## Auth and secrets (minimal)

- Workload Identity Pool: `hpci-workload-identity-pool`
- Service Account: `grunner-ci-sa@${{ secrets.PROJECT_ID }}.iam.gserviceaccount.com`
- Secrets: `atlas_service_account_private_key` (+ optional `atlas_service_account_pub_key`)
- Workflow steps (inline):
  - `google-github-actions/auth@v2` with `project_id`, `workload_identity_provider`, `service_account`.
  - Fetch SSH private key into `~/.ssh/atlas_key` and chmod 600.

## Canonical workflow shapes (for reference)

- Feature tests (stable triggers dev on NCI):
  - `qxub exec --config stable-testing.yaml --platform gadi --execdir /scratch/a56/qxub-work/qxub -- '<activate qxub-dev; pytest>'`
- Remote execution self-test (dev override on runner):
  - Checkout + `pip install -e .` then `qxub exec --config dev-testing.yaml --platform gadi -- hostname`
- Production usage in other repos:
  - `runs-on: qxub-runner` then `qxub exec --platform gadi -- <command>`

## Implementation checklist

1) Reuse HPCI infra (WIF, SA, SSH secret) — already provisioned.
2) Build/publish single stable runner image on `main` merges.
3) Ensure NCI `qxub-dev` exists and is writable under `/scratch/a56/qxub-work/qxub`.
4) Migrate workflows to `runs-on: qxub-runner`; inline the two auth steps.
5) Provide the two platform config files to operators (stable/dev).
6) Validate with a real `qxub exec` to `gadi` (queue `copyq` is fine for smoke).

## Operations

- Provisioning: attach `qxub-runner` label; confirm `qxub --version` and `gcloud --version` present.
- Updates (runner image): rebuild on `main`; roll runners gradually.
- Updates (NCI dev env): run the sync job to re-pin to commit and re-install `-e`.
- Monitoring: watch job latencies vs historical; regressions usually correlate to `pip install -e .` usage.
- Backup plan: if dev override breaks, fall back to stable (remove `pip install -e .`).

## Troubleshooting quick refs

- OpenSSL errors: ensure runner is 20.04-based; avoid installing system OpenSSL 3.x.
- Auth failures: re-check WIF provider string, SA email, and workload identity binding.
- SSH denied: verify key present at `~/.ssh/atlas_key` (600) and host alias for `gadi` resolves.
- Stale dev env: re-run sync steps (fetch/checkout, activate, reinstall editable).

## Success metrics

- >50% reduction in CI walltime for equivalent jobs.
- Zero OpenSSL-mismatch incidents across releases.
- No container rebuilds in routine workflows.
- Simple, copy-pastable auth steps across repos.
