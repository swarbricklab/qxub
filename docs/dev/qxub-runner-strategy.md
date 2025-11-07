# qxub runner strategy: operator notes (concise)

Aligns with the CI runbook. Operate pre-configured self-hosted runners for fast, reliable qxub execution using a single stable image and inline auth.

## Runner model

- Label: `qxub-runner`
- Image: `ghcr.io/swarbricklab/qxub-runner:stable` (Ubuntu 20.04; OpenSSL 1.1.1)
- Preinstalled: qxub (stable), gcloud SDK, SSH defaults and RNG setting
- Modes:
  - Stable (default): use preinstalled qxub
  - Dev override: checkout + `pip install -e .` only when needed

## Auth pattern (inline)

Use workload identity + SSH key directly in workflows (no custom setup action):

```yaml
- name: Authenticate to GCP
  uses: google-github-actions/auth@v2
  with:
    project_id: ${{ secrets.PROJECT_ID }}
    workload_identity_provider: ${{ secrets.WORKLOAD_IDENTITY_PROVIDER }}
    service_account: ${{ secrets.SERVICE_ACCOUNT }}
- name: Get SSH key
  run: |
    gcloud secrets versions access latest --secret="atlas_service_account_private_key" > ~/.ssh/atlas_key
    chmod 600 ~/.ssh/atlas_key
```

Infra reused: WIF pool `hpci-workload-identity-pool`, SA `grunner-ci-sa`, secrets `atlas_service_account_private_key` (+ optional pub key).

## NCI environments and configs

- Stable env: `/g/data/a56/software/qsub_tools/` → conda `qxub`
- Dev env: `/scratch/a56/qxub-work/qxub/` → conda `qxub-dev`
- Keep dev env in sync: `git fetch && git checkout <sha>`, `conda activate qxub-dev`, `pip uninstall qxub -y && pip install -e .`

Config vs platform definitions (critical): configs are flexible overlays referencing fixed platform definitions by their canonical name. Do not invent new platform names for scenarios. Use two config files (stable vs dev) that both reference the same platform (e.g., `gadi`) and change only `remote.conda_init`.

## Canonical workflow shapes

- Feature tests (stable triggers dev on NCI):
  - `qxub exec --config stable-testing.yaml --platform gadi --execdir /scratch/a56/qxub-work/qxub -- '<activate qxub-dev; pytest>'`
- Remote execution self-test (dev override on runner):
  - checkout + `pip install -e .`, then `qxub exec --config dev-testing.yaml --platform gadi -- hostname`
- Production usage (other repos):
  - `runs-on: qxub-runner` then `qxub exec --platform gadi -- <command>`

## Operator checklist

1) Build/publish stable runner image on `main` merges.
2) Deploy/label runners as `qxub-runner`; verify `qxub --version` and `gcloud --version`.
3) Provide/maintain two config files (stable/dev) referencing the same platform name.
4) Validate with a `copyq` smoke: `qxub exec --platform gadi -- hostname`.
5) Monitor latency; regressions often correlate with dev overrides (`pip install -e .`).
6) Roll image updates gradually; resync NCI dev env on relevant commits.

## Troubleshooting quick refs

- OpenSSL errors → ensure image is 20.04; avoid system OpenSSL 3.x
- Auth failures → re-check WIF provider, SA email, identity bindings
- SSH denied → key at `~/.ssh/atlas_key` (600), `gadi` host alias resolves
- Stale dev env → re-run sync and reinstall editable

## Notes

- Prefer inline auth over a bespoke setup action to reduce maintenance.
- All variability lives in configs; platform definition names are immutable and referenced verbatim.
