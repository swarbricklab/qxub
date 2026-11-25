# Test clusters in containers

Each directory here is a one-node cluster for one scheduler: the controller,
one compute node and anything the scheduler needs (Slurm's accounting
database, PBS's PostgreSQL), in containers. `tests/cluster/` runs real jobs on
them with qxub. The main test suite stays on the fake commands in
`tests/fakes/`: it's fast and runs anywhere. These tests check what the fakes
can't, such as real command-line options and output, a job's whole life on a
compute node, and how qx learns that the scheduler stopped a job.

| Directory | Scheduler | Notes |
|---|---|---|
| `slurm/` | Slurm 24.11 (Debian 13) | `slurmdbd`'s database is a separate `mariadb:10.11` container, as slurmdbd 24.11 can't use Debian 13's MariaDB 11.8 |
| `pbs/` | OpenPBS 23.06, built from source (Rocky Linux 8) | Plain PBS, without Gadi's job-log epilogue (no `Exit Status:` line), like most PBS sites; job history is on, so `qstat -fx` knows finished jobs |

The `Cluster tests` workflow (`.github/workflows/cluster-tests.yml`) runs
both on every push and pull request. It builds each image from its
directory, reusing cached layers, so OpenPBS is compiled only when `ci/pbs/`
changes. Pushes to `main` whose tests pass publish the images as
`ghcr.io/swarbricklab/qxub-test-slurm` and `ghcr.io/swarbricklab/qxub-test-pbs`
(tagged `latest` and with the commit), along with the build cache later runs
start from.

## Running them yourself

From the repository root, with Docker running:

```bash
docker compose -f ci/slurm/compose.yaml up -d --build --wait
```

```bash
docker compose -f ci/slurm/compose.yaml exec -T -u tester slurm /qxub/ci/run-cluster-tests.sh slurm
```

```bash
docker compose -f ci/slurm/compose.yaml down
```

For PBS, replace `slurm` with `pbs` throughout. Building the PBS image
compiles OpenPBS, which takes several minutes the first time. To use a
published image instead (if the package is public, or after
`docker login ghcr.io`), pull and tag it, then start without `--build`:

```bash
docker pull ghcr.io/swarbricklab/qxub-test-pbs:latest
```

```bash
docker tag ghcr.io/swarbricklab/qxub-test-pbs:latest qxub-test-pbs
```

```bash
docker compose -f ci/pbs/compose.yaml up -d --no-build --wait
```

The repository is mounted read-only at `/qxub`. `run-cluster-tests.sh` copies
it to the `tester` user's home, installs it there in a virtual environment and
runs `pytest tests/cluster` with `QXUB_TEST_CLUSTER` set, so only that
scheduler's tests run. Elsewhere they are skipped. Pytest arguments go after
the scheduler's name (`… slurm -k walltime`).

`/usr/local/bin/scheduler-logs` in each container prints the daemons' recent
log lines.

## What they don't cover

- **Memory limits:** these aren't enforced, since the containers have no
  cgroups. Walltime and cancellation are, because the scheduler itself
  handles them.
- **Site setup:** Gadi's projects and storage, Pawsey's accounts and other
  site specifics need the real systems.
- **Memory use:** Slurm samples it every 5 seconds, so a job shorter than
  that shows none.
