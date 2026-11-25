#!/bin/bash
# Run tests/cluster against the scheduler in this container (ci/README.md).
# Usage, as the `tester` user inside a test cluster container:
#   /qxub/ci/run-cluster-tests.sh slurm [pytest args...]
# The repository is mounted read-only at /qxub; it is copied, so the
# checkout on the host isn't written to.
set -euo pipefail
scheduler="$1"
shift

rm -rf "$HOME/qxub"
cp -r /qxub "$HOME/qxub"
if [ ! -x "$HOME/venv/bin/python" ]; then
    # QXUB_TEST_PYTHON: the image's Python 3.10+, if python3 is older
    "${QXUB_TEST_PYTHON:-python3}" -m venv "$HOME/venv"
fi
"$HOME/venv/bin/pip" install -q -e "$HOME/qxub" pytest

cd "$HOME/qxub"
QXUB_TEST_CLUSTER="$scheduler" exec "$HOME/venv/bin/python" -m pytest tests/cluster -ra "$@"
