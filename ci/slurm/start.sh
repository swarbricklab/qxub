#!/bin/bash
# Start the one-node cluster, then wait (or run the given command).
set -euo pipefail

if [ ! -s /etc/munge/munge.key ]; then
    dd if=/dev/urandom bs=1 count=1024 of=/etc/munge/munge.key 2>/dev/null
fi
chown munge:munge /etc/munge/munge.key /run/munge
chmod 400 /etc/munge/munge.key
runuser -u munge -- munged

# The database is the `db` container (compose.yaml)
until slurmdbd 2>/dev/null && sleep 2 && sacctmgr -n list cluster >/dev/null 2>&1; do
    pkill slurmdbd 2>/dev/null || true
    sleep 2
done
sacctmgr -i add cluster qxub >/dev/null 2>&1 || true
sacctmgr -i add account testing Cluster=qxub >/dev/null 2>&1 || true
sacctmgr -i add user tester Account=testing >/dev/null 2>&1 || true

slurmctld
slurmd
until sinfo -h -o %T 2>/dev/null | grep -qE '^(idle|mixed|allocated)'; do
    scontrol update NodeName=slurm State=RESUME >/dev/null 2>&1 || true
    sleep 1
done
echo "Slurm is up"

if [ $# -gt 0 ]; then exec "$@"; fi
exec sleep infinity
