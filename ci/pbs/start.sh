#!/bin/bash
# Start the one-node PBS cluster, then wait (or run the given command).
set -euo pipefail
host=$(hostname -s)

sed -i "s/^PBS_SERVER=.*/PBS_SERVER=$host/" /etc/pbs.conf
mkdir -p /var/spool/pbs/mom_priv
cat > /var/spool/pbs/mom_priv/config <<CONFIG
\$clienthost $host
# Copy job output back with cp: one host, one filesystem
\$usecp *:/ /
# Check walltime every few seconds, not every two minutes
\$min_check_poll 5
\$max_check_poll 10
CONFIG

/etc/init.d/pbs start
export PATH=/opt/pbs/bin:$PATH
until qstat -B >/dev/null 2>&1; do sleep 1; done

qmgr -c "create node $host" 2>/dev/null || true
qmgr <<QMGR
set server job_history_enable = True
set server scheduler_iteration = 2
set server flatuid = True
create queue normal queue_type = execution
set queue normal enabled = True
set queue normal started = True
create queue express queue_type = execution
set queue express enabled = True
set queue express started = True
set server default_queue = normal
QMGR
until pbsnodes -a 2>/dev/null | grep -q "state = free"; do sleep 1; done
echo "PBS is up"

if [ $# -gt 0 ]; then exec "$@"; fi
exec sleep infinity
