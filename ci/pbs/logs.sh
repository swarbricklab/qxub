#!/bin/bash
# The PBS daemons' recent log lines, for a failed CI run.
for d in server_logs sched_logs mom_logs; do
    f=$(ls -t /var/spool/pbs/$d/* 2>/dev/null | head -1)
    [ -n "$f" ] && { echo "== $f"; tail -n 80 "$f"; }
done
