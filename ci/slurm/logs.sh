#!/bin/bash
# The Slurm daemons' recent log lines, for a failed CI run.
for f in /var/log/slurm/*.log; do echo "== $f"; tail -n 80 "$f"; done
