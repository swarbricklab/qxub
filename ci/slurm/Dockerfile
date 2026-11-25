# A one-node Slurm cluster for qxub's cluster tests (tests/cluster).
#
# slurmctld, slurmd, slurmdbd (so sacct works) and munge all run in this
# one container; its database is the `db` container (compose.yaml), as
# slurmdbd 24.11 can't use Debian 13's MariaDB 11.8. Jobs run as the
# `tester` user, who also runs qxub, so qxub's status files are on the same
# filesystem as the jobs. See ci/README.md.
FROM debian:trixie

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        slurm-wlm slurmdbd munge \
        python3 python3-venv python3-pip git ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*

RUN useradd -m -s /bin/bash tester \
    && mkdir -p /etc/slurm /var/spool/slurmctld /var/spool/slurmd /var/log/slurm /run/munge \
    && chown slurm:slurm /var/spool/slurmctld /var/log/slurm

COPY slurm.conf slurmdbd.conf cgroup.conf /etc/slurm/
COPY start.sh /usr/local/bin/start-slurm
COPY logs.sh /usr/local/bin/scheduler-logs
RUN chown slurm:slurm /etc/slurm/slurmdbd.conf && chmod 600 /etc/slurm/slurmdbd.conf \
    && chmod +x /usr/local/bin/start-slurm /usr/local/bin/scheduler-logs

CMD ["/usr/local/bin/start-slurm"]
