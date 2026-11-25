# A one-node OpenPBS cluster for qxub's cluster tests (tests/cluster).
#
# The server, scheduler, MoM (execution host) and PBS's own PostgreSQL all
# run in this one container, built from OpenPBS's source. Jobs run as the
# `tester` user, who also runs qxub. Unlike Gadi's PBS, nothing appends an
# `Exit Status:` line to the job log. See ci/README.md.
FROM rockylinux:8
LABEL org.opencontainers.image.source=https://github.com/swarbricklab/qxub
LABEL org.opencontainers.image.description="qxub's one-node OpenPBS test cluster"

ARG OPENPBS_VERSION=23.06.06

RUN dnf install -y dnf-plugins-core epel-release \
    && dnf config-manager --set-enabled powertools \
    && dnf install -y \
        gcc gcc-c++ make libtool autoconf automake rpm-build \
        hwloc-devel libX11-devel libXt-devel libedit-devel libical-devel \
        ncurses-devel perl postgresql-devel postgresql-contrib python3-devel \
        tcl-devel tk-devel swig expat-devel openssl-devel libXext libXft \
        cjson-devel \
        expat libedit postgresql-server python3 sendmail sudo tcl tk \
        libical cjson hostname procps-ng which git \
    && dnf clean all

RUN curl -sSL "https://github.com/openpbs/openpbs/archive/refs/tags/v${OPENPBS_VERSION}.tar.gz" \
        | tar xz -C /tmp \
    && cd "/tmp/openpbs-${OPENPBS_VERSION}" \
    && ./autogen.sh >/dev/null \
    && ./configure --prefix=/opt/pbs >/dev/null \
    && make -j"$(nproc)" >/dev/null \
    && make install >/dev/null \
    && cd / && rm -rf "/tmp/openpbs-${OPENPBS_VERSION}" \
    && /opt/pbs/libexec/pbs_postinstall \
    && chmod 4755 /opt/pbs/sbin/pbs_iff /opt/pbs/sbin/pbs_rcp \
    && sed -i 's/^PBS_START_MOM=0/PBS_START_MOM=1/' /etc/pbs.conf

# qxub needs Python 3.10 or later; Rocky 8's python3 is 3.6
RUN dnf install -y python3.11 python3.11-pip && dnf clean all
ENV QXUB_TEST_PYTHON=python3.11

RUN useradd -m -s /bin/bash tester
COPY start.sh /usr/local/bin/start-pbs
COPY logs.sh /usr/local/bin/scheduler-logs
RUN chmod +x /usr/local/bin/start-pbs /usr/local/bin/scheduler-logs
ENV PATH=/opt/pbs/bin:$PATH

CMD ["/usr/local/bin/start-pbs"]
