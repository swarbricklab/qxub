ARG PYTHON_VERSION=3.10
ARG QXUB_VERSION=main

FROM mambaorg/micromamba:1.5.5 AS builder

ARG PYTHON_VERSION
ARG QXUB_VERSION

COPY <<ABSCONDA_ENV /tmp/env.yaml
name: qxub
channels:
- conda-forge
dependencies:
- python=${PYTHON_VERSION}
- click
- tailer=0.4.1
- rich
- omegaconf
- pyyaml
- requests
- pip
- pip:
  - git+https://github.com/swarbricklab/qxub.git@${QXUB_VERSION}
ABSCONDA_ENV

RUN micromamba install -y -n base --channel conda-forge conda-pack && \
    micromamba create -y -n qxub -f /tmp/env.yaml && \
    micromamba run -n base conda-pack -p /opt/conda/envs/qxub -o /tmp/absconda-env.tar.gz && \
    mkdir -p /tmp/absconda-env && \
    tar -xzf /tmp/absconda-env.tar.gz -C /tmp/absconda-env && \
    rm /tmp/absconda-env.tar.gz

FROM debian:bookworm-slim AS runtime
COPY --from=builder /tmp/absconda-env/ /opt/conda/envs/qxub/
RUN if [ -x "/opt/conda/envs/qxub/bin/conda-unpack" ]; then \
        PATH="/opt/conda/envs/qxub/bin:${PATH}" /opt/conda/envs/qxub/bin/conda-unpack && \
        rm -rf /opt/conda/envs/qxub/conda-meta/history; \
    fi

ENV CONDA_DEFAULT_ENV=qxub
ENV CONDA_PREFIX=/opt/conda/envs/qxub
ENV PATH=/opt/conda/envs/qxub/bin:/opt/conda/bin:${PATH}

ENTRYPOINT ["qxub"]
