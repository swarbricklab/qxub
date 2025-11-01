# Dockerfile for qxub-compatible GitHub Actions runner
# Solves OpenSSL 3.x compatibility issues with hpci-scripts SSH keys
# Build version: 2024-11-01

# Use Ubuntu 20.04 which has OpenSSL 1.1.1 (compatible with hpci-scripts keys)
# This avoids the OpenSSL 3.x libcrypto compatibility issues
FROM ubuntu:20.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Core tools
    curl \
    wget \
    git \
    unzip \
    jq \
    # SSH client with OpenSSL 1.1.1 (compatible version)
    openssh-client \
    # Python ecosystem
    python3 \
    python3-pip \
    python3-venv \
    # Google Cloud SDK dependencies
    apt-transport-https \
    ca-certificates \
    gnupg \
    lsb-release \
    # Additional utilities
    file \
    && rm -rf /var/lib/apt/lists/*

# Verify OpenSSH/OpenSSL versions for compatibility debugging
RUN echo "=== Installed SSH/SSL versions ===" && \
    ssh -V && \
    openssl version && \
    echo "================================"

# Install Google Cloud SDK (needed for hpci-scripts)
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    rm -rf /var/lib/apt/lists/*

# Set up Python environment
RUN python3 -m pip install --upgrade pip setuptools wheel

# Create working directory
WORKDIR /workspace

# Pre-install common Python packages for faster CI
RUN pip install \
    pyyaml \
    click \
    requests \
    packaging \
    # Add other packages qxub commonly needs
    && pip cache purge

# Create symlink for python (some tools expect 'python' not 'python3')
RUN ln -sf /usr/bin/python3 /usr/bin/python

# Set up SSH directory
RUN mkdir -p /root/.ssh && chmod 700 /root/.ssh

# Set environment variables for optimal SSH behavior
ENV SSH_USE_STRONG_RNG=0

# Default working directory for GitHub Actions
WORKDIR /github/workspace

# Metadata
LABEL org.opencontainers.image.title="qxub GitHub Actions Runner"
LABEL org.opencontainers.image.description="Ubuntu 20.04 container with OpenSSL 1.1.1 for hpci-scripts SSH key compatibility"
LABEL org.opencontainers.image.version="1.0"
LABEL org.opencontainers.image.authors="Swarbrick Lab"
