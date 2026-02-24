FROM ubuntu:24.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    build-essential \
    cmake \
    gfortran \
    git \
    libhdf5-dev \
    libhdf5-openmpi-dev \
    libopenmpi-dev \
    python3-pip \
    python3-setuptools \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python packages
# Ubuntu 24.04 requires --break-system-packages for pip install outside of a venv
RUN pip3 install --no-cache-dir --break-system-packages \
    fypp \
    sphinx \
    sphinx-rtd-theme

# Set environment variables for MPI to run as root (common in CI)
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Set default workdir
WORKDIR /workspace
