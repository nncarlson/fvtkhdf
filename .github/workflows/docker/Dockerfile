FROM ubuntu:24.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    gfortran \
    cmake \
    python3-pip \
    libhdf5-dev \
    libopenmpi-dev \
    libhdf5-openmpi-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install python packages globally
# Ubuntu 24.04 requires --break-system-packages for global pip install
RUN pip3 install --break-system-packages fypp sphinx furo

# Set environment variables for MPI
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

CMD ["/bin/bash"]
