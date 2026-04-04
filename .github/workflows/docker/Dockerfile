FROM ubuntu:24.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive
ARG FPM_VERSION=0.13.0
ARG VTK_VERSION=9.6.1

RUN apt-get update && apt-get install -y \
    gfortran \
    cmake \
    python3-pip \
    libhdf5-dev \
    libopenmpi-dev \
    libhdf5-openmpi-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages globally. Ubuntu 24.04 requires
# --break-system-packages for global pip installs.
# Keep fypp in the image for both CMake and fpm workflows, and add VTK so the
# optional reader smoke tests run in CI.
RUN python3 -m pip install --break-system-packages \
    fypp \
    sphinx \
    furo \
    "vtk==${VTK_VERSION}"

# Install fpm in the image so the workflow does not need to download it on
# each run.
RUN python3 -c "import urllib.request; urllib.request.urlretrieve('https://github.com/fortran-lang/fpm/releases/download/v${FPM_VERSION}/fpm-${FPM_VERSION}-linux-x86_64-gcc-12', '/usr/local/bin/fpm')" \
    && chmod +x /usr/local/bin/fpm \
    && fpm --version

# Set environment variables for MPI
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

CMD ["/bin/bash"]
