FROM mambaorg/micromamba:1.5-focal

LABEL version="0.0.1" maintainer="Nicolas Vannieuwkerke <nicolas.vannieuwkerke@ugent.be>"

RUN micromamba install -y --name base -c conda-forge -c bioconda -c defaults \
    python==3.12 \
    tabix

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
