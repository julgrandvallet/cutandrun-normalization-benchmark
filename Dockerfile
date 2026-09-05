FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="cutandrun-normalization-benchmark"
LABEL org.opencontainers.image.source="https://github.com/julgrandvallet/cutandrun-normalization-benchmark"
LABEL org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml \
 && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH

COPY --chown=$MAMBA_USER:$MAMBA_USER bin/ /opt/pipeline/bin/
USER root
RUN chmod +x /opt/pipeline/bin/*
USER $MAMBA_USER
ENV PATH=/opt/pipeline/bin:$PATH

# Fail the build, not a six-hour pipeline run, if the analysis code is broken.
RUN scale_factors.py --demo && compare_methods.py --demo \
 && Rscript -e 'library(edgeR); library(DESeq2); library(optparse); cat("R deps OK\n")'
