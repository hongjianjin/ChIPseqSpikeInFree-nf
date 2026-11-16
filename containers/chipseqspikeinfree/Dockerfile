# Container for the ChIPseqSpikeInFree-nf pipeline: R + Bioconductor + ChIPseqSpikeInFree
#
# Built by .github/workflows/build-container.yml and pushed to
#   ghcr.io/hongjianjin/chipseqspikeinfree:<version>
#
# Manual build:
#   docker build --build-arg CSIF_REF=1.2.6 -t ghcr.io/hongjianjin/chipseqspikeinfree:1.2.6 .
#
# Singularity on the St. Jude HPC:
#   singularity pull chipseqspikeinfree_1.2.6.sif docker://ghcr.io/hongjianjin/chipseqspikeinfree:1.2.6

FROM bioconductor/bioconductor_docker:RELEASE_3_20

# git tag of stjude/ChIPseqSpikeInFree to install
ARG CSIF_REF=1.2.6

LABEL org.opencontainers.image.title="chipseqspikeinfree" \
      org.opencontainers.image.description="R environment with ChIPseqSpikeInFree for the ChIPseqSpikeInFree-nf Nextflow pipeline" \
      org.opencontainers.image.source="https://github.com/hongjianjin/ChIPseqSpikeInFree-nf" \
      org.opencontainers.image.licenses="GPL-3.0"

# Use only the container's R libraries. Without this, an R_LIBS/R_LIBS_USER exported on the host
# (e.g. a personal R library in $HOME) leaks in through Singularity and R loads packages
# compiled for another R version ("libRlapack.so: cannot open shared object file").
ENV R_LIBS=/usr/local/lib/R/site-library:/usr/local/lib/R/library \
    R_LIBS_USER=/dev/null \
    R_LIBS_SITE=/usr/local/lib/R/site-library

# procps provides `ps`, which Nextflow uses to collect task metrics
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'BiocManager::install(c( \
        "GenomicAlignments", "GenomicRanges", "IRanges", "Rsamtools" \
    ), ask = FALSE, update = FALSE)' \
    && Rscript -e 'install.packages(c("remotes", "optparse"), repos = "https://cloud.r-project.org")'

RUN Rscript -e "remotes::install_github('stjude/ChIPseqSpikeInFree@${CSIF_REF}', dependencies = c('Depends', 'Imports', 'LinkingTo'), upgrade = 'never')" \
    && Rscript -e 'library(ChIPseqSpikeInFree); stopifnot(file.exists(system.file("scripts/run_ChIPseqSpikeInFree.R", package = "ChIPseqSpikeInFree")))'

CMD ["R"]
