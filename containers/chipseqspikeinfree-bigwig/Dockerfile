# Container for ChIPseqSpikeInFree-bigwig-nf: the ChIPseqSpikeInFree image (R, Bioconductor,
# ChIPseqSpikeInFree 1.2.6) plus Python with pybigtools, which reads bigWig files.
#
# Built by .github/workflows/build-container.yml and pushed to
#   ghcr.io/<owner>/chipseqspikeinfree-bigwig:<tag>
#
# Manual build:
#   docker build -t ghcr.io/hongjianjin/chipseqspikeinfree-bigwig:0.1.0 .
#
# Singularity on the St. Jude HPC:
#   singularity pull chipseqspikeinfree-bigwig_0.1.0.sif docker://ghcr.io/hongjianjin/chipseqspikeinfree-bigwig:0.1.0
#
# NOT BUILT OR TESTED YET: written on a host without Docker. Check with
#   docker run --rm <image> python3 -c "import pybigtools, numpy, requests; print(pybigtools.__file__)"

ARG BASE=ghcr.io/hongjianjin/chipseqspikeinfree:1.2.6
FROM ${BASE}

LABEL org.opencontainers.image.title="chipseqspikeinfree-bigwig" \
      org.opencontainers.image.description="R + ChIPseqSpikeInFree + Python + pybigtools for the ChIPseqSpikeInFree-bigwig-nf Nextflow pipeline" \
      org.opencontainers.image.source="https://github.com/hongjianjin/ChIPseqSpikeInFree-nf" \
      org.opencontainers.image.licenses="GPL-3.0"

# A venv keeps pip away from the system Python (works on Ubuntu releases that block system pip)
RUN apt-get update \
    && apt-get install -y --no-install-recommends python3 python3-venv \
    && rm -rf /var/lib/apt/lists/*

# pybigtools 0.3.0 is what bin/csif_bigwig_bin.py was tested with; numpy and requests are
# ranges (tested with numpy 2.4.6, requests 2.34.2, Python 3.11)
RUN python3 -m venv /opt/csif-py \
    && /opt/csif-py/bin/pip install --no-cache-dir "pybigtools==0.3.0" "numpy>=1.24,<3" "requests>=2.28"

ENV PATH=/opt/csif-py/bin:${PATH}

RUN python3 -c "import pybigtools, numpy, requests"

CMD ["R"]
