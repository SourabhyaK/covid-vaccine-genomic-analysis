# Comparative Genomic Analysis of COVID-19 Vaccine Platforms
# Uses the official Bioconductor image so R + Biostrings come preinstalled -
# no manual BiocManager wrestling like the bare-metal WSL setup required.

FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install the Bioconductor R packages this project needs
RUN R -e 'BiocManager::install(c("Biostrings", "msa"), update = FALSE, ask = FALSE)'

# Install Python 3 + pip (the base image is Debian-based)
RUN apt-get update && apt-get install -y python3 python3-pip && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /pipeline
COPY . /pipeline

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Default: run the test suite, then the full analysis (Python prototype + R)
CMD bash -c "\
    cd tests && pytest -v && cd .. && \
    echo '--- Python prototype ---' && \
    cd src && python3 build_proteins.py && python3 align.py && cd .. && \
    echo '--- R / Biostrings analysis ---' && \
    cd src && Rscript analysis.R && cd .. \
    "
