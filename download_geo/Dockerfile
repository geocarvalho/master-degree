# Pull miniconda3 image
FROM bioconductor/bioconductor_docker:devel

# Create /download_geo as directory and set it as work directory
WORKDIR /download_geo

# Install packages
RUN R -e "install.packages('BiocManager'); library('BiocManager'); BiocManager::install('GEOquery'); BiocManager::install('minfi'); BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg19'); BiocManager::install('IlluminaHumanMethylation450kmanifest')"

# Open bash
CMD ["/bin/bash"]