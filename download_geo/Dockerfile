# Pull miniconda3 image
FROM bioconductor/bioconductor_docker:devel

# Create /download_geo as directory and set it as work directory
WORKDIR /download_geo

# Install packages
RUN R -e "install.packages('BiocManager'); library('BiocManager'); BiocManager::install('FlowSorted.Blood.450k'); BiocManager::install('GEOquery'); BiocManager::install('FlowSorted.Blood.450k'); BiocManager::install('minfi'); BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg19'); BiocManager::install('IlluminaHumanMethylation450kmanifest'); BiocManager::install('limma'); BiocManager::install('missMethyl')"

# Correct the Blas error
RUN git clone https://github.com/bmbolstad/preprocessCore.git && cd preprocessCore/ && R CMD INSTALL --configure-args="--disable-threading" .

# Open bash
CMD ["/bin/bash"]
