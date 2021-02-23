# Download dataset
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109381/suppl//GSE109381_RAW.tar

# Build and run docker
docker build -t master-degree/download_geo .
docker container run --rm -it  -v $(pwd):/download_geo  master-degree/download_geo

# https://www.bioconductor.org/help/docker/
# docker run -it -e PASSWORD=bioc -v $(pwd):/download_geo bioconductor/bioconductor_docker:devel bash

