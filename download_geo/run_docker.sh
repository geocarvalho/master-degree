docker build -t master-degree/download_geo .
docker container run --rm -it  -v $(pwd):/download_geo  master-degree/download_geo

# https://www.bioconductor.org/help/docker/
# docker run -it -e PASSWORD=bioc -v $(pwd):/download_geo bioconductor/bioconductor_docker:devel bash