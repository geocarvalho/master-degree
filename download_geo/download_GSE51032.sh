#/bin/bash

# Download RAW data
mkdir GSE51032; cd GSE51032
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51032/suppl//GSE51032_RAW.tar

# Download matrix
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51032/matrix/GSE51032_series_matrix.txt.gz