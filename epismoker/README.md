# Install EpiSmokEr

```
> conda install -c conda-forge r-devtools
> R
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")
Sys.setenv("TAR" = "internal")
```

## It didn't work so I did it by hand 
```
> cp /tmp/RtmpjcrsLK/file1d74371525912/EpiSmokEr_0.1.0.tar.gz .
> R
install("/home/watson/george/master-degree/GSE109381/epismoker/EpiSmokEr")
```
