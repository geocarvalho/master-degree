# Download datasets

```
Rscript download_dataset.R  
```

# Arquivos

* [GSE72251 outputs](https://drive.google.com/file/d/1NcBOA6BtLO0kUWKibwNwbgbfUTbkLbJY/view?usp=sharing)

```
# zip
tar -czvf GSE72251_outputs.tar.gz GSE72251/
# unzip
tar -xvf GSE72251_outputs.tar.gz
```

---

# Atualizando o projeto

* O download dos dados foi feito de acordo com o notion, se necessario e possivel usar o PyMethyl para isso. Ele pareceu bem mais facil de lidar

## Usando pymethyl

```
docker pull joshualevy44/methylnet:0.1
docker pull joshualevy44/pymethylprocess:0.1.3
cd /home/labbe-x/Documentos/github/master-degree/download_geo/GSE51032
sudo docker container run --rm -it -v $(pwd):/pymethyl/gse/ joshualevy44/pymethylprocess:0.1.3
```

## Para utilizar o docker no servidor de Wilson

```
sudo docker container run --rm -it  -v $(pwd):/download_geo  download_geo
```

## Cell estimation feito pra homem e mulheres de mama

* Script `estimate_cell_count.R`, resultado:

GSE51032/GSE51032_cell_estimation.csv

> 660 linhas, com 659 pacientes

## Epismoker

https://github.com/geocarvalho/master-degree/tree/main/epismoker

* Resultado disponivel para os 659 pacientes

./GSE51032/epismoker_SSt_GSE109381.csv

## Criar classes automaticas com n automatico 

* O script create_classes.py ja faz isso pro dados fenotipicos
* Com esse arquivo vou criar o desenho da hipotese
* Funcionando vamos usar nos valores beta para criar classes e poder usar no SSDP+