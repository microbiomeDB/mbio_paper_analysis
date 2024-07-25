FROM rocker/rstudio:4.3.2

USER root
RUN apt-get update && apt-get install -y \
    libglpk-dev \
	libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    libgsl-dev 


# The following from Rserve dockerfile
### Rserve
RUN R -e "install.packages('Rserve', version='1.8-9', repos='http://rforge.net')"

### CRAN
RUN R -e "install.packages('bit64')"
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('jsonlite')"
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('Rcpp')"
RUN R -e "install.packages('readr')"
RUN R -e "install.packages('rlist')"

### Bioconductor
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('SummarizedExperiment')"
RUN R -e "BiocManager::install('DESeq2')"

## Profiling and development help packages
RUN R -e "print('installing devtools2')"
RUN R -e "install.packages('devtools')"

RUN R -e "install.packages('vegan')"

RUN R -e "remotes::install_github('zdk123/SpiecEasi','v1.1.1', upgrade_dependencies=F)" 
RUN R -e "BiocManager::install('Maaslin2')"

## Our packages
# RUN R -e "remotes::install_github('VEuPathDB/veupathUtils', 'v2.5.5', upgrade_dependencies=F)"
# RUN R -e "remotes::install_github('VEuPathDB/plot.data', 'v5.3.1', upgrade_dependencies=F)"
RUN R -e "remotes::install_github('microbiomeDB/microbiomeComputations', 'v5.1.3', upgrade_dependencies=F)"
RUN R -e "remotes::install_github('microbiomeDB/MicrobiomeDB')"
RUN R -e "remotes::install_github('microbiomeDB/microbiomeData')"

COPY . /home/rstudio/Documents

EXPOSE 8787

CMD ["/init"]