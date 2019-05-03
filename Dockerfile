FROM rocker/tidyverse
LABEL Maintainer="Xindi Guo <xindi.guo@sagebase.org>"

RUN apt-get install -y net-tools
RUN apt-get update -qq && apt-get -y install libffi-dev
RUN Rscript -e "install.packages('synapser', repos=c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install()"
RUN Rscript -e "BiocManager::install(c('monocle'))"
RUN Rscript -e "BiocManager::install(c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))"
RUN Rscript -e "devtools::install_github('cole-trapnell-lab/garnett')"

COPY bin/run-ic.R /usr/local/bin/
RUN chmod a+x /usr/local/bin/run-ic.R

COPY . garnett
WORKDIR garnett

RUN R CMD INSTALL .
