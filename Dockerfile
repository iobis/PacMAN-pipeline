FROM continuumio/miniconda3

WORKDIR /code

RUN apt-get update && apt-get install bsdmainutils
RUN conda install -n base -c conda-forge mamba
RUN mamba create -c conda-forge -c bioconda -n snakemake snakemake
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN echo "source activate snakemake" > ~/.bashrc
