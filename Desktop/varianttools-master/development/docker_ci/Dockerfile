#
# Docker image for variant tools
#
FROM    continuumio/miniconda3

MAINTAINER Bo Peng <bpeng@mdanderson.org>

RUN     apt-get update
RUN     apt-get -y install swig  gcc g++ build-essential bzip2 libbz2-dev libz-dev curl git vim libblas-dev liblapack-dev  libcurl4-openssl-dev libssl-dev


RUN	conda update python
RUN	pip install numpy scipy tables cython

RUN     conda install -c conda-forge hdf5 blosc gsl libboost
ENV 	LD_INCLUDE_PATH=/opt/conda/include/


WORKDIR /home/bpeng
RUN     git clone http://github.com/vatlab/VariantTools  VariantTools

WORKDIR /home/bpeng/VariantTools
RUN     git fetch
RUN		git checkout f74ee0c66e042f55d82c2a67d14c20e054e57597
RUN     python setup.py install

ENV     HOME /home/bpeng
RUN     mkdir /home/bpeng/temp

# download hg19 reference genome and refGene database
# WORKDIR /home/bpeng/temp
RUN     touch temp.vcf
RUN     vtools init test --build hg19
RUN     vtools import temp.vcf
RUN     vtools use refGene

WORKDIR /home/bpeng
RUN     rm -rf temp

RUN     mkdir /home/bpeng/temp

# download hg18 reference genome and refGene database
WORKDIR /home/bpeng/temp
RUN     touch temp.vcf
RUN     vtools init test --build hg18
RUN     vtools import temp.vcf
RUN     vtools use refGene

WORKDIR /home/bpeng
RUN     rm -rf temp
