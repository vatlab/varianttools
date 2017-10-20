#
# Docker image for variant tools
#
FROM    continuumio/miniconda3

MAINTAINER Bo Peng <bpeng@mdanderson.org>

RUN     apt-get update
RUN     apt-get -y install swig  gcc g++ build-essential bzip2 libbz2-dev libz-dev curl git vim libblas-dev liblapack-dev

RUN	conda update python
RUN	conda install cython
RUN	pip install numpy scipy tables
WORKDIR /home/bpeng
RUN     git clone http://github.com/vatlab/VariantTools VariantTools
WORKDIR /home/bpeng/VariantTools
RUN     git pull
RUN     python setup.py install

ENV     HOME /home/bpeng
RUN     mkdir /home/bpeng/temp

# download hg19 reference genome and refGene database
WORKDIR /home/bpeng/temp
RUN     touch temp.vcf
RUN     vtools init test --build hg19
RUN     vtools import temp.vcf
RUN     vtools use refGene

WORKDIR /home/bpeng
RUN     rm -rf temp
