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
RUN     git clone http://github.com/vatlab/VariantTools  VariantTools

WORKDIR /home/bpeng/VariantTools
RUN     git fetch
RUN		git checkout v3
RUN     git pull
RUN     python setup.py install

#ENV     HOME /home/bpeng
#RUN     mkdir /home/bpeng/temp

# download hg19 reference genome and refGene database
# WORKDIR /home/bpeng/temp
#RUN     touch temp.vcf
#RUN     vtools init test --build hg19
#RUN     vtools import temp.vcf
#RUN     vtools use refGene

#WORKDIR /home/bpeng
#RUN     rm -rf temp

#RUN     mkdir /home/bpeng/temp

# download hg18 reference genome and refGene database
# WORKDIR /home/bpeng/temp
#RUN     touch temp.vcf
#RUN     vtools init test --build hg18
#RUN     vtools import temp.vcf
#RUN     vtools use refGene

#WORKDIR /home/bpeng
#RUN     rm -rf temp

# RUN     mkdir /home/bpeng/temp
# WORKDIR /home/bpeng/boost
# RUN		wget http://downloads.sourceforge.net/project/boost/boost/1.49.0/boost_1_49_0.tar.gz
# RUN		tar -P -xvf boost_1_49_0.tar.gz boost_1_49_0/boost boost_1_49_0/libs/iostreams boost_1_49_0/libs/regex boost_1_49_0/libs/filesystem boost_1_49_0/libs/detail # boost_1_49_0/libs/system
# RUN     rm -rf boost_1_49_0.tar.gz

# WORKDIR /home/bpeng
# RUN     rm -rf temp



