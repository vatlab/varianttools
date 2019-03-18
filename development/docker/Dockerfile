#
# Docker image for variant tools
#
FROM    continuumio/miniconda3

MAINTAINER Bo Peng <bpeng@mdanderson.org>

RUN     apt-get update
RUN     apt-get -y install swig  gcc g++ build-essential bzip2 libbz2-dev libz-dev curl git vim libblas-dev liblapack-dev

RUN	conda update conda python
RUN	conda install cython numpy pandas scipy pytables
WORKDIR /home/bpeng
RUN     pip install variant_tools

#ENV     HOME /home/bpeng
#RUN     mkdir /home/bpeng/temp

# download hg19 reference genome and refGene database
#WORKDIR /home/bpeng/temp
#RUN     touch temp.vcf
#RUN     vtools init test --build hg19
#RUN     vtools import temp.vcf
#RUN     vtools use refGene

#WORKDIR /home/bpeng
#RUN     rm -rf temp

#RUN     mkdir /home/bpeng/temp

# download hg18 reference genome and refGene database
#WORKDIR /home/bpeng/temp
#RUN     touch temp.vcf
#RUN     vtools init test --build hg18
#RUN     vtools import temp.vcf
#RUN     vtools use refGene

#WORKDIR /home/bpeng
#RUN     rm -rf temp




