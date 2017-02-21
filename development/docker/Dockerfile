#
# Docker image for variant tools
#
FROM    mdabioinfo/simupop:latest

MAINTAINER Bo Peng <bpeng@mdanderson.org>

WORKDIR /home/bpeng
RUN     git clone http://github.com/BoPeng/VariantTools VariantTools
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
