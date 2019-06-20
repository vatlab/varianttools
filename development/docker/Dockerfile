#
# Docker image for variant tools
#
FROM    continuumio/miniconda3

MAINTAINER Bo Peng <bpeng@mdanderson.org>

RUN     apt-get update
RUN     apt-get -y install swig  gcc g++ build-essential bzip2 libbz2-dev libz-dev curl git vim libblas-dev liblapack-dev

RUN	conda update python
RUN	conda install cython

RUN wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.gz \
 && tar -zxvf hdf5-1.10.5.tar.gz \
 && cd hdf5-1.10.5 \
 && ./configure --prefix=/usr/local --enable-build-mode=production --enable-threadsafe --disable-hl \
 && make -j4 \
 && make install \
 && make clean

RUN wget http://download.zeromq.org/zeromq-4.0.3.tar.gz \
  && tar -zxvf zeromq-4.0.3.tar.gz 

RUN wget -O boost_1_49_0.tar.gz "http://downloads.sourceforge.net/project/boost/boost/1.49.0/boost_1_49_0.tar.gz?r=&ts=1435893980&use_mirror=iweb" \
 && tar -xf boost_1_49_0.tar.gz boost_1_49_0/boost boost_1_49_0/libs/iostreams boost_1_49_0/libs/regex boost_1_49_0/libs/filesystem boost_1_49_0/libs/detail boost_1_49_0/libs/system 

RUN conda install pytables scipy
WORKDIR /home/bpeng
RUN     git clone http://github.com/vatlab/VariantTools VariantTools
WORKDIR /home/bpeng/VariantTools
RUN     git pull
RUN mv /zeromq-4.0.3 ./src
RUN mv /boost_1_49_0 ./src

RUN     python setup.py install
# https://community.paperspace.com/t/storage-and-h5py-pytables-e-g-keras-save-weights-issues-heres-why-and-how-to-solve-it/430
# HDF5 locking issues
ENV     HDF5_USE_FILE_LOCKING FALSE


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

RUN     mkdir /home/bpeng/temp

# download hg18 reference genome and refGene database
WORKDIR /home/bpeng/temp
RUN     touch temp.vcf
RUN     vtools init test --build hg18
RUN     vtools import temp.vcf
RUN     vtools use refGene


WORKDIR /home/bpeng
RUN     rm -rf temp




