FROM ubuntu:23.04

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y python3 python3-pip git

RUN rm /usr/lib/python3.11/EXTERNALLY-MANAGED
RUN pip install --upgrade pip
RUN pip install cython==3.0.8

RUN git clone https://github.com/spvensko/SNAF.git
WORKDIR SNAF
RUN python3 setup.py build_ext --inplace
