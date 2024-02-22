FROM ubuntu:20.04

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y python3 python3-pip git

RUN pip install --upgrade pip
RUN pip install --break-system-packages tensorflow cython==3.0.8 pandas numpy==1.22 anndata matplotlib
RUN pip install protobuf==3.20.*
COPY SNAF-0.7.0-py3-none-any.whl /
RUN pip install /SNAF-0.7.0-py3-none-any.whl

ENV NUMBA_CACHE_DIR=/tmp
