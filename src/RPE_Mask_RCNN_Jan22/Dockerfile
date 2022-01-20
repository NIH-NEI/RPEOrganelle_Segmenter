FROM nvidia/cuda:10.0-cudnn7-runtime-ubuntu18.04

RUN mkdir /app
ADD mrcnn /app/mrcnn/
ADD swigs /app/swigs/
COPY ./requirements.txt /app/requirements.txt

WORKDIR /app
ENV PYTHONPATH=/app:/app/packages

RUN apt-get update
RUN apt-get install ffmpeg libsm6 libxext6 -y
RUN apt-get -y install python3-pip

COPY ./*.py ./*.sh /app

RUN pip3 install --upgrade pip
RUN pip3 install --no-cache-dir -r requirements.txt -t /app/packages/

RUN apt-get -y install swig
RUN sh build_swig.sh

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64

RUN apt-get clean

VOLUME /rpemrcnn

ENTRYPOINT [ "python3" ]
