FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
    bash \
    bowtie2 \
    samtools \
    bedtools
    
WORKDIR /app
