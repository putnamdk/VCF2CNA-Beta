#!/bin/bash

BASE_DIR=$PWD
tar -jxvf gb19.tar.bz2
tar -jxvf gb38.tar.bz2
cd gc
tar -jxvf GRCh37.gc.tar.bz2
tar -jxvf GRCh38.gc.tar.bz2
cd ../mapability
tar -jxvf hg19.map.tar.bz2
tar -jxvf hg38.map.tar.bz2
