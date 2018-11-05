#! /bin/bash

dirname=${1:-mcguire}
scene=${2:-*}
format=${3:-obj}

echo ${dirname}/yocto-${format}/${scene}/*.${format}

for filename in ${dirname}/yocto-${format}/${scene}/*.${format}; do
    echo ../yocto-gl/bin/yitrace --double-sided ${filename} 
    ../yocto-gl/bin/yview --double-sided --eyelight ${filename}
done
