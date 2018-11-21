#!/bin/bash
set -x #echo on
mkdir -p $2
mkdir -p $2/meshes
mkdir -p $2/models
mkdir -p $2/textures
./bin/yscnproc -o $3 $1 && ./bin/ytrace -o $4 $7 $3 && ./bin/yimdiff -o $6 -t 0.01 $4 $5
