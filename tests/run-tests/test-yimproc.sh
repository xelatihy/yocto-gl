#!/bin/bash
set -x #echo on
./bin/yimproc -o $2 $1 && ./bin/yimdiff $2 $3
