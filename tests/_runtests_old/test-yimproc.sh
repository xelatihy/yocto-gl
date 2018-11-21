#!/bin/bash
set -x #echo on
./bin/yimproc -o $2 $1 && ./bin/yimdiff -o $4 $2 $3
