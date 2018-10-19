#!/bin/bash
set -x #echo on
./bin/ytrace -o $2 $4 $1 && ./bin/yimdiff -t 0.01 $2 $3
