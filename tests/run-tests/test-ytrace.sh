#!/bin/bash
set -x #echo on
./bin/ytrace -o $2 $5 $1 && ./bin/yimdiff -o $4 -t 0.01 $2 $3
