#!/bin/bash
set -x #echo on
./bin/ytrace -o $2 $4 $1 && ./bin/yimdiff $2 $3
