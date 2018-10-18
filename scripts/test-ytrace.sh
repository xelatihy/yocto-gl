#!/bin/bash
set -x #echo on
./bin/ytrace -o tests/run-tests/ytrace/output/$2 tests/textures/$1 && ./bin/yimdiff tests/run-tests/ytrace/output/$2 tests/run-tests/ytrace/result/$2
