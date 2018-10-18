#!/bin/bash
set -x #echo on
./bin/yimproc -o tests/run-tests/yimproc/output/$2 tests/textures/$1 && ./bin/yimdiff tests/run-tests/yimproc/output/$2 tests/run-tests/yimproc/result/$2
