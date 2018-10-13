# Yocto/GL scripts

This directory contains scripts quickly hacked up to simplify our workflow on
Yocto/GL. These scripts are mostly used on OsX and not debugged on other 
platforms. So you should consider them as just hacks.

All scripts should are intended to be executed from the project root directory.
For example, to build just run `./scripts/release.sh`.

- `build-release.sh`: build in release mode
- `build-xcode.sh`: run cmake and launch xcode on it
- `build-format.sh`: run clang-format on the repo
- `build-clean.sh`: remove the build and bin directories
- `_deprecated.py`: script not used often containing deprecated functionality 
    that has not been ported yet
