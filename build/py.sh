#!/bin/sh
# Run python from the nTop_pkgs env with its DLL dirs on PATH.
E=/c/Users/Benji/anaconda3/envs/nTop_pkgs
export PATH="$E:$E/Library/mingw-w64/bin:$E/Library/usr/bin:$E/Library/bin:$E/Scripts:$E/bin:$PATH"
exec "$E/python.exe" "$@"
