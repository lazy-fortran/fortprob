#!/bin/sh
set -eu
# FortSym is a development dependency only.
fortsym_dir=${FORTSYM_DIR:-$HOME/code/lazy-fortran/fortsym}
expected=063008070086aaf2d78f6e9d9a8fddff44c74a46
test "$(git -C "$fortsym_dir" rev-parse HEAD)" = "$expected"
mkdir -p build/codegen
gfortran -I"$fortsym_dir/build/include" gen/generate.f90 \
    -L"$fortsym_dir/build/lib" -Wl,-rpath,"$fortsym_dir/build/lib" \
    -lfortsym -lfortgen -o build/codegen/generate
build/codegen/generate
