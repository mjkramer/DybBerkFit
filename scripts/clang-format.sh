#!/bin/bash

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..
cd $BASE

srcfiles=$(find ShapeFit toySpectra -regex '.*\.\(C\|cc\|cxx\|cpp\|h\|H\|hh\|hxx\|hpp\)$')
clang-format -i --assume-filename=foo.cpp --style=file $srcfiles
