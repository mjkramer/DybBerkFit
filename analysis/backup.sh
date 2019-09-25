#!/bin/bash

tag=$1; shift

mkdir -p bak/$tag
mv *.pdf *.root pics bak/$tag
mkdir pics
