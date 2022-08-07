#!/bin/bash

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/..

echo "Cleaning..."
find "$BASE" -regex '.*\.\(so\|d\|pcm\)$' | xargs rm -f
