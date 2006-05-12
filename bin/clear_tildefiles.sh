#!/bin/sh

echo "Removing files *~"
find . -name \*~ -print -exec rm \{\} \;
