#!/bin/bash


FILES="$PWD/*.txt"
for f in ${FILES}; do
  mv "${f}" "${f##*_}"
done

