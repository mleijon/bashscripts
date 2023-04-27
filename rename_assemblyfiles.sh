#!/bin/bash

dir="/sfuSTORE1/micke/S.Aureus/SA1-mastit-SE/unicycler/*"
odir="/sfuSTORE1/micke/S.Aureus/SA1-mastit-SE/unicycler/unicycler-assembly-files"
for f in $dir; do
   if [ -d "$f" ]; then
      name=${f##*/}
      cp "$f"/assembly.fasta "$odir"/"$name".fasta
   fi
done
