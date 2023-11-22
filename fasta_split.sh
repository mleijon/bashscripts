#!/bin/bash

clear
#############################################Input Control#######################################################
Help()
{
   # Display Help
   echo
   echo "#############################################################################################"
   echo "# 'fasta_split.sh' splits an input fastfile (arg. 1) into a user defined number (arg. 2) of #"
   echo "# approx. equally sized files without splitting any fasta records between files. Optionally #"
   echo "# the prefix for the output files can be given (arg. 3), which otherwise defaults to:       #"
   echo "# 'split_fa.'                                                                               #"
   echo "#                                                                                           #"
   echo "# Syntax: fasta_split.sh arg1 arg2 [arg3]                                                   #"
   echo "# options:                                                                                  #"
   echo "# -h     Print this Help.                                                                   #"
   echo "#############################################################################################"
   echo
}
re='^[0-9]+$'
if [ $# -lt 2 ]; then
   echo "Filename and number of splits must be provided."
   Help
   exit; 1
elif [ ! -f "$1" ]; then
   echo "\"$1\" not found."; exit 1
elif ! [[ $2 =~ $re ]]; then
   echo "\"$2\" is not an integer."; exit 1
fi
if [ $# -eq 2 ]; then
   prefix="split_fa"
else
   prefix=$3
fi
for f in "$PWD"/"$prefix"*;do
  rm -f "$f"
done
#################################################################################################################
echo -n "Splitting \"$1\" into $2 files..."
rec_cnt=$(grep -c '>' "$1")
recs_file=$(awk -v a="$rec_cnt" -v b="$2" 'BEGIN{print int(a/b)}')
split -dt'>' -l"$recs_file" --additional-suffix=.fasta "$1" "$prefix"
clear
echo -n "Correcting output files..."
files=$PWD/$prefix"*"
for f in $files; do
   sed '/^>$/d' "$f"|
   awk '{ if ( substr($1,1,1) != ">" && NR == 1 ) print ">"$0; else print $0;}' - > "${f/$prefix/_$prefix}"
   mv "${f/$prefix/_$prefix}" "$f"
done
if [ "$2" -lt 10 ]; then
   cat "$f" >> "$(awk -v file="$f" -v nr="$2" 'BEGIN{sub(/.\.fasta/, nr - 1 ".fasta", file); print file}')"
   rm "$f"
elif [ "$2" -eq 10 ]; then
   cat "$f" >> "$prefix""09.fasta"
   rm "$f"
else
   : 
fi
clear
echo -e "\aDone!"

