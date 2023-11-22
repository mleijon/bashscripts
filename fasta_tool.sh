#!/bin/bash

clear

Help()
{
   # Display Help
   echo
   echo "##########################################HELP###############################################"
   echo "# 'fasta_tools.sh' with the flag '-s' splits an input fasta file (arg. 1) into a user       #"
   echo "# defined number (arg. 2) of approx. equally sized files without breaking any fasta records #"
   echo "# between files. Optionally the prefix for the output files can be given (arg. 3), which    #"
   echo "# otherwise defaults to: 'split_fa'. With the flag '-e' a number of fasta sequence records  #"
   echo "# (arg 2) will be extracted from the input file (arg. 1). The arguments arg. 1-3 must be    #"
   echo "#  given after the flags (-s or -e). The maximum number of splits is 99.                    #"
   echo "#                                                                                           #"
   echo "# Syntax: fasta_split.sh [-s] [-e] [-h] arg1 arg2 [arg3]                                    #"
   echo "# options:                                                                                  #"
   echo "# -h     Print this Help.                                                                   #"
   echo "# -s     Splits fasta file (arg. 1) into a given (arg. 2) number of files. An output file   #"   
   echo "#        prefix (arg. 3) may optionally be given.                                           #"
   echo "# -e     Extract a number (arg. 2) of sequence records from the input fasta file (arg. 1)   #"
   echo "#############################################################################################"
   echo
   read -n 1 -s
}

split_file()
{
   if [ "$2" -gt  99 ]; then
      echo "Maximally 99 splits."; exit 1
   fi
   if [ $# -eq 2 ]; then
      prefix="split_fa"
   else
      prefix=$3
   fi
   for f in "$PWD/""$prefix"*;do
     rm -f "$f"
   done
   rec_cnt=$(grep -c '>' "$1")
   recs_file=$(awk -v a="$rec_cnt" -v b="$2" 'BEGIN{print int(a/b)}')
   split -dt'>' -a2 -l"$recs_file" --additional-suffix=.fasta "$1" "$prefix"
   files=$PWD/$prefix"*"
   for f in $files; do
      sed '/^>$/d' "$f"|
      awk '{ if ( substr($1,1,1) != ">" && NR == 1 ) print ">"$0; else print $0;}' - > "${f/$prefix/_$prefix}"
      mv "${f/$prefix/_$prefix}" "$f"
   done
   if [ "$2" -lt 10 ]; then
      cat "$f" >> "$(awk -v file="$f" -v nr="$2" 'BEGIN{sub(/.\.fasta/, nr - 1 ".fasta", file); print file}')"
   elif [ "$2" -eq 10 ]; then
      cat "$f" >> "$prefix""09.fasta"
   else
      cat "$f" >> "$(awk -v file="$f" -v nr="$2" 'BEGIN{sub(/..\.fasta/, nr - 1 ".fasta", file); print file}')"
   fi
   rm "$f"
   echo -e "Splitting completed!"
}

extract_seqs()
{
if [ "$2" -ge "$(grep -c '>' < "$1")" ]; then
   echo "Nr of sequence records in fasta file are less or equal to $2. No splitting."
   exit 0
fi
csplit -zs "$1" /"$(grep -m$(( $2 + 1 )) '>' "$1" |tail -n1)"/
mv xx00 "$(awk -v file="$1" -v nr="$2" 'BEGIN{sub(/\.fasta/, "_" nr ".fasta" , file); print file}')"
rm xx01
echo "Extraction completed!"
}

#############################################Input Control#######################################################
do_splits=false; do_extract=false
while getopts ":seh" opt; do
   case ${opt} in
      s)
         do_splits=true
         ;;
      e)
         do_extract=true
         ;;
      h)
         Help
         exit 0
         ;;
     \?)
         echo "Invalid option: $OPTARG" 1>&2
         Help
         exit 1
         ;;
   esac
done

shift $((OPTIND - 1))

re='^[0-9]+$'
if [ $# -lt 2 ]; then
   echo "At least two arguments must be provided"
   Help
   exit 1
elif [ ! -f "$1" ]; then
   echo "\"$1\" not found."
   Help
   exit 1 
elif ! [[ $2 =~ $re ]]; then
   echo "\"$2\" is not an integer."
   Help
   exit 1
fi

if [[ $do_splits == true && $do_extract == true ]]; then
   # shellcheck disable=SC2086
   split_file "$1" "$2" "$3"
   extract_seqs "$1" "$2"
elif [[ $do_splits ==  false && $do_extract == false ]]; then
   echo "Either of the flags -s and -e is required."
   Help
   exit 1
elif [[ $do_splits == true ]]; then
  split_file "$1" "$2"
else
  extract_seqs "$1" "$2"
fi
#################################################################################################################
