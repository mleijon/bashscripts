#!/bin/bash

########################################################################################################################
###############################################DEFAULT VALUES###########################################################
###SPADES######
nr_p=$(grep -c 'processor' /proc/cpuinfo)                  # Nr of threads set to all (used by spades)
mem_kb=$(grep MemTotal /proc/meminfo|awk '{print $2;}')
m=$((mem_kb*8/10485760))                                   # Memory limit for spades set to 80% of system RAM
flag=""                                                    # Spades operation mode defaults to standard
###DIAMOND####
ddb="/ssd2/classify/nr"                                    # diamond db path
                                                           # Not input parameter must be changed in th script
c=1
b=$m/6                                                     # Diamond block size (set to 1/6 of memory)
if ((b > 12))
then
b=12                                                       # but limited to 12
fi
##############
param_err=false
seq="c"                                                    # Sequence type set to contigs
usearch="/ssd2/classify/usearch11.0.667_i86linux64"
########################################################################################################################
while getopts ":m:p:b:c:f:s:r:h" opt; do
  case ${opt} in
    m)
      case $OPTARG in
        0.1|0.2|0.3|0.4|0.5|0.6|0.7|0.8|0.9|1.0)
          fact=$(dc<<<"$OPTARG 10 * p")
          fact=${fact%.*}
          m=$((mem_kb*fact/10485760))
          ;;
        *)
          echo "Invalid -m"
          param_err=true 
          ;;
      esac
      ;;
    p)
      if [ "$OPTARG" -eq "$OPTARG" ] 2> /dev/null && [ "$OPTARG" -le "$nr_p" ]; then
        nr_p=$OPTARG
      else
        echo "Invalid -p"
        param_err=true 
      fi
      ;;
    f)
      case $OPTARG in
        "isolate"|"sc"|"meta"|"plasmid"|"metaplasmid"|"metaviral"|"bio"|"rna"|"rnaviral"|"iontorrent")
          flag=$OPTARG
          ;;
        *)
          echo "Invalid -spades_flag"
          param_err=true
          ;;
      esac
      ;;
    b)
      if [[ $OPTARG =~ ^([0-9]+|([0-9]*[.][0-9]+))$ ]] && ! [[ $OPTARG =~  ^(0*(\.0+))$ ]]; then
        b=$OPTARG
      else
        echo "invalid -b"
        param_err=true
      fi 
      ;;
    c)
      case $OPTARG in 
        1|2|3|4)
          c=$OPTARG
          ;;
        *) 
          echo "Invalid -c"
          param_err=true
          ;;
      esac 
      ;;
    s)
      case $OPTARG in
        "r"|"c")
          seq=$OPTARG
          ;;
        *)
          echo "Invalid -s"
          param_err=true
          ;;
      esac
      ;;
    r)
    case $OPTARG in
      "n"|"y")
        remove=$OPTARG
        ;;
      *)
        echo "Invalid -r"
        param_err=true
        ;;
    esac
    ;;
    h)
      echo "script usage: ./$(basename "$0")"
      echo " [-s <Sequence type: reads (r) or spades contigs (c)>, default: c.]"
      echo " [-r <Remove intermediate files (yes(y)/no(n), default: y>"
      echo " [-m <fraction of memory used by spades (0.1-1.0 step 0.1)>, default: 0.8]"
      echo " [-p <nr of threads used by spades>, default: tall threads available]"
      echo " [-f <flag for spades operation mode (isolate, sc, meta, plasmid," 
      echo "     metaplasmid, metaviral, bio, rna, rnaviral, iontorrent)>, default: standard]"
      echo " [-b <diamond block size>, default: min(0.8*ram/6, 12)]"
      echo " [-c <diamond number of chunks (1, 2, 3 or 4)>, default: 1]"
      echo " [-h <Display this help text>]"
      exit 0
      ;;
    \?)
      echo "Invalid option: $OPTARG" 1>&2
      exit 1
      ;;
    :)
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      exit 1
      ;;
  esac
done
if $param_err; then
  exit 1
fi
########################################################################################################################
FILES="./fa/*_R1_*.fastq.gz"
##############################################TRIMMING WITH TRIMMOMATIC#################################################
for f in $FILES; do
  if [ ! -d ./trimmed ]; then
        mkdir ./trimmed
  fi
  out_dir="./trimmed"
  name=${f##*/fa/}
  name=${name%_L*}
  fi_out=$out_dir/$name".fastq.gz"
  wait;trimmomatic PE -threads "$nr_p" -phred33 -quiet -basein "$f" -baseout \
  "$fi_out" -summary "$out_dir"/'summary.log'  SLIDINGWINDOW:4:15 MINLEN:75
done
########################################################################################################################
FILES="./trimmed/*_1P*.fastq.gz"
################################################ASSEMBLY WITH SPADES####################################################
if [ "$seq" == "c" ]; then
  for f in $FILES; do
    name=${f##*/trimmed/}
    name=${name%_1P*}
    outdir=./"$flag"spades/$name
    if [ "$flag" == "meta" ]; then
      wait;spades.py -o "$outdir" "--""$flag" -1 "$f" -2 "${f/_1P/_2P}" -t "$nr_p" -m $m
    else
      if [ "$flag" == "" ]; then
        wait;spades.py -o "$outdir" -1 "$f" -2 "${f/_1P/_2P}" -s "${f/_1P/_1U}" -s "${f/_1P/_2U}" -t "$nr_p" -m $m
      else
        wait;spades.py -o "$outdir" "--""$flag" -1 "$f" -2 "${f/_1P/_2P}" -s "${f/_1P/_1U}" -s "${f/_1P/_2U}" -t "$nr_p" \
        -m $m
      fi
    fi
  done

########################################################################################################################
else
  if [ ! -d ./merged ]; then
        mkdir ./merged
  else
        rm -r ./merged
        mkdir ./merged
  fi
  FILES="./trimmed/*_1P*.fastq.gz"
  for f in $FILES; do
    name=${f##*/trimmed/}
    name=${name%_1P*}
    outdir=./merged
    { cat "$f"; cat "${f/_1P/_2P}"; cat "${f/_1P/_1U}"; cat "${f/_1P/_2U}"; } >> "$outdir"/"$name".fastq.gz
    wait;gunzip "$outdir"/"$name".fastq.gz
############################################DEREPLICATION WITH USEARCH##################################################
    wait;"$usearch" -fastx_uniques "$outdir"/"$name"".fastq" -fastaout "$outdir/$name"'_uq.fasta' -sizeout \
    -relabel Uniq -strand both
    wait; rm "$outdir"/"$name"'.fastq'
    wait;gzip "$outdir"/"$name""_uq.fasta"
  done
fi
########################################################################################################################
############################################CLASSIFICATION WITH DIAMOND#################################################
if [ "$seq" == "c" ]; then
  FILES="./""$flag""spades/*/contigs.fasta"
else
  FILES="./merged/*.fasta.gz"
fi
for f in $FILES; do
  if [ "$seq" == "c" ]; then
    outname=${f##*spades/}
    outname=${outname%/*}
  else
    outname=${f##*merged/}
    outname=${outname%.fasta.gz*}
  fi
  if [ ! -d ./diamond ]; then
    mkdir ./diamond
  fi
  wait;diamond blastx -d $ddb -q "$f" -o ./diamond/"$outname".daa --max-target-seqs 1 --evalue 1E-5 --outfmt 6 \
  qseqid full_qseq evalue staxids sscinames sskingdoms skingdoms sphylums -b "$b" -c "$c" 2>  \
  ./diamond/"$outname""_diamond.log";wait
########################################################################################################################
#################################################OUTPUT USING AWK#######################################################
  awk -v FS='\t' '{if ($6 == "Viruses") {gsub(/ /,"_");print ">"$1":"$5"\n"$2}}' > \
  ./diamond/"$outname""_viruses.fa" ./diamond/"$outname".daa
  awk -v FS='\t' '{if ($6 == "Eukaryota") {gsub(/ /,"_");print ">"$1":"$5"\n"$2}}' > \
  ./diamond/"$outname""_eukaryota.fa" ./diamond/"$outname".daa
  awk -v FS='\t' '{if ($6 == "Bacteria") {gsub(/ /,"_");print ">"$1":"$5"\n"$2}}' > \
  ./diamond/"$outname""_bacteria.fa" ./diamond/"$outname".daa
  awk -v FS='\t' '{if (!($6 == "Viruses") && !($6 =="Eukaryota") && !($6=="Bacteria"))
  {gsub(/ /,"_");print ">"$1":"$5"\n"$2}}' > ./diamond/"$outname""_other.fa" ./diamond/"$outname".daa
done
if [ "$seq" == 'r' ]; then
  FILES="./diamond/*.fa"
  for f in $FILES; do
    sed 's/;:/;/' "$f" > "$f".tmp
    mv "$f".tmp "$f"
  done
fi
if [ "$remove" == 'y' ]; then
  rm -r ./trimmed ./merged/ ./"$flag"spades 2> /dev/null
fi
gzip ./diamond/*.daa
