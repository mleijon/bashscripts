#!/usr/bin/env bash

clear
cd resources||exit

function download () {
    echo -e "Downloading: ${1##*/}"
    wget -q ftp://ftp.ncbi.nlm.nih.gov/"$1".md5
    wget ftp://ftp.ncbi.nlm.nih.gov/"$1"
    echo -ne "Done!\nVerifying download..."
    if md5sum --status -c "${1##*/}".md5; then
      echo -e "\n${1##*/} download OK!"
    else
      echo -e "\n${1##*/} download failed! Exits."
      exit $?
    fi
}

read -r filesize <<< "$(wget --spider ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz 2>&1 |\
grep "SIZE"|awk -v FS=' ' '{printf("%.0f", $5/2^30);}')"
echo -en "This script creates a diamond database based on GenBank nr partition.
The download size of nr.gz is $filesize GB. Continue (y/n)?"; read -r cnt
case $cnt in
  ([nN]) exit;;
  ([yY]) :;;
  (*) echo "Unrecognized entry. Exits."; exit;;
esac

if [ -f  nr.dmnd ]; then
  echo -e "Updating diamond (nr.dmnd) database.\n Removing old db..."
  rm -f nr.dmnd nodes.dmp names.dmp nr.gz taxdump.tar.gz
  rm -f GB_Release_Number prot.accession2taxid.FULL.gz
  echo -e "Done!\n"
else
  echo "Creating diamond (nr.dmnd) database."
fi

wget -q ftp://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number
echo "Current GenBank version is $(< GB_Release_Number)"
download pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
download pub/taxonomy/taxdump.tar.gz
download blast/db/FASTA/nr.gz
tar zxf taxdump.tar.gz nodes.dmp names.dmp
source "$CONDA_PREFIX"/etc/profile.d/conda.sh
conda list --name "svamv1" &>/dev/null
if [ ! "$?" -eq 0 ]; then
  echo "Installing svamv1 conda environment:"
	conda env create -f ../workflow/envs/env.yaml
  echo $'Done!\n'
fi
conda activate svamv1
diamond makedb --in nr.gz --db nr --taxonmap prot.accession2taxid.FULL.gz --taxonnodes nodes.dmp --taxonnames names.dmp
