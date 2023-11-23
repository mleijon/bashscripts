#!/bin/bash

tmpfile=$(mktemp)
sed '/^>/s/$/@@/' "$1"|sed -z 's/\n//g'|sed 's/@@/\n/g'|sed 's/>/\n>/g'|sed '/^$/d' > "$tmpfile"
nr_rows=$(sed '/^$/d' "$tmpfile"|wc -l)
nr_seqs=$((nr_rows/2))
if [ "$2" -gt $nr_seqs ]; then
	echo "Only $nr_seqs sequences."
	exit
fi
dseq=$(tail -n $((2*(nr_seqs-"$2"+1))) "$tmpfile"|head -2|tail -1)
dseq+=$(seq -s@ $((1+"$2"-${#dseq}%"$2"))|tr -d "[:digit:]")
echo "$dseq"|sed -r "s/(.{$3})/>split\n\1\\n/g"|tr -d "@"|\
awk '{sub("split", "split_" ++i); print; if (substr($1,1,1) != ">") --i}' -
