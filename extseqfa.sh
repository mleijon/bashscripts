#!/bin/bash

nr_rows=$(sed '/^$/d' $1|wc -l)
nr_seqs=$((nr_rows/2))
if [ $2 -gt $nr_seqs ]; then
	echo "Only $nr_seqs sequences."
	exit
fi
sed '/^$/d' $1|tail -n $((2*(nr_seqs-$2+1)))|head -2 


#sed -r "s/(.{5})/\1\\n/g"
