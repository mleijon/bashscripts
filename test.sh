#!/usr/bin/env bash
dirs="/ngs/mikael/CLINF/2*/"
for d in  $dirs ; do
scp -r ${d}/ decypher@172.17.1.169:micke/data/clinf/reindeer/$(basename ${d})/;wait
done