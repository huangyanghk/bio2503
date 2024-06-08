#!/bin/bash
rm SRR_list.txt
for i in `seq 58 81`
do
        echo SRR21768${i} >> SRR_list.txt
done
prefetch --option-file SRR_list.txt
