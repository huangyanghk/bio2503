#!/bin/bash
for i in `seq 58 81`
do
hisat2 -x ht2_files/HFTH1 -U SRR21763${i}/SRR21763${i}.out.fastq.gz -S  SRR21763${i}/SRR21763${i}.sam
done
