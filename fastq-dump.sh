#!/bin/bash
for i in `seq 58 81`
do
        fastq-dump --gzip SRR21763${i}/SRR21763${i}.sra --outdir SRR21763${i}
done
