#!/bin/bash
for i in `seq 58 81`
do
        fastp -i SRR21763${i}/SRR21763${i}.fastq.gz -o SRR21763${i}/SRR21763${i}.out.fastq.gz --json SRR21763${i}/SRR21763${i}.json --html SRR21763${i}/SRR21763${i}.html
done
