#!/bin/bash
for i in `seq 58 81`
do
        samtools view -b SRR21763${i}/SRR21763${i}.sam > SRR21763${i}/SRR21763${i}.bam
        samtools sort SRR21763${i}/SRR21763${i}.bam -o SRR21763${i}/SRR21763${i}.sort.bam
done
