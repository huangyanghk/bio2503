#!/bin/bash
GTF="/home/wsy/bio2503/project/upperstream/genomic.gtf"
for i in `seq 58 81`
do
        featureCounts -t exon -g gene_id -a $GTF -o SRR21763${i}/SRR21763${i}.all.txt SRR21763${i}/SRR21763${i}.sort.bam        cat SRR21763${i}/SRR21763${i}.all.txt | cut -f1,7 > SRR21763${i}/SRR21763${i}.count.txt
done
