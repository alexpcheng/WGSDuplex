#!/bin/bash

module load java

BAM=$1
OUTBAM=$2
CHROM=$3

BASENAME=$(dirname $OUTBAM)/$(basename $OUTBAM .bam)
tmpdir=./tmpdir

samtools view -hb -f3 $BAM $CHROM > $BASENAME.step1.pe_mapped.bam

fgbio --tmp-dir $tmpdir GroupReadsByUmi -i $BASENAME.step1.pe_mapped.bam \
                        -o $BASENAME.step2.pe_mapped_groupbyumi.bam \
                        -m 0 \
                        -s paired \
                        -e 1

python workflow/scripts/GenerateSEBamWithUmiTagsPDXIllumina.py --input_bam $BASENAME.step2.pe_mapped_groupbyumi.bam --output_bam $OUTBAM

rm $BASENAME.step1.pe_mapped.bam
rm $BASENAME.step2.pe_mapped_groupbyumi.bam