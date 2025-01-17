#!/bin/bash

BAM=$1
OUTPATH=$2
BASENAME=$OUTPATH/$(basename $BAM .bam)
REF=$3
tmpdir=$4

######## STEP 1 CREATE UNMAPPED BAM THAT IS SIMULATED PE ##########

python workflow/scripts/CreateSyntheticPairedEndBam.py $BAM $BASENAME.step1.pe_mapped.bam

######## STEP 2 GROUP BY UMI
fgbio --tmp-dir $tmpdir GroupReadsByUmi -i $BASENAME.step1.pe_mapped.bam \
                        -o $BASENAME.step2.pe_mapped_groupbyumi.bam \
                        -m 0 \
                        -s paired \
                        -e 1 

######## STEP 3 COLLECT STATS AND CREATE CONSENSUS
fgbio --tmp-dir $tmpdir CollectDuplexSeqMetrics -i $BASENAME.step2.pe_mapped_groupbyumi.bam \
                                -o $BASENAME.step3.pe_mapped_groupbyumi.duplexseqmetrics.txt

fgbio --tmp-dir $tmpdir CallDuplexConsensusReads -i $BASENAME.step2.pe_mapped_groupbyumi.bam \
                                -o $BASENAME.step3.pe_consensus.bam \
                                --threads 16 \
                                --consensus-call-overlapping-bases false 

######## STEP 4 FILTER CONSENSUS
fgbio --tmp-dir $tmpdir -Xmx25G FilterConsensusReads -s true \
                                    -i $BASENAME.step3.pe_consensus.bam \
                                    -o $BASENAME.step4.pe_consensus_filt.bam \
                                    -M 2 1 1 \
                                    -r $REF \
                                    -N 0 \
                                    -E 1 -e 1 -n 1

######## STEP 5 GET R1 BACK
samtools view -h -b -f 64 $BASENAME.step4.pe_consensus_filt.bam > $BASENAME.step5.r1_consensus_filt.bam

######## STEP 6 ALIGN & INDEX
samtools fastq $BASENAME.step5.r1_consensus_filt.bam \
    | bwa mem -K 100000000 \
                -p \
                -v 3 \
                -t 8 \
                -Y $REF \
                /dev/stdin \
    | fgbio ZipperBams -i /dev/stdin \
                        -u $BASENAME.step5.r1_consensus_filt.bam \
                        -r $REF \
    | samtools sort -@8 -O bam \
    > $BASENAME.step6.r1_consensus_filt_mapped.bam
    
samtools index $BASENAME.step6.r1_consensus_filt_mapped.bam

rm $BASENAME.step1.pe_mapped.bam
rm $BASENAME.step3.pe_consensus.bam
rm $BASENAME.step4.pe_consensus_filt.bam
rm $BASENAME.step5.r1_consensus_filt.bam