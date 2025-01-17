#!/bin/bash

GROUPUMIBAM=$1
OUTPATH=$2
BASENAME=$OUTPATH/$3
REF=$4
tmpdir=$5

                    
fgbio --tmp-dir $tmpdir CallMolecularConsensusReads -i $GROUPUMIBAM \
                                -o $BASENAME.step5.pe_consensus.bam \
                                -M 1 \
                                --threads 16 \
                                --consensus-call-overlapping-bases false 

fgbio --tmp-dir $tmpdir -Xmx25G FilterConsensusReads \
                                    -i $BASENAME.step5.pe_consensus.bam \
                                    -o $BASENAME.step5.pe_consensus_filt.bam \
                                    -M 2 \
                                    -r $REF \
                                    -N 0 \
                                    -E 1 -e 1 -n 1

samtools view -h -b -f 64 $BASENAME.step5.pe_consensus_filt.bam > $BASENAME.step5.r1_consensus_filt.bam

######## STEP 6
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
samtools view $BASENAME.step6.r1_consensus_filt_mapped.bam | wc -l > $BASENAME.step6.r1_consensus_filt_mapped.bam.sscs.counts.txt
rm $BASENAME.step5.pe_consensus.bam
rm $BASENAME.step5.pe_consensus_filt.bam
rm $BASENAME.step5.r1_consensus_filt.bam