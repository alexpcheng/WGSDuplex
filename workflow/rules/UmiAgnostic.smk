rule AnnotateUmiAgnosticFeatureMap:
    input:
        FM = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.featuremap.vcf",
        lofreq = "results/WgsDuplex/RawData/AlleleFractions/{sample_id}/{chrom}/{sample_id}.{chrom}.lofreq_unfiltered.vcf"
    output:
        tsv = "results/WgsDuplex/UmiAgnostic/{sample_id}/{chrom}/{sample_id}.{chrom}.annotatedfeaturemap.tsv"
    threads: 1
    resources: mem = "199G"
    shell:
        """
        python workflow/scripts/AnnotateUmiAgnosticFeatureMap.py --FeatureMap {input.FM} --lofreq_vcf {input.lofreq} --output {output.tsv}
        """

rule ReformatUmiAgnostic:
    input:
        "results/WgsDuplex/UmiAgnostic/{sample_id}/{chrom}/{sample_id}.{chrom}.annotatedfeaturemap.tsv"
    output:
        'results/WgsDuplex/UmiAgnostic/{sample_id}/{chrom}/{sample_id}.{chrom}.umiagno.annotated.reformated.bed'
    threads: 1
    resources: mem = "3G"
    params: gnomad = 'resources/hg38/blacklists/gnomad.genomes.v3.1.2.sites.{chrom}.AF0.0001.vcf.gz'
    shell:
        """
        tail -n +2 {input} | awk '{{$2=$2-1" "$2; print $0}}' | tr ' ' '\t' > {input}.tmp.bed
        bash workflow/scripts/blacklist_and_add_gnomad_for_bed.sh {input}.tmp.bed {input}.tmp2.bed {input}.tmp3.bed {params.gnomad}
        echo -e "chrom\tpos\tref\talt\taltID\tRX\tX_CIGAR\tX_EDIST\tX_FC1\tX_FC2\tX_FILTERED_COUNT\tX_FLAGS\tX_LENGTH\tX_MAPQ\tX_READ_COUNT\tX_RN\tX_SCORE\trq\trs\tPIR\tPIR2\tDP\tAF\tSB\tref_top\tref_bot\talt_top\talt_bot\tgnomad_chrom\tgnomad_pos\tgnomad_id\tgnomad_ref\tgnomad_alt" > {output}
        cut -f1,3-34 {input}.tmp3.bed | awk '{{if ($29==".") print $0}}' | tr ' ' '\t' >> {output}
        rm {input}.tmp.bed {input}.tmp3.bed 
        """

rule DenoiseUmiAgnostic:
    input:
        'results/WgsDuplex/UmiAgnostic/{sample_id}/{chrom}/{sample_id}.{chrom}.umiagno.annotated.reformated.bed'
    output:
        "results/WgsDuplex/UmiAgnostic/{sample_id}/{chrom}/{sample_id}.{chrom}.umiagno.annotated.reformated.denoised.tsv"
    shell:
        """
        python workflow/scripts/DenoisingPipelineUmiAgnostic.py {input} {output}
        """

rule UmiAgnosticGenerateTriN:
    input:
        "results/WgsDuplex/UmiAgnostic/{sample_id}/{chrom}/{sample_id}.{chrom}.umiagno.annotated.reformated.denoised.tsv"
    output:
        a = 'results/WgsDuplex/UmiAgnostic/{sample_id}/{chrom}/{sample_id}.{chrom}.umiagno_denoised.trinucleotides.txt',
        b = 'results/WgsDuplex/UmiAgnostic/{sample_id}/{chrom}/{sample_id}.{chrom}.umiagno_denoised.trinucleotides.pdf'
    conda: "wgsumiR"
    resources:
        mem = "80G"
    shell:
        """
        Rscript workflow/scripts/GenerateUmiAgnosticTriN.R {input} {wildcards.sample_id} {output.a} {output.b}
        """

rule CollectAndCollapseUmiAgnostic:
    input:
        expand('results/WgsDuplex/UmiAgnostic/{{sample_id}}/{chrom}/{{sample_id}}.{chrom}.umiagno_denoised.trinucleotides.txt', chrom = chromosomes)
    output:
        'results/WgsDuplex/UmiAgnostic/{sample_id}/{sample_id}.umiagno_denoised.trinucleotides.txt'
    conda: "wgsumiR"
    shell:
        """
        Rscript workflow/scripts/AggregateChromLevelTrinucleotides.R {output} {wildcards.sample_id} {input}
        """

# rule BladderTumorFractions:
#     input:
#         bam = "/gpfs/commons/projects/high_depth_liquid_biopsy/alcheng/UG_MANUSCRIPT/data/bams/ultima_duplex/{sample_id}/{sample_id}.bam",
#         vcf = "results/Misc/BladderTumors/{sample_id}.somatic.snvs.clonal.vcf"
#     output:
#         md_bam = temp("results/{sample_id}.md.bam"),
#         md_bai = temp("results/{sample_id}.md.bam.bai"),
#         bam = temp("results/Misc/Bladder/{sample_id}/{sample_id}_int_{sample_id}/{sample_id}_int_{sample_id}.denoised.bam"),
#         featuremap = temp("results/Misc/Bladder/{sample_id}/{sample_id}_int_{sample_id}/{sample_id}_int_{sample_id}.denoised.featuremap.vcf"),
#         tf = "results/Misc/Bladder/{sample_id}/{sample_id}_int_{sample_id}/{sample_id}_int_{sample_id}.denoised.tumor_fraction.tsv"
#     resources: mem="90G"
#     params:
#         ref = "/gpfs/commons/projects/high_depth_liquid_biopsy/parameter_data/Homo_sapiens_assembly38.fasta"
#     threads: 8
#     shell:
#         """
#         samtools calmd -b -@ {threads} {input.bam} {params.ref} > {output.md_bam}
#         samtools index {output.md_bam}

#         python workflow/scripts/WgsNoDuplexUltimaTumorInformedTF.py \
#             --input_bam {output.md_bam} \
#             --output_bam {output.bam} \
#             --output_featuremap {output.featuremap} \
#             --tumor_vcf {input.vcf} \
#             --output_tf {output.tf} \
#             --sample_id {wildcards.sample_id}_int_{wildcards.sample_id}
#         """
