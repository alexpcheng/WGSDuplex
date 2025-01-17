human_chromosomes = [f'human_chr{str(i)}' for i in range(1,23)] + ['human_chrX', 'human_chrY'] 
mouse_chromosomes = [f'mouse_chr{str(i)}' for i in range(1,20)] + ['mouse_chrX', 'mouse_chrY']
PDX_samples = ['NYGC_DB_RN_PDX_duplexUDI', 'NYGC_DB_SS_PDX_duplexUDI', 'NYGC_RRR_PDX_duplexUDI', 'NYGC_ZA30_PDX_duplexUDI']

#human_chromosomes = ['human_chr21']
#mouse_chromosomes = ['mouse_chr18']
#PDX_samples = ['NYGC_ZA30_PDX_duplexUDI']

rule all:
    input:
        expand("results/PDX/UltimaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.snvaccuracy.txt", sample_id = PDX_samples, chrom = human_chromosomes+mouse_chromosomes),
        expand("results/PDX/IlluminaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.snvaccuracy.txt", sample_id = PDX_samples, chrom = human_chromosomes+mouse_chromosomes),
        
        expand("results/PDX/UltimaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.indelaccuracy.txt", sample_id = PDX_samples, chrom = human_chromosomes+mouse_chromosomes),
        expand("results/PDX/IlluminaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.indelaccuracy.txt", sample_id = PDX_samples, chrom = human_chromosomes+mouse_chromosomes),
        expand("results/UltimaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.indelaccuracy.fgbio.singlestrand.fam2.txt", sample_id = PDX_samples, chrom = human_chromosomes+mouse_chromosomes)

rule DuplexCollapsePDX:
    input:
        bam = "/gpfs/commons/projects/high_depth_liquid_biopsy/PDX_crams_from_UG/5bp_trim/trimmed_bams/{sample_id}.{chrom}.bam",
        bwa_hg38_index = "resources/hg38_mm39/hg38_mm39.fasta.sa",
        bwa_hg38_fasta = "resources/hg38_mm39/hg38_mm39.fasta"
    output:
        bam = "results/PDX/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.step6.r1_consensus_filt_mapped.bam",
        groupbyumi_bam = "results/PDX/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.step2.pe_mapped_groupbyumi.bam"
    params:
        output_path = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/",
        tmpdir = "./tmpdir/"
    threads: 8
    resources: mem="80G"
    shell:
        """
        bash workflow/scripts/FgbioDuplexPipeline.sh {input.bam} {params.output_path} {input.bwa_hg38_fasta} {params.tmpdir}
        """

rule AddDuplexBamTagsPDX:
    input:
        original_bam = "/gpfs/commons/projects/high_depth_liquid_biopsy/PDX_crams_from_UG/5bp_trim/trimmed_bams/{sample_id}.{chrom}.bam",
        duplex_bam = "results/PDX/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.step6.r1_consensus_filt_mapped.bam",
        groupbyumi_bam = "results/PDX/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.step2.pe_mapped_groupbyumi.bam"
    output:
        tagged_bam = "results/PDX/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.SE.tagged.bam"
    threads: 1
    resources: mem="200G"
    shell:
        """
        python workflow/scripts/GenerateSEBamWithDuplexTags.py \
            --dscs_bam {input.duplex_bam} \
            --groupbyumi_bam {input.groupbyumi_bam} \
            --output_bam {output.tagged_bam} \
            --original_bam {input.original_bam}
        """

rule SnvAccuracyPDX:
    input:
        md_bam = "results/PDX/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.SE.tagged.bam"
    output:
        "results/PDX/UltimaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.snvaccuracy.txt"
    resources: mem = "16G"
    shell:
        """
        python workflow/scripts/SnvErrorsPDX.py --bamfile {input.md_bam} --output {output}
        """

rule IndelAccuracyPDX:
    input:
        md_bam = "results/PDX/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.SE.tagged.bam"
    output:
        one = "results/PDX/UltimaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.indelaccuracy.txt"
    resources: mem = "64G"
    shell:
        """
        python workflow/scripts/IndelErrorsPDX.py --bamfile {input.md_bam} --output {output.one}
        """

rule SnvAccuracyIlluminaPDX:
    input:
        bam = "/gpfs/commons/projects/high_depth_liquid_biopsy/alcheng/UG_MANUSCRIPT/data/bams/illumina/PDX/{sample_id}.bam"
    output:
        bam = "results/PDX/IlluminaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.bam",
        txt="results/PDX/IlluminaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.snvaccuracy.txt"
    resources: mem = "16G"
    shell:
        """
        bash workflow/scripts/ReformatIlluminaToUltima.sh {input.bam} {output.bam} {wildcards.chrom}
        python workflow/scripts/SnvErrorsPDX.py --bamfile {output.bam} --output {output.txt}
        """

rule IndelAccuracyIlluminaPDX:
    input:
       bam = "results/PDX/IlluminaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.bam"
    output:
        one = "results/PDX/IlluminaBenchmarking/Indel/{sample_id}/{sample_id}.{chrom}.indelaccuracy.txt"
    resources: mem = "64G"
    shell:
        """
        python workflow/scripts/IndelErrorsPDX.py --bamfile {input.bam} --output {output.one}
        """