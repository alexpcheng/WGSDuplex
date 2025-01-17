rule DuplexCollapse:
    input:
       bam = "data/{sample_id}.{chrom}.bam",
       bwa_hg38_index = "resources/hg38/Homo_sapiens_assembly38.fasta.sa",
       bwa_hg38_fasta = "resources/hg38/Homo_sapiens_assembly38.fasta"
    output:
        bam = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.step6.r1_consensus_filt_mapped.bam",
        groupbyumi_bam = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.step2.pe_mapped_groupbyumi.bam"
    params:
        output_path = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/",
        tmpdir = "./tmpdir/"
    threads: 8
    resources: mem="80G"
    shell:
        """
        bash workflow/scripts/FgbioDuplexPipeline.sh {input.bam} {params.output_path} {input.bwa_hg38_fasta} {params.tmpdir}
        """

rule AddDuplexBamTags:
    input:
        original_bam = "data/{sample_id}.{chrom}.bam",
        duplex_bam = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.step6.r1_consensus_filt_mapped.bam",
        groupbyumi_bam = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.step2.pe_mapped_groupbyumi.bam"
    output:
        tagged_bam = temp("results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.SE.tagged.bam")
    threads: 1
    resources: mem="160G"
    shell:
        """
        python workflow/scripts/GenerateSEBamWithDuplexTags.py \
            --dscs_bam {input.duplex_bam} \
            --groupbyumi_bam {input.groupbyumi_bam} \
            --output_bam {output.tagged_bam} \
            --original_bam {input.original_bam}
        """
        
rule Calmd:
    input:
        tagged_bam = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.SE.tagged.bam"
    output:
        md_bam = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.SE.mdtagged.bam"
    threads: 16
    params:
        ref = "resources/hg38/Homo_sapiens_assembly38.fasta"
    resources: mem="4G"
    shell:
        """
        samtools calmd -b -@ {threads} {input.tagged_bam} {params.ref} > {output.md_bam}
        samtools index {output.md_bam}
        """

rule DuplexFeatureMap:
    input:
        tagged_bam = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.SE.mdtagged.bam"
    output:
        vcf = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.featuremap.vcf"
    resources: mem="40G"
    shell:
        """
        bash workflow/scripts/DuplexFeatureMap.sh {input.tagged_bam} {output.vcf} {wildcards.chrom}
        """

rule AnnotateDuplexFeatureMap:
    input:
        FM = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.featuremap.vcf",
        lofreq = "results/WgsDuplex/RawData/AlleleFractions/{sample_id}/{chrom}/{sample_id}.{chrom}.lofreq_unfiltered.vcf"
    output:
        tsv = 'results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.annotatedfeaturemap.tsv',
        dup = 'results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.annotatedfeaturemap.duplexonly.tsv'
    threads: 1
    resources: mem="80G"
    shell:
        """
        python workflow/scripts/AnnotateDuplexFeatureMap.py --FeatureMap {input.FM} --lofreq_vcf {input.lofreq} --output {output.tsv} --duplex_only_output {output.dup}
        """

rule DuplexCollapseML:
    input:
        FM= 'results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.annotatedfeaturemap.tsv',
        model= 'resources/duplex_collapsing_model05152023.pkl'
    output:
        'results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.annotatedfeaturemap.mlduplex.tsv'
    resources: mem="64G"
    shell:
        """
        python workflow/scripts/DuplexCollapseML.py --AnnotatedFeatureMap {input.FM} --output {output} --model {input.model}
        """

rule ReannotateMLDuplex:
    input:
        expand('results/WgsDuplex/Duplex/{{sample_id}}/{chrom}/{{sample_id}}.{chrom}.annotatedfeaturemap.mlduplex.tsv', chrom = chromosomes)
    output:
        'results/WgsDuplex/Duplex/{sample_id}/{sample_id}.mlduplex.annotated.allchr.tsv'
    threads: 1
    params: path = 'results/WgsDuplex/Duplex/{sample_id}'
    resources: mem = "4G"
    shell:
        """
        python workflow/scripts/ReannotateMLDuplex.py --sample_id {wildcards.sample_id} --path {params.path} --output {output}
        """

rule Reformat:
    input:
        'results/WgsDuplex/Duplex/{sample_id}/{sample_id}.mlduplex.annotated.allchr.tsv'
    output:
        'results/WgsDuplex/Duplex/{sample_id}/{sample_id}.mlduplex.annotated.allchr.reformated.bed'
    threads: 1
    params: path = 'results/WgsDuplex/Duplex/{sample_id}',
            gnomad = 'resources/hg38/blacklists/gnomad.genomes.v3.1.2.sites.chr21.AF0.0001.vcf.gz'
    resources: mem = "4G"
    shell:
        """
        tail -n +2 {input} | awk '{{$2=$2-1" "$2; print $0}}' | tr ' ' '\t' > {input}.tmp.bed
        bash workflow/scripts/blacklist_and_add_gnomad_for_bed.sh {input}.tmp.bed {input}.tmp2.bed {input}.tmp3.bed {params.gnomad}
        echo -e "chrom\tpos\tref\talt\taltID\tMI\tMI_pos\tRX\tX_CIGAR\tX_EDIST\tX_FC1\tX_FC2\tX_FILTERED_COUNT\tX_FLAGS\tX_LENGTH\tX_MAPQ\tX_READ_COUNT\tX_RN\tX_SCORE\trq\trs\taD\taE\tbD\tbE\tcD\tcE\tcU\tcS\tdR\tdS\tcC\tFM_passed\tnum_N\tduplex_length\tduplex_PIR\tduplex_PIR2\tduplex_bp\tDP\tAF\tSB\tref_top\tref_bot\talt_top\talt_bot\tA_top\tA_bot\tT_top\tT_bot\tC_top\tC_bot\tG_top\tG_bot\tN_top\tN_bot\tduplex_bp2\tpir\tpir2\ttriN\tgnomad_chrom\tgnomad_pos\tgnomad_id\tgnomad_ref\tgnomad_alt" > {output}
        cut -f1,3-65 {input}.tmp3.bed | awk '{{if ($60==".") print $0}}' | tr ' ' '\t' >> {output}
        rm {input}.tmp.bed {input}.tmp3.bed
        """

rule DenoiseAndGenerateTrinucleotideCounts:
    input:
        'results/WgsDuplex/Duplex/{sample_id}/{sample_id}.mlduplex.annotated.allchr.reformated.bed'
    output:
        a = 'results/WgsDuplex/Duplex/{sample_id}/{sample_id}.duplex_denoised.trinucleotides.txt',
        b = 'results/WgsDuplex/Duplex/{sample_id}/{sample_id}.duplex_denoised.trinucleotides.pdf'
    conda: "wgsumiR"
    resources:
        mem = "32G"
    shell:
        """
        Rscript workflow/scripts/DenoiseAndGenerateTriN.R {input} {wildcards.sample_id} {output.a} {output.b}
        """