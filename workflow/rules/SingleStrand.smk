rule SingleStrandCollapse:
    input:
        groupbyumi_bam = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.step2.pe_mapped_groupbyumi.bam",
        bwa_hg38_fasta = "resources/hg38/Homo_sapiens_assembly38.fasta"
    output:
        bam="results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.step6.r1_consensus_filt_mapped.bam",
    params:
        output_path = "results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/",
        tmpdir = "./tmpdir/"
    threads: 8
    shell:
        """
        bash workflow/scripts/FgbioSingleStrandPipeline.sh {input.groupbyumi_bam} {params.output_path} {wildcards.sample_id}.{wildcards.chrom} {input.bwa_hg38_fasta} {params.tmpdir}
        """        

rule AddSingleStrandBamTags:
    input:
        #original_bam = "/gpfs/commons/projects/high_depth_liquid_biopsy/hg38_crams_from_UG/5bp_trim/trimmed_bams/{sample_id}.{chrom}.bam",
        original_bam = "/gpfs/commons/projects/high_depth_liquid_biopsy/alcheng/UG_MANUSCRIPT/data/bams/ultima_duplex/{sample_id}/{sample_id}.{chrom}.bam",
        sscs_bam = "results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.step6.r1_consensus_filt_mapped.bam",
        groupbyumi_bam = "results/WgsDuplex/Duplex/{sample_id}/{chrom}/{sample_id}.{chrom}.step2.pe_mapped_groupbyumi.bam"
    output:
        tagged_bam = "results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.SE.tagged.bam"
    threads: 1
    resources: mem="80G"
    shell:
        """
        python workflow/scripts/GenerateSEBamWithSingleStrandTags.py \
            --sscs_bam {input.sscs_bam} \
            --groupbyumi_bam {input.groupbyumi_bam} \
            --output_bam {output.tagged_bam} \
            --original_bam {input.original_bam}
        """
    
rule SingleStrandFeatureMap:
    input:
        tagged_bam = "results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.SE.tagged.bam"
    output:
        vcf = "results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.featuremap.vcf"
    resources: mem="30G"
    shell:
        """
        bash workflow/scripts/SingleStrandFeatureMap.sh {input.tagged_bam} {output.vcf} {wildcards.chrom}
        """

rule AnnotateSingleStrandFeatureMap:
    input:
        FM = 'results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.featuremap.vcf',
        lofreq = "results/WgsDuplex/RawData/AlleleFractions/{sample_id}/{chrom}/{sample_id}.{chrom}.lofreq_unfiltered.vcf"
    output:
        tsv = 'results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.annotatedfeaturemap.sscsonly.tsv'
    threads: 1
    resources: mem="80G"
    shell:
        """
        python workflow/scripts/AnnotateSingleStrandFeatureMap.py --FeatureMap {input.FM} --lofreq_vcf {input.lofreq} --output {output.tsv} --sscs_only_output {output.tsv}
        """

rule ReformatSingleStrand:
    input:
        "results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.annotatedfeaturemap.sscsonly.tsv"
    output:
        'results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.sscsonly.annotated.reformated.bed'
    threads: 1
    params: gnomad = 'resources/hg38/blacklists/gnomad.genomes.v3.1.2.sites.{chrom}.AF0.0001.vcf.gz'
    resources: mem = "3G"
    shell:
        """
        tail -n +2 {input} | awk '{{$2=$2-1" "$2; print $0}}' | tr ' ' '\t' > {input}.tmp.bed
        bash workflow/scripts/blacklist_and_add_gnomad_for_bed.sh {input}.tmp.bed {input}.tmp2.bed {input}.tmp3.bed {params.gnomad}
        echo -e "chrom\tpos\tref\talt\taltID\tMI\tMI_pos\tRX\tX_CIGAR\tX_EDIST\tX_FC1\tX_FC2\tX_FILTERED_COUNT\tX_FLAGS\tX_LENGTH\tX_MAPQ\tX_READ_COUNT\tX_RN\tX_SCORE\trq\trs\tcD\tcE\tcU\tcS\tdR\tdS\tcC\tFM_passed\tnum_N\tsscs_length\tsscs_PIR\tsscs_PIR2\tsscs_bp\tDP\tAF\tSB\tref_top\tref_bot\talt_top\talt_bot\tgnomad_chrom\tgnomad_pos\tgnomad_id\tgnomad_ref\tgnomad_alt" > {output}
        cut -f1,3-47 {input}.tmp3.bed | awk '{{if ($42==".") print $0}}' | tr ' ' '\t' >> {output}
        rm {input}.tmp.bed {input}.tmp3.bed 
        """

rule DenoiseAndGenerateSingleStrandTrinucleotideCounts:
    input:
        'results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.sscsonly.annotated.reformated.bed'
    output:
        a = 'results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.singlestrand_denoised.trinucleotides.txt',
        b = 'results/WgsDuplex/SingleStrand/{sample_id}/{chrom}/{sample_id}.{chrom}.singlestrand_denoised.trinucleotides.pdf'
    conda: "wgsumiR"
    resources:
        mem = "80G"
    shell:
        """
        Rscript workflow/scripts/DenoiseAndGenerateSingleStrandTriN.R {input} {wildcards.sample_id} {output.a} {output.b}
        """

rule CollectAndCollapseSingleStrand:
    input:
        expand('results/WgsDuplex/SingleStrand/{{sample_id}}/{chrom}/{{sample_id}}.{chrom}.singlestrand_denoised.trinucleotides.txt', chrom = chromosomes)
    output:
        'results/WgsDuplex/SingleStrand/{sample_id}/{sample_id}.singlestrand_denoised.trinucleotides.txt'
    conda: "wgsumiR"
    shell:
        """
        Rscript workflow/scripts/AggregateChromLevelTrinucleotides.R {output} {wildcards.sample_id} {input}
        """