rule CramToBam:
    input:
        cram = "data/{sample_id}.cram"
    output:
        bam = "data/{sample_id}.bam",
        bai = "data/{sample_id}.bam.bai"
    threads:
        16
    shell:
        """
        samtools view -@ {threads} -h -O bam {input.cram} > {output.bam}
        samtools index {output.bam}
        """

rule SplitBamByChromosome:
    input:
        bam = "data/{sample_id}.bam"
    output:
        bam = "data/{sample_id}.{chrom}.bam"
    threads:
        16
    resources:
        mem= "1G"
    shell:
        """
        samtools view -@ {threads} -h -O bam {input.bam} {wildcards.chrom} > {output.bam}
        samtools index {output.bam}
        """

rule Flagstat:
    input:
        bam = "data/{sample_id}.{chrom}.bam"
    output:
        flagstat = "results/WgsDuplex/RawData/Flagstats/{sample_id}/{chrom}/{sample_id}.{chrom}.flagstat"
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        """
        
rule AlleleFrequency:
    input:
        bam = "data/{sample_id}.{chrom}.bam",
        hg38 = 'resources/hg38/Homo_sapiens_assembly38.fasta'
    output:
        vcf = "results/WgsDuplex/RawData/AlleleFractions/{sample_id}/{chrom}/{sample_id}.{chrom}.lofreq_unfiltered.vcf.gz",
        unzip = "results/WgsDuplex/RawData/AlleleFractions/{sample_id}/{chrom}/{sample_id}.{chrom}.lofreq_unfiltered.vcf"
    threads: 16
    resources: mem = "8G"
    conda: 'AlleleFrequency'
    shell:
        """
        lofreq call-parallel --pp-threads {threads} -f {input.hg38} --no-default-filter -A -B -a 1 -b 1 -o {output.vcf} {input.bam}
        gunzip -c {output.vcf} > {output.unzip}
        """