rule IndexHg38:
    input:
        hg38 = 'resources/hg38/Homo_sapiens_assembly38.fasta'
    output:
        fai = 'resources/hg38/Homo_sapiens_assembly38.fasta.fai',
        sa = 'resources/hg38/Homo_sapiens_assembly38.fasta.sa'
    shell:
        """
        samtools faidx {input.hg38}
        bwa index {input.hg38}
        """