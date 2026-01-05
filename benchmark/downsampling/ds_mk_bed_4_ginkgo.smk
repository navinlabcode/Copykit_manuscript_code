samples, = glob_wildcards("fulldata/{sample}.bam")
samtools_path="/volumes/seq/code/3rd_party/samtools/samtools-1.10/samtools"

rule all:
    input:
        expand('fulldata_bed/{sample}.bed', sample=samples),
        # expand('downsampled_sort_1M_bed/{sample}.bed', sample=samples),
        expand('downsampled_sort_750k_bed/{sample}.bed', sample=samples),
        expand('downsampled_sort_500k_bed/{sample}.bed', sample=samples),
        expand('downsampled_sort_250k_bed/{sample}.bed', sample=samples),
        expand('downsampled_sort_125k_bed/{sample}.bed', sample=samples),
        expand('downsampled_sort_75k_bed/{sample}.bed', sample=samples),
        expand('downsampled_sort_50k_bed/{sample}.bed', sample=samples),
        expand('downsampled_sort_25k_bed/{sample}.bed', sample=samples),
        expand('downsampled_sort_10k_bed/{sample}.bed', sample=samples),

rule downsampling_data:
    input:
        "fulldata/{sample}.bam"
    output:
        f1 = 'fulldata_bed/{sample}.bed',
        f2 = temp('fulldata/{sample}_rmdup.bam')
    shell:
      """
        /volumes/seq/code/3rd_party/git/sambamba-0.8.1/sambamba-0.8.1-linux-amd64-static markdup -r -t {threads} {input} {output.f2}
        bamToBed -i {output.f2} > {output.f1}
      """

# rule downsampling_1M:
#     input:
#         "downsampled_sort_1M/{sample}.bam"
#     output:
#         'downsampled_sort_1M_bed/{sample}.bed'
#     shell:
#         "bamToBed -i {input} > {output}"

rule downsampling_750k:
    input:
        "downsampled_sort_750k/{sample}.bam"
    output:
        'downsampled_sort_750k_bed/{sample}.bed'
    shell:
        "bamToBed -i {input} > {output}"

rule downsampling_500k:
    input:
        "downsampled_sort_500k/{sample}.bam"
    output:
        'downsampled_sort_500k_bed/{sample}.bed'
    shell:
        "bamToBed -i {input} > {output}"

rule downsampling_250k:
    input:
        "downsampled_sort_250k/{sample}.bam"
    output:
        'downsampled_sort_250k_bed/{sample}.bed'
    shell:
        "bamToBed -i {input} > {output}"

rule downsampling_125k:
    input:
        "downsampled_sort_125k/{sample}.bam"
    output:
        'downsampled_sort_125k_bed/{sample}.bed'
    shell:
        "bamToBed -i {input} > {output}"

rule downsampling_75k:
    input:
        "downsampled_sort_75k/{sample}.bam"
    output:
        'downsampled_sort_75k_bed/{sample}.bed'
    shell:
        "bamToBed -i {input} > {output}"

rule downsampling_50k:
    input:
        "downsampled_sort_50k/{sample}.bam"
    output:
        'downsampled_sort_50k_bed/{sample}.bed'
    shell:
        "bamToBed -i {input} > {output}"
        
rule downsampling_25k:
    input:
        "downsampled_sort_25k/{sample}.bam"
    output:
        'downsampled_sort_25k_bed/{sample}.bed'
    shell:
        "bamToBed -i {input} > {output}"

rule downsampling_10k:
    input:
        "downsampled_sort_10k/{sample}.bam"
    output:
        'downsampled_sort_10k_bed/{sample}.bed'
    shell:
        "bamToBed -i {input} > {output}"
