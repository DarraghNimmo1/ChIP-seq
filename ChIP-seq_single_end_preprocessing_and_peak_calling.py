####################################################
## Darragh Nimmo
## Trinity College Dublin
## April 2021
# Snakemake workflow for single-end ChIP-seq preprocessing and peak calling.
# Adjust for peak calling (e.g. broad/narrow) as necessary.
####################################################


###############################################################################################
###IMPORTS and variables
###############################################################################################

import os

import functools

configfile: 'config.yaml'

READS_DIR = config['reads_dir']

GENOME = config['genome']

BLACK_LIST = config['black_list'] 

PICARD_JAR = config['picard_jar']

CONTROL = config['control']


directory_function = functools.partial(os.path.join, config['results'])
BAM_DIR = directory_function('Bam')
PEAK_DIR = directory_function('Peak')



###############################################################################################
###Rules
###############################################################################################

rule all:
    input:
        expand(os.path.join(BIGWIG_DIR, '{sample}.bigwig'), sample = SAMPLE)
        expand(os.path.join(PEAK_DIR, '{sample}_processed_peaks.bed'), sample = SAMPLE)
        

#rule trim_adapters:
#    input:
#        reads = os.path.join(READS_DIR, '{sample}.fq.gz')
#    output:
#        reads = os.path.join(READS_DIR, '{sample}_trimmed.fq.gz')
#    params:
#        dir = READS_DIR
#        extra = "-a AGATCGGAAGAG"
#    message:
#        "Trimming adapters from the sequencing reads"
#    shell:
#        """
#        trim_galore -o {params.dir} {params.extra}  {input.reads}
#        """

rule alignment:
    input:
        reads = os.path.join(READS_DIR, '{sample}_trimmed.fq.gz')
    output:
        bam = os.path.join(BAM_DIR, '{sample}_mapped.bam')
    message:
        "Aligning reads to the reference genome"
    shell:
        """
        bowtie2 --very-sensitive -x /index/hg38/bowtie2/hg38.fa  -U {input.reads} --threads 48 | samtools view -@ 48 -bS - > {output}
        """

rule coordinate_sort_index_1:
    input:
        bam = rules.alignment.output.bam
    output:
        bam = os.path.join(BAM_DIR,  '{sample}_mapped.sorted.bam'),
        index = os.path.join(BAM_DIR, '{sample}_mapped.sorted.bam.bai')
    message:
        "Coordinate sorting and indexing the alignment file"
    shell:
        """
        samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index -@ 48 {output.bam}
        """
    
rule mark_duplicates:
    input:
        bam=rules.coordinate_sort_index_1.output.bam
    output:
        bam= os.path.join(BAM_DIR, SAMPLE+'_md.bam'),
        metrics= os.path.join(BAM_DIR, SAMPLE+"_md.txt")
    params:
	   picard = PICARD_JAR
    message:
        "Removing duplicates from the bam file."
    shell:
        """
        java -jar {params.picard} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES true
        """

rule coordinate_sort_index_2:
    input:
        bam = rules.mark_duplicates.output.bam
    output:
        bam = os.path.join(BAM_DIR,  '{sample}_md.sorted.bam'),
        index = os.path.join(BAM_DIR, '{sample}_md.sorted.bam.bai')
    message:
        "Coordinate sorting and indexing the alignment file"
    shell:
        """
        samtools sort -@ 48 -o {output.bam} {input.bam}; samtools index -@ 48 {output.bam}
        """

rule convert_bigwig_and_normalize:
    input:
        bam = rules.coordinate_sort_index_2.output.bam
    output:
        bigwig = os.path.join(BIGWIG_DIR, '{sample}.bigwig')
    params:
        extra = "-of bigwig -bs 10 --normalizeUsing CPM --effectiveGenomeSize 2913022398"
    shell:
        """
        bamCoverage {params.extra} -b {input.bam} -o {output.bigwig} -p 48
        """

rule call_peaks:
    input:
        bam = rules.coordinate_sort_index_2.output.bam
    output:
        BroadPeak = os.path.join(PEAK_DIR, '{sample}_peaks.broadPeak')
    params:
        dir = PEAK_DIR
        control = CONTROL
    messaage:
        "Peak calling."
    shell:
        """
        macs2 callpeak -t {input.bam} -f BAM -c {params.control} -g hs --outdir {params.dir} -n {wildcards.sample}
        """
        
rule preprocess_peaks:
	input:
		BroadPeak = rules.call_peaks.output.BroadPeak,
	output:
		bed = os.path.join(PEAK_DIR, '{sample}_processed_peaks.bed')
	params:
		BlackList = BLACK_LIST
	message:
		"Preprocessing peaks."
	shell:
		"""
		bedtools subtract -a {input.BroadPeak} -b {params.BlackList} > {output.bed}
		"""
