################################################
#ChIP-seq processing script
#Darragh Nimmo
#November 2021
#################################################

import os
import functools

configfile: 'config.yaml'

READS_DIR = config['reads']

SAMPLE = config['sample']

INPUT_BAM = config['input_bam']

EXCLUSION_HG19 = config['exclusion_list']

directory_function = functools.partial(os.path.join, config['results'])
BAM_DIR = directory_function('Bam')
BIGWIG_DIR = directory_function('BigWig')
PEAK_DIR = directory_function('Peak')

#ls *fq.gz |while read line ; do bbduk.sh in="$line" out=trimmed/"$line" ref=/usr/local/bin/resources/adapters.fa k=23 mink=11 hdist=1 ktrim=r int=f ; done

rule all:
        input:
                expand(os.path.join(PEAK_DIR, '{sample}_peaks.preprocessed.bed'), sample = SAMPLE),
                expand(os.path.join( BIGWIG_DIR, '{sample}_marked_dups.sorted.bw'), sample = SAMPLE)


rule align_to_genome:
    input:
        read =  os.path.join(READS_DIR, 'MDAMB468_{sample}.fq.gz')
    output:
        bam = temp(os.path.join(BAM_DIR, '{sample}_mapped.bam'))
    shell:
        """
        bowtie2 --very-sensitive -x /index/hg19/bowtie2/hg19  -U {input.read} --threads 24 | samtools view -@ 24 -bS - > {output.bam}
        """
rule coordinate_sort_index_1:
    input:
        bam = rules.align_to_genome.output.bam
    output:
        bam = temp(os.path.join(BAM_DIR, '{sample}_mapped.sorted.bam' )),
        txt = temp(os.path.join(BAM_DIR, '{sample}_mapped.sorted.bam.bai'))
    shell:
        """
        samtools sort -@ 24 -o {output.bam} {input.bam}; samtools index -@ 24 {output.bam}
        """

rule remove_duplicates:
    input:
        bam = rules.coordinate_sort_index_1.output.bam
    output:
        bam = temp(os.path.join(BAM_DIR, '{sample}_marked_dups.bam')),
        txt = temp(os.path.join(BAM_DIR, '{sample}_marked_dups.txt'))
    shell:
        """
        java -jar /home/darragh/picard.jar MarkDuplicates -I {input.bam} -O {output.bam} -M {output.txt} --REMOVE_DUPLICATES true
        """

rule coordinate_sort_index_2:
    input:
        bam = rules.remove_duplicates.output.bam
    output:
        bam = os.path.join(BAM_DIR, '{sample}_marked_dups.sorted.bam' ),
        txt = os.path.join(BAM_DIR, '{sample}_marked_dups.sorted.bam.bai' )
    shell:
        """
        samtools sort -@ 24 -o {output.bam} {input.bam}; samtools index -@ 24 {output.bam}
        """
rule make_bigwigs:
    input:
        bam = rules.coordinate_sort_index_2.output.bam
    output:
        bigwig = os.path.join( BIGWIG_DIR, '{sample}_marked_dups.sorted.bw')
    shell:
        """
        bamCoverage -of bigwig -p 24 --normalizeUsing CPM -b {input.bam} -o {output.bigwig}
        """

rule call_peaks:
    input:
        bam = rules.coordinate_sort_index_2.output.bam,
        input = INPUT_BAM
    output:
        peak = os.path.join(PEAK_DIR, '{sample}_peaks.narrowPeak')
    params:
        name = SAMPLE,
        dir = PEAK_DIR
    shell:
        """
        macs2 callpeak -t {input.bam} -c {input.input} -q 0.05 -n {params.name} -f BAM -g hs --outdir {params.dir}
        """
#--broad --broad-cutoff 0.05
#Set Q value

rule preprocess_peaks:
        input:
                peak = rules.call_peaks.output.peak
        output:
                peak = os.path.join(PEAK_DIR, '{sample}_peaks.preprocessed.bed')
        params:
                exclude = EXCLUSION_HG19
        shell:
                """
                bedtools intersect -v -a {input.peak} -b {params.exclude} | grep -P 'chr[\dXY]+[ \t]' > {output.peak}
                """
