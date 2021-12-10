
################################################
#ChIP-seq processing script
#Darragh Nimmo
#November 2021
#################################################

import os
import glob

file_dir = "/home/darragh/ChIP-seq_2/data/concat/PDAC"

file_list = []

cell_label = "HPEC"
#assay_label = "H3K27me3"

configfile: 'config.yml'


#CELL = config['cell_type']
ASSAY = config['assay_type']

rule all:
    input:
        expand(os.path.join(file_dir,cell_label,'{assay}', cell_label+'_{assay}_fastqc.html'), assay = ASSAY),
        expand(os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}.sorted.bam'), assay = ASSAY),
        expand(os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}_maked_dup.bam'), assay = ASSAY),
        expand(os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}_maked_dup.sorted.bam'), assay = ASSAY),
        expand(os.path.join(file_dir,cell_label,'{assay}','Bigwig','{assay}_maked_dup.sorted.bw'), assay = ASSAY)




rule concatenate_lanes:
    run:
        for assay_label in assay_list:
            if not os.path.exists('/home/darragh/ChIP-seq_2/data/concat/'+cell_label+'/'+assay_label):
                os.makedirs('/home/darragh/ChIP-seq_2/data/concat/'+cell_label+'/'+assay_label)#

                os.system('cp /home/darragh/ChIP-seq_2/data/'+cell_label+'*'+assay_label+'*/*.fastq.gz' ' /home/darragh/ChIP-seq_2/data/concat/'+cell_label+'/'+ assay_label)#

            os.system('zcat /home/darragh/ChIP-seq_2/data/concat/'+cell_label+'/'+assay_label+'/*.fastq.gz | pigz >  /home/darragh/ChIP-seq_2/data/concat/'+cell_label+'/'+assay_label+'/'+cell_label+'_'+assay_label+'.fq.gz')
            os.system('rm /home/darragh/ChIP-seq_2/data/concat/'+cell_label+'/'+assay_label+'/*.fastq.gz')


rule check_quality:
    input:
        raw_read = os.path.join(file_dir,cell_label,'{assay}', cell_label+'_R2_{assay}.fq.gz')
    output:
        html_file = os.path.join(file_dir,cell_label,'{assay}', cell_label+'_{assay}_fastqc.html')
    shell:
        """
        fastqc -t 48 {input.raw_read}
        """

rule grouped_quality:
    run:
        if not os.path.exists(os.path.join(file_dir,cell_label,"Quality")):
            os.makedirs(os.path.join(file_dir,cell_label,"Quality"))
        for assay in assay_list:
            os.replace(os.path.join(file_dir,cell_label,assay,cell_label+'_'+assay+'_fastqc.html'), os.path.join(file_dir,cell_label,"Quality",cell_label+'_'+assay+'_fastqc.html'))
            os.replace(os.path.join(file_dir,cell_label,assay,cell_label+'_'+assay+'_fastqc.zip'), os.path.join(file_dir,cell_label,"Quality",cell_label+'_'+assay+'_fastqc.zip'))
        multiqc .


rule align_to_genome:
    input:
        read = os.path.join(file_dir,cell_label,'{assay}',cell_label+'_R2_{assay}.fq.gz')
    output:
        bam = os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}_mapped.bam')
    shell:
        """
        bowtie2 --very-sensitive -x /index/hg38/bowtie2/hg38.fa  -U {input.read} --threads 6 | samtools view -@ 6 -bS - > {output.bam}
        """
rule coordinate_sort_index_1:
    input:
        bam = rules.align_to_genome.output.bam
    output:
        bam = os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}.sorted.bam'),
        txt = os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}.sorted.bam.bai')
    shell:
        """
        samtools sort -@ 6 -o {output.bam} {input.bam}; samtools index -@ 6 {output.bam}
        """

rule remove_duplicates:
    input:
        bam = os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}.sorted.bam')
    output:
        bam = os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}_maked_dup.bam'),
        txt = os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}_maked_dup.txt')
    shell:
        """
        java -jar /home/darragh/bin/picard.jar MarkDuplicates -I {input.bam} -O {output.bam} -M {output.txt} --REMOVE_DUPLICATES true
        """

rule coordinate_sort_index_2:
    input:
        bam = os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}_maked_dup.bam')
    output:
        bam = os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}_maked_dup.sorted.bam'),
        txt = os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}_maked_dup.sorted.bam.bai')
    shell:
        """
        samtools sort -@ 6 -o {output.bam} {input.bam}; samtools index -@ 6 {output.bam}
        """

rule make_bigwigs:
    input:
        bam = os.path.join(file_dir,cell_label,'{assay}','Bam','{assay}_maked_dup.sorted.bam')
    output:
        bigwig = os.path.join(file_dir,cell_label,'{assay}','Bigwig','{assay}_maked_dup.sorted.bw')
    shell:
        """
        bamCoverage -of bigwig -p 6 --normalizeUsing CPM -b {input.bam} -o {output.bigwig}
        """

rule call_peaks:
    input:
        bam = os.path.join(file_dir,'{cell}',assay_label,'Bam',assay_label+'_maked_dup.sorted.bam'),
        input = os.path.join(file_dir,'{cell}',"Input",'Bam',"Input"+'.sorted.bam')
    output:
        peak = os.path.join(file_dir,'{cell}',assay_label,'Peak',assay_label+'_peaks.broadPeak')
    params:
        name = assay_label,
        dir = os.path.join(file_dir,'{cell}',assay_label,'Peak')
    shell:
        """
        macs2 callpeak -t {input.bam} -c {input.input} -q 0.05 -n {params.name} -f BAM -g hs --outdir {params.dir} --broad --broad-cutoff 0.01
        """


#bamCoverage -of bigwig -p 4 --normalizeUsing CPM -b reads.bam -o {output.bigwig}
