#conda activate covidpipe
#All single samples July 2020, ~8k samples
#snakemake --snakefile pipeRAW_notemp_SINGLE.snk.py -j 48 --config ifq=/workspace/datasets/sars_cov_2/cov-UK/raw/single/ fa=src/NC_045512.fasta out=/workspace/datasets/sars_cov_2/cov-UK/out_singleJuly/ -np

import subprocess, sys, os, glob 
from os.path import join
from os.path import basename
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import pandas as pd

#############################################################################
#                                                                           #
# Description                                                               #
#                                                                           #
#############################################################################

#https://covid19.galaxyproject.org/genomics/4-variation/#what-s-the-point
# 1. Filter reads with mapping quality of at least 20, that were mapped as proper pairs
# 2. Map all reads against COVID-19 reference NC_045512.2 using bwa mem
# 3. Mark duplicate reads with picard markduplicates
# 4. Perform realignments using lofreq viterbi
# 5. Call variants using lofreq call
# 6. Annotate variants using snpeff against database created from NC_045512.2 GenBank file
# 7. Convert VCFs into tab delimited dataset
# 8. Merge all tab delimited vcfs
# 9. Annotate the tab delimited vcf using VEP annotation.

# Input parameters  ------------------------------------------------------------------------
#

IFQ = config["ifq"]

REF = config["fa"]

workspace= config["out"]

threads=1

#snpEff DB
#snpDB=""


## Functions -------------------------------------------------------------------

plasy=[]
with open ("src/GISAID_homoplasy_Mask.txt",'r') as fi:
        for line in fi:
                if line[0]!="p":
                        plasy.append(line.split("\t")[0])

#Check file-names
#SAMPLES=[]
#for file in glob.glob(IFQ+"*.fastq.gz"):
#  SAMPLES.append('.'.join(file.split("/")[-1].split(".")[:-1]))
        
# Rules ------------------------------------------------------------------------
# 

#PAIRED READS
SAMPLES, = glob_wildcards(IFQ+"/{sample}.fastq.gz")

rule all:
    input:
        #expand(workspace+"qualtrim/{sample}.R1.paired.fastq", sample=SAMPLES),
        #expand(workspace+"qualtrim/{sample}.R2.paired.fastq", sample=SAMPLES),
        #expand(workspace+'alignments/{sample}.bam', sample=SAMPLES),
        #expand(workspace+'alignments/{sample}.dup.bam', sample=SAMPLES),
        #expand(workspace+'alignments/{sample}.dup.realg.bam', sample=SAMPLES),
        expand(workspace+'variants/{sample}.dup.realg.vcf', sample=SAMPLES),
        expand(workspace+'variants/{sample}.dup.realg.ann.vcf', sample=SAMPLES),
        expand(workspace+'variants/{sample}.dup.realg.ann.html', sample=SAMPLES),
        expand(workspace+'variants/{sample}.dup.realg.ann.txt', sample=SAMPLES),
        workspace+'allmuta.vcf',
        workspace+'allmuta_annot.vcf'

#QUALITY FILTER
# rule filter:
#    input:
#         #fastq=join(FASTA_DIR, PATTERN)
#         faR1=expand(IFQ+"{{sample}}_{pair}.fastq.gz", pair=["1"]),
#         faR2=expand(IFQ+"{{sample}}_{pair}.fastq.gz", pair=["2"])
#    output:
#         #files=workspace+'qualtrim/{sample}.fq'
#         R1out=workspace+"qualtrim/{sample}.R1.paired.fastq",
#         R2out=workspace+"qualtrim/{sample}.R2.paired.fastq",
#         R1out_unpaired=workspace+"qualtrim/{sample}.R1.unpaired.fastq",
#         R2out_unpaired=workspace+"qualtrim/{sample}.R2.unpaired.fastq"
#    shell:"""
#       trimmomatic PE -phred33 {input.faR1} {input.faR2} {output.R1out} {output.R1out_unpaired} {output.R2out} {output.R2out_unpaired} ILLUMINACLIP:/home/mgrau/apps/Trimmomatic-0.39/adapters/TruSeq2-SE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 MINLEN:100 HEADCROP:10 TRAILING:5
#   """ 

rule filter:
   input:
        fa=expand(IFQ+"{{sample}}.fastq.gz", sample=SAMPLES)
   output:
        faout=workspace+"qualtrim/{sample}.fastq"
   shell:"""
        cutadapt --minimum-length 20 --cut 5 --cut -5 -o {output.faout} {input.fa}
        """

#Mapping reads to the new refs. TODO change bwa -> minimap2
rule mapping:
    input:
        refFasta=REF,
        fain=workspace+"qualtrim/{sample}.fastq",
    output:
        sam=workspace+'alignments/{sample}.sam',
        bam=workspace+'alignments/{sample}.bam'
    params:
        threads=threads
    shell:"""
        rgid=$(echo {input.fain}  | md5sum | cut -d " " -f1)
        rgpu=${{rgid}}.PU
	    bwa index {input.refFasta}
        bwa mem -R "@RG\\tID:$rgid\\tPL:illumina\\tPU:$rgpu\\tSM:{input.fain}" {input.refFasta} -t {params.threads} {input.fain} > {output.sam} 
        samtools view -bT {input.refFasta} {output.sam} | samtools sort > {output.bam}
        rm {input.fain}
     """

#Mark duplicates
rule markDuplicates:
    input:
        bam=workspace+'alignments/{sample}.bam',
        sam=workspace+'alignments/{sample}.sam'
    output:
        bam=workspace+'alignments/{sample}.dup.bam',
        metrics=workspace+'alignments/{sample}.metrics.txt'
    params:
    shell:"""
        picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics}
        rm {input.bam}
        rm {input.sam}
    """

#lofreq viterbi
rule realign:
    input:
        refFasta=REF,
        bam=workspace+'alignments/{sample}.dup.bam'
    output:
        bam=workspace+'alignments/{sample}.dup.realg.bam'
    shell:"""
        lofreq viterbi -o {output.bam} -f {input.refFasta} {input.bam}
        rm {input.bam}
    """

#Call variants
rule calling:
    input:
        refFasta=REF,
        bam=workspace+'alignments/{sample}.dup.realg.bam'
    output:
        vcf=workspace+'variants/{sample}.dup.realg.vcf'
    shell:"""
        lofreq call -o {output.vcf} -f {input.refFasta} {input.bam}   
        rm {input.bam}
    """

rule annotate:
    input:
        vcf=workspace+'variants/{sample}.dup.realg.vcf'
    output:
        vcf=workspace+'variants/{sample}.dup.realg.ann.vcf',
        html=workspace+'variants/{sample}.dup.realg.ann.html'
    params:
        html='{sample}.dup.realg.ann.html',
    shell:"""
        #with our DB
        java -Xmx4g -jar /home/mgrau/apps/snpEff/snpEff.jar sarscov2 -s {output.html} {input.vcf} > {output.vcf}
        #with snpeff DB
        #java -Xmx4g -jar /home/mgrau/apps/snpEff/snpEff.jar NC_045512.2 {input.vcf} > {output.vcf}
    """

rule extract:
    input:
        vcf=workspace+'variants/{sample}.dup.realg.ann.vcf'
    output:
        txt=workspace+'variants/{sample}.dup.realg.ann.txt'
    params:
        html='{sample}.dup.realg.ann.html'
    shell:"""
        java -Xmx4g -jar /home/mgrau/apps/snpEff/SnpSift.jar extractFields {input.vcf} CHROM POS REF ALT DP AF > {output.txt}
      """

#Including warnings
#java -Xmx4g -jar /home/mgrau/apps/snpEff/SnpSift.jar extractFields {input.vcf} CHROM POS REF ALT DP AF SB DP4 "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].EFFECT" "EFF[*].GENE" "EFF[*].CODON" > {output.txt}

rule mergeVCF:
    input:
        txt=expand(workspace+'variants/{sample}.dup.realg.ann.txt', sample= SAMPLES)
    output:
        vcf=workspace+'allmuta.vcf'
    params:
        txt=workspace+'variants/'
    run:
        letters = ['A', 'T', 'C', 'G']
        df_muts = pd.DataFrame(columns = ["CHROM","POS","REF","ALT","DP","AF",'SAMPLE'])
        df_muts.to_csv(output.vcf, sep="\t",index=False)
	for file in list(glob.glob(params.txt+"*.dup.realg.ann.txt")):
          df_muts = pd.read_csv(file, sep ='\t')
          #discard indels
          df_muts = df_muts[(df_muts['REF'].isin(letters))&(df_muts['ALT'].isin(letters))]
          #dicard homplasy
          df_muts = df_muts[(~df_muts['POS'].isin(plasy)) & (df_muts['POS']>10) & (df_muts['POS'] < 29804)]
          df_muts["SAMPLE"] = os.path.splitext(os.path.basename(file))[0].split("_")[0]
          df_muts.to_csv(output.vcf, mode='a', header=False,sep="\t",index=False)

rule annotateMUTAS:
    input:
        vcf=workspace+'allmuta.vcf'	
    output:
        an=workspace+'allmuta_annot.vcf'
    params:
        vep_annot="src/annotated_mutations_VEP.tsv"
    run:
        dfvep = pd.read_csv(params.vep_annot,sep="\t")
        dfvep = dfvep[['pos', 'ref', 'alt', "csqn_type"]]
        single_muts = pd.read_csv(input.vcf,sep="\t")
        single_muts["mutation"] = single_muts['REF'].astype(str)+single_muts["POS"].astype(str) + single_muts["ALT"].astype(str)
        df = pd.merge(single_muts,dfvep, how='left',left_on=['POS','REF','ALT'],right_on=['pos','ref','alt']).dropna().drop_duplicates()
        df['branch']=""
        df = df[['mutation','branch', 'ref','alt','pos','csqn_type','AF', 'SAMPLE']]
        df.to_csv(output.an,sep="\t",index=False)
