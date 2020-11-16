#conda activate covidpipe

#dry run
#snakemake --snakefile pipeRAW_fromVCF.snk.py -j 48 --config ifq=/input/vcf out=/outfolder/ -np

import subprocess, sys, os, glob 
from os.path import join
from os.path import basename
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

#############################################################################
#                                                                           #
# Description                                                               #
#                                                                           #
#############################################################################

# Extract and annotate variants from a folder containing VCF files.
# Dependences: snpEff, SnpSift, a VEP txt file containing all the possible mutations annotated and a GISAID_homoplasy_Mask.txt containing all the homoplasy positions.


# Input parameters  ------------------------------------------------------------------------
#

IFQ = config["ifq"]

workspace= config["out"]

threads=1

## Functions -------------------------------------------------------------------

plasy=[]
with open ("data/GISAID_homoplasy_Mask.txt",'r') as fi:
        for line in fi:
                if line[0]!="p":
                        plasy.append(line.split("\t")[0])
        
# Rules ------------------------------------------------------------------------
# 

SAMPLES, = glob_wildcards(IFQ+"/{sample}_samp.lofreq.bed.vcf.gz")

rule all:
    input:
        expand(workspace+'variants/{sample}_samp.lofreq.bed.ann.vcf', sample=SAMPLES),
        expand(workspace+'variants/{sample}_samp.lofreq.bed.ann.html', sample=SAMPLES),
        expand(workspace+'variants/{sample}_samp.lofreq.bed.ann.txt', sample=SAMPLES),
        workspace+'allmuta.vcf',
        workspace+'allmuta_annot.vcf'


rule annotate:
    input:
        vcf=IFQ+'/{sample}_samp.lofreq.bed.vcf.gz'
    output:
        vcf=workspace+'variants/{sample}_samp.lofreq.bed.ann.vcf',
        html=workspace+'variants/{sample}_samp.lofreq.bed.ann.html'
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
        vcf=workspace+'variants/{sample}_samp.lofreq.bed.ann.vcf'
    output:
        txt=workspace+'variants/{sample}_samp.lofreq.bed.ann.txt'
    params:
        html='{sample}_samp.lofreq.bed.ann.html'
    shell:"""
        java -Xmx4g -jar /home/mgrau/apps/snpEff/SnpSift.jar extractFields {input.vcf} CHROM POS REF ALT DP AF > {output.txt}
      """

#Including warnings
#java -Xmx4g -jar /home/mgrau/apps/snpEff/SnpSift.jar extractFields {input.vcf} CHROM POS REF ALT DP AF SB DP4 "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].EFFECT" "EFF[*].GENE" "EFF[*].CODON" > {output.txt}

rule mergeVCF:
    input:
        txt=expand(workspace+'variants/{sample}_samp.lofreq.bed.ann.txt', sample= SAMPLES)
    output:
        vcf=workspace+'allmuta.vcf'
    params:
        txt=workspace+'variants/'
    run:
        letters = ['A', 'T', 'C', 'G']
        df_muts = pd.DataFrame(columns = ["CHROM","POS","REF","ALT","DP","AF",'SAMPLE'])
        df_muts.to_csv(output.vcf, sep="\t",index=False)
	for file in list(glob.glob(params.txt+"*_samp.lofreq.bed.ann.txt")):
          df_muts = pd.read_csv(file, sep ='\t')
          #discard indels
          df_muts = df_muts[(df_muts['REF'].isin(letters))&(df_muts['ALT'].isin(letters))]
          #dicard homoplasy
          df_muts = df_muts[(~df_muts['POS'].isin(plasy)) & (df_muts['POS']>10) & (df_muts['POS'] < 29804)]
          df_muts["SAMPLE"] = os.path.splitext(os.path.basename(file))[0].split(".")[0]
          df_muts.to_csv(output.vcf, mode='a', header=False,sep="\t",index=False)

rule annotateMUTAS:
    input:
        vcf=workspace+'allmuta.vcf'	
    output:
        an=workspace+'allmuta_annot.vcf'
    params:
        vep_annot="data/annotated_mutations_VEP.tsv"
    run:
        dfvep = pd.read_csv(params.vep_annot,sep="\t")
        dfvep = dfvep[['pos', 'ref', 'alt', "csqn_type"]]
        single_muts = pd.read_csv(input.vcf,sep="\t")
        single_muts["mutation"] = single_muts['REF'].astype(str)+single_muts["POS"].astype(str) + single_muts["ALT"].astype(str)
        df = pd.merge(single_muts,dfvep, how='left',left_on=['POS','REF','ALT'],right_on=['pos','ref','alt']).dropna().drop_duplicates()
        df['branch']=""
        df = df[['mutation','branch', 'ref','alt','pos','csqn_type','AF', 'SAMPLE']]
        df.to_csv(output.an,sep="\t",index=False)

