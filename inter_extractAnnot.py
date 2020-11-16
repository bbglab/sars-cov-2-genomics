#Extract and annotate mutations from Nextstrain JSON output.
#Two options, 1)Extract all mutations OR 2) only mutations from branches.
#By default, it extracts ALL mutations.

#How to call
#python inter_extractAnnot.py input output1 output2
#python inter_extractAnnot.py tree.json listMutations.csv listMutationsAnn.csv 

import json,sys,re, copy
import numpy as np
import sys, os, json
from bgreference import refseq 
import pandas as pd
from scipy import stats
from tqdm import tqdm
import scipy.stats as st
import statistics 

#extract list of mutations
with open (sys.argv[1]) as json_file:
		data=json.load(json_file)

#######################
##					 ##
## Extract mutations ##
##					 ##
#######################

mydi={}
def rr (d, father):
	for k in d.keys():
		#For some reason, there are nodes without branch_attrs ('NODE_0000000')
		if 'branch_attrs' in k:
			mydi[d['name']]={"father":str(father),"mutations":json.dumps(d['branch_attrs']['mutations']),"mutationsDict":d['branch_attrs']['mutations']}
		if k == 'children':
			for ch in d['children']:
				rr(ch, d['name'])
	return mydi

mydi=rr(data['tree'],0)

#One option is keeping only mutations in branches 
#(only the ones that appears as father somewhere)
# parents = [mydi[d]['father'] for d in mydi.keys()]
# res="mutation\tclade\n"
# for x in mydi.keys():
# 	if len(mydi[x]['mutations'])>2 and x in parents:
# 		for elem in mydi[x]['mutationsDict']['nuc']:
# 			res+=elem+"\n"
# 		#res+=json.dumps(mydi[x]['mutationsDict']['nuc'])+"\n"
# print(res[:-1])

#Another option is keeping all mutations (branches ands leafs)
contaBranch=0
res="mutation\tclade\n"
for x in mydi.keys():
	if len(mydi[x]['mutations'])>2:
		contaBranch+=1
		for elem in mydi[x]['mutationsDict']['nuc']:
			res+=elem+"\n"
with open (sys.argv[2],'w') as fo:
	fo.write(res[:-1])
#print(res[:-1])

print("Mutations succesfully extracted from the tree.\nStarting annotation...")

#######################
##					 ##
## 	  Annotate 	     ##
##					 ##
#######################

#list of mutations non-syn pos-ref-alt
dfvep = pd.read_csv("data/annotated_mutations_VEP.tsv",sep="\t") 
dfvep = dfvep[['pos', 'ref', 'alt', "csqn_type"]]

dfplas = pd.read_csv("data/GISAID_homoplasy_Mask.txt",sep="\t") 

#Prepare output for spectra 
single_vcf = pd.read_csv(sys.argv[2],sep="\t") 
single_vcf ['POS'] = single_vcf ['mutation'].apply(lambda x : int(''.join(x[1:-1])))
#Remove homoplasy
single_muts = single_vcf[~single_vcf['POS'].isin(dfplas['position_Start'])]
single_muts = single_muts[(single_muts['POS']>55) & (single_muts['POS']<29804)]
single_muts = single_muts.drop(columns=['POS'])
single_muts.to_csv(sys.argv[2],sep="\t",index=False)

# #Prepare output for expected
res="mutation\tbranch\tref\talt\tpos\tcsqn_type\n"
for index, row in tqdm(single_muts.iterrows()):
	ref=row[0][0]
	alt=row[0][-1]
	pos=int(row[0][1:-1])
	mutation=row[0]
	#include unique because the first ones appears multiple times.
	csqn_type=dfvep[(dfvep.pos==pos)&(dfvep.ref==ref)&(dfvep.alt==alt)].csqn_type.unique()
	res+=str(mutation)+"\t\t"+str(ref)+"\t"+str(alt)+"\t"+str(pos)+"\t"+"".join(csqn_type)+"\n"

with open (sys.argv[3],'w') as fo:
	fo.write(res)
