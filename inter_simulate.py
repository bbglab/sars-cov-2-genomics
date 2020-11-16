#1) Count observed syn mutations and simulate (10**6 times) weighted random mutations along the genome (with the observed' size).
#2) Calculate the mean number of mutations per chunk and perform a g-test. 

#how to run
#python inter_simulate.py input output chunk_size
#python inter_simulate.py listMutationsAnn.csv sim.csv 1000

import json,sys,re, copy
import numpy as np
import sys, os, json
from bgreference import refseq 
import pandas as pd
from scipy import stats
from tqdm import tqdm
import scipy.stats as st
import statistics 

jsfile="data/prob_syn_full_genome.json"
#number of simulations
N=10**6
outfile=sys.argv[2]
chunk_size=int(sys.argv[3])

def coronavirus(pos, size=1):
    return refseq("sarscov2", "0", pos, size)

def load_signature_dictionary (fjson):
    with open(fjson) as json_file: 
        return json.load(json_file) 

def get_observed_mutation_count():
    df = pd.read_csv(sys.argv[1],sep="\t")
    df=df[df["csqn_type"].isin(["upstream_gene_variant","downstream_gene_variant","synonymous_variant"])]   
    list_positions = df["pos"].values
    #return np.unique(list_positions)
    return list_positions

def count_mutation_chunk(start, pos):
    return len(np.sort(pos[(pos >= start) & (pos <= (start+chunk_size))]))

start_position = 2
chunk_start_position = 2
last_position = 29902
chunk_last_position = 29902
bases = set(["A","C","T","G"])
simulated_probs = []
simulated_position = []
observed_mutations = get_observed_mutation_count() # Get observed mutations to simulate
N_observed_mutations = len(observed_mutations)
print("Observed syn mutations: "+str(N_observed_mutations))
d_signature = load_signature_dictionary(jsfile)

#list of mutations non-syn pos-ref-alt
dfnonsyn = pd.read_csv("data/annotated_mutations_VEP.tsv",sep="\t") 
dfnonsyn=dfnonsyn[dfnonsyn["csqn_type"].isin(["stop_gained","missense_variant","start_lost"])] 
dfnonsyn = dfnonsyn[['pos', 'ref', 'alt']]

#list of homoplasy positions
dfplas = pd.read_csv("data/GISAID_homoplasy_Mask.txt",sep="\t") 
dfplas = dfplas[['position_Start']]

#Only syn
for position in range(start_position,last_position+1,1): # Iterate over the virus genome
    ref_triplet = coronavirus(position-1,3)
    ref = ref_triplet[1]
    alts =  bases - set(ref)
    for alt in alts: # for each possible mutation of the base
        #filter 1. ONLY simulated non-syn.
        #filter 2. non homoplasy
        if dfnonsyn[(dfnonsyn.pos==position)&(dfnonsyn.ref==ref)&(dfnonsyn.alt==alt)].empty and (position not in dfplas and position >55 and position < 29804):
                #syn
                prob = d_signature[ref_triplet+"_"+alt]
                simulated_probs.append(prob)
        else:
            #nonsyn
            simulated_probs.append(0)
        simulated_position.append(position)

print(len(simulated_position))
print(len(simulated_probs))

# Normalize simulated probabilities
#simulated_probs /= sum(simulated_probs)
simulated_probs = list(map(lambda x: x/sum(simulated_probs), simulated_probs))

# Now the dictionaries are ready, let randomize 
#simulated_position #array of positions
background = np.random.choice(simulated_position, size=(N, N_observed_mutations), p=simulated_probs, replace=True) 

# Now in background we have the simulations let's start counting per each chunk of interest
list_results = []

idchunk=0
for startp in tqdm(range (chunk_start_position,chunk_last_position,chunk_size)):
    observed_count_chunk=count_mutation_chunk(startp, observed_mutations)
    # Now the simulated
    chunk_simulation = []
    for simulation in background:
        chunk_simulation.append(count_mutation_chunk(startp, simulation))
    mean_count_simulated = np.nanmean(chunk_simulation) # This is the mean number of simulated mutations in the chunk! let's compare it with the observed ones
    #confidence interval
    (cih,cil) = st.t.interval(0.95, len(chunk_simulation)-1, loc=np.mean(chunk_simulation), scale=st.sem(chunk_simulation))
    std=np.std(chunk_simulation)
    a = observed_count_chunk
    b = N_observed_mutations - observed_count_chunk
    c = mean_count_simulated
    d = N_observed_mutations - mean_count_simulated  
    if a > 0:
        odds_ratio = (a/b) / (c/d)
        u,p_value = stats.power_divergence(f_obs=[a, b], f_exp=[c, d], lambda_="log-likelihood") # perform g-test
        list_results.append([startp,startp+chunk_size,idchunk, observed_count_chunk, mean_count_simulated, odds_ratio, u, p_value, cil, cih, std]) # save results
        df_results = pd.DataFrame(list_results,columns=["start_chunk","end_chunk","id_chunk","n_observed","mean_simulated","odds_ratio","u","pvalue", "cintervL", "cintervH","stdev"]) 
        df_results.to_csv(outfile,sep="\t",index=False)
    idchunk+=1