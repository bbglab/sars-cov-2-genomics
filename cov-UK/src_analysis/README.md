 README

## Analysis ##

Once the snakemake pipe finishes, there are scripts and notebooks available to perform the analysis. The pipeline output, a tab-delimited file with the annotated variants (`allmuta-annot.csv`), is required for the following analyses:

1. Spectra. 
   Notebook: `spectra_intra.ipynb`. 
   Input: `../output/paired/allmuta_annot_paired.csv`
   Output: Spectra plots

2. Expected-mutation-rate. 
   Python script for simulations: `mutrate_simulate.py`. It generates the simulations using the `../output/paired/allmuta_annot_paired.csv` file. Example call:
   `python mutrate simulate.py allmuta_annot_paired.csv segment_size start_position end_position`
   Once the simulation is complete, you can visualize the plots using the notebook `mutrate_genome.ipynb`.
   
3. Signatures. 
   Notebook to create the matrix. `sigProfiler.ipynb`.
   Input: `../output/paired/allmuta_annot_paired.csv`
   Output: Matrix required to run [SigProfilerJulia](`https://bitbucket.org/bbglab/sigprofilerjulia/src/master/`)

