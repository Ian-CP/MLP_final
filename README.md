# MLP_final (introduction goes here)
Modular Linux Production, final version (1.0)
<br/>

# Requirements
* Linux Environment

* Python > 3.7<br/>
* OpenMM 8.1.1<br/>
* Haddock3 3.0 (dictates biopython version)<br/>
* biopython 1.78<br/>
* freesasa 2.2.1<br/>
* MDAnalysis 2.3.0<br/>
* mdTraj 1.10<br/>
* pdbfixer 1.9.0<br/>
* pandas 1.5.1<br/>
<br/>

If coarse-grained protocol is to be used<br/>
* martini_openmm 0.1<br/>
* vermouth (Martinizer) 0.10.0<br/>
<br/>

# Examples
NF1 - Homodimer <br/>
```
$ python modular_linux_production_AA6.py -i input -o output -f XLMS_NF1_filtered.csv --protein NF1 --reject_interspecies HUMAN
```
<br/>

<br/>PREX1 - Single protein<br/>
```
$ python modular_linux_production_AA6.py -i input -o output --protein PREX1 --reject_interspecies HUMAN --skip_haddock --skip_simulation2
```
<br/>

<br/>TSC1/TSC2/TBCD7 - multimer <br/>
```
$ python modular_linux_production_AA6.py -i input -o output -f TSC_WIPI3_multimer_peptides.csv --reject_interspecies HUMAN --skip_simulation2
```
<br/>

# Usage
Basic inputs: <br/>
```$ python modular_linux_production_AA4.py -i <input_folder> -o <output_folder>```<br/><br/>
-i -> Directory. Contains the CSV files with the crosslinked results.<br/>
-o -> Directory. Will contain all of the outputs from the code, separated by each pair of linked proteins<br/>

<br/>
* Notes<br/>
If a specific PDB/Fasta files are to be used, simply place the files ```Protein_Species.pdb``` or ```Protein_Species.fa``` under the ```pdb_fa_files/``` folder.<br/>
Same is available for the docked structure by placing the file ```Protein1_Protein2.pdb``` under ```OpenMM_sim``` (```CG_OpenMM``` for the coarse-grained protocol) <br/>
