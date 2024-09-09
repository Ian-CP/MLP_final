# MLP_final (introduction goes here)
Modular Linux Production, final version (1.0)
<br/>

# Requirements
* Linux - Ubuntu > 18.0<br/>

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

# Usage
Basic inputs: <br/>
```>> $ python modular_linux_production_AA4.py -i <input_folder> -o <output_folder>```<br/><br/>
-i -> Directory. Contains the CSV files with the crosslinked results.<br/>
-o -> Directory. Will contain all of the outputs from the code, separated by each pair of linked proteins<br/>
<br/>
Other inputs: <br/>
-h -> Shows help info inside the terminal.<br/>
-f -> CSV file name. Will only analyze the named file. Must be present in input folder<br/>
-x -> Crosslinker name. Used to manually change all crosslinkers.<br/>
-m -> Choice of \<SWISS\> or \<AF\>. Will determine the preferred database to fetch PDB files. Either from SwissModel/Uniprot (SWISS) or AlphaFold repository (AF).<br/>
<br/>
Advanced Inputs: <br/>
--remove_buried -> Boolean, default True. If True, will use FreeSaSa to determine which residues are solvent innaccessible and remove them from consideration. <br/>
--reject_interspecies -> Boolean/String, defaults False. If True, will reject protein pairs that are of different species. If a species is named, will only consider pairs of that species.<br/>
--protein -> String. Will limit the analysis to protein pairs that contain at least 1 of the desired protein.<br/>
--change_species -> String. Will manually overwrite the identified species.<br/>
--crosslinker_size -> Integer. Will manually overwrite the crosslinker arm distance. Will not change the name of the crosslinker<br/>
--continue_from -> \<Protein1_Protein2\>. Will skip protein pairs until reaching the desired one, then will continue the analysis normally. Useful in cases where analysis was prematurelly stopped. <br/>
<br/>
--fix_partial_prot -> Boolean, default True. If False, will attempt to reconstruct the whole protein using PDBFixer. Not recommended<br/>
--skip_simulation1 -> Boolean, default False. If True, will skip checking interprotein crosslinks. Automatically done in intraprotein cases.<br/>
--skip_haddock -> Boolean, default False. If True, will skip the docking portion of the code using haddock. Either provide your own docked structure under OpenMM_sim/\<Protein1_Protein2.pdb\> or skip second simulation.<br/>
--skip_simulation2 -> Boolean, default False. If True, will skip checking for interprotein crosslinks.<br/>
<br/>
* Notes<br/>
If a specific PDB/Fasta files are to be used, simply place the files ```Protein_Species.pdb``` or ```Protein_Species.fa``` under the ```pdb_fa_files/``` folder.<br/>
Same is available for the docked structure by placing the file ```Protein1_Protein2.pdb``` under ```OpenMM_sim``` (```CG_OpenMM``` for the coarse-grained protocol) <br/>
