# ligand_screener
Scripts to generate data for actives and run the ligand screener validation pipeleine. 

Description of scripts

1. find_actives_from_chembl.py - Search the ChEMBL database and create the datasets of targets (json file). 
2. filter_murcko_scaffold.py - Read the json files of the targets and 
  a) apply the lead like filters
  b)remove chembl molecule that is present in PDB overlays (pdb_filter.py)
  c)generate murcko scaffold
  d)does scaffold based clustering
  e)does scaffold filtering
  f)writes out smiles file of actives for each targets and files with scaffold cluster ids
3. decoy_to_sdf.py and active_tosdf.py - generate 3D sdf from smiles using rdkit
4. generate_confs_on_demand.py - generate conformers for each sdf molecule using ccdc conformer generator
5. calculate_scores.py - calculate different enrichment scores and create ROC curve
