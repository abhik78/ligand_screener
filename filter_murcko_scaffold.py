import json
import os
import argparse
import logging
logging.basicConfig(level=logging.INFO)
import csv
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


def read_json_file(filename, choice):
    '''
    reads a json file name in this format
    {"CHEMBL21": ["Nc1ccc(cc1)S(=O)(=O)N", 628.0]} output smiles and activity dictionary based on choice
    smiles_dict = 'CHEMBL189352': 'COc1ccc2c(cnn2n1)c3ccnc(Nc4ccc(cc4)C#N)n3'
    activity_dict = 'CHEMBL359794': 1.0
    :param filename: a json file name
    :param choice: either 'smiles_dict' or 'activity_dict'
    :return: smiles or activity dictionary
    '''
    smiles_dict = {}
    activity_dict = {}
    with open(filename, 'r') as json_data:
        doc = json.load(json_data)
        for row in doc:
            smiles_dict[row] = doc[row][0]
            activity_dict[row] = doc[row][1]
    if choice == 'smiles_dict':
        return smiles_dict
    elif choice == 'activity_dict':
        return activity_dict


def apply_lead_like_filters(smiles_dict):
    '''Apply lead like filtering, exclude structures
    AlogP > 4.5
    mol wt > 450 g/mmol

    :param smiles_dict: {'CHEMBL12345' : 'c1ccccc1OC'}
    :return: filtered smiles dict
    '''
    new_dict = {}
    for k, v in smiles_dict.items():
        rdkit_mol = Chem.MolFromSmiles(v)
        if rdkit_mol:
            if Crippen.MolLogP(rdkit_mol) < 4.5 or Descriptors.ExactMolWt(rdkit_mol) < 450:
                new_dict[k] = v
    return new_dict


def getMurckoScaffold(smiles_dict):
    '''Reads a smile dictionary in this format
    'CHEMBL189352': 'COc1ccc2c(cnn2n1)c3ccnc(Nc4ccc(cc4)C#N)n3'
    Returns a dictionary of Murcko scaffolds with the corresponding molecules
    'Cc1n[nH]c2ccc(cc12)c3cncc(OC[C@@H](N)Cc4ccccc4)c3': 'CHEMBL379218'
    :param smiles_file: smiles dictionary
    :return: dictionary of scaffolds and chembl_id
    '''
    """ """
    smiles_list = smiles_dict.values()
    chembl_id_list = smiles_dict.keys()

    mols_list = [Chem.MolFromSmiles(x) for x in smiles_list]

    scaffolds = {}
    for mol, chembl_id in zip(mols_list, chembl_id_list):

        try:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold = Chem.MolToSmiles(core)


        except Exception as e:
            print("rdkit could not read {}".format(chembl_id))
        if scaffold in scaffolds:
            scaffolds[scaffold].append(chembl_id)
        else:
            scaffolds[scaffold] = []
            scaffolds[scaffold].append(chembl_id)


    return scaffolds

def get_cluster_ids_for_actives(scaffold_dict):

    '''
    Cluster the active molecules based on their scaffolds and returns
    a dictionary of chembl_id and cluster_id in this format
    {'CHEMBL3589744': 648, 'CHEMBL3589808': 648}
    :param scaffold_dict: dictionary of scaffolds
    :return: dictionary of active molecule with cluster ids
    '''

    cluster_ids = {}

    for rank, (key, values) in enumerate(scaffold_dict.items(), 1):
        for value in values:
            cluster_ids[value] = rank

    return cluster_ids

def SelectActives(mols, activities):
    '''
    :param mols: name of the molecule or ChEMBL id
    :param activities: dictionary of activity
    :return: Returns the list of molecules for each scaffold sorted based on their activities
    '''
    scaffold_activities = {}
    for mol in mols:
        scaffold_activities[mol] = activities[mol]

    return sorted(scaffold_activities, key=scaffold_activities.get)


def mt100scaffolds(scaffolds, activities, smiles):
    '''
    Returns the most active molecule for each scaffold in this format
    'N#Cc1ccc(Nc2nccc(n2)c3cnn4ncccc34)cc1': 'CHEMBL359794'
    :param scaffolds: dictionary of scaffolds
    :param activities: dictionary of activities
    :param smiles: dictionary of smiles
    :return: a dictionary of actives
    '''
    actives = {}
    for scaffold, mol_names in scaffolds.items():
        best = SelectActives(mol_names, activities)[0]
        actives[smiles[best]] = best
    return actives


def lt100scaffolds(scaffolds, activities, smiles):
    '''The number of highest affinity ligands taken from each Murcko scaffold are increased
    until we achieve at least 100 ligands or until all ligands are included
    'N#Cc1ccc(Nc2nccc(n2)c3cnn4ncccc34)cc1': 'CHEMBL359794'
    :param scaffolds: dictionary of scaffolds
    :param activities: dictionary of activities
    :param smiles: dictionary of smiles
    :return: dictionary of actives
    '''
    """ """
    actives = {}
    k = 0
    while len(actives) <= 100:
        for scaffold, mol_names in scaffolds.items():
            if len(mol_names) > k:
                best = SelectActives(mol_names, activities)[k]
                actives[smiles[best]] = best
        k += 1
        if len(actives) == len(smiles):
            break

    return actives

def write_csv(my_dict, filename):
    '''
    csv writer used to write two columns .smi file in order to DUD-E to read
    :param my_dict: dictionary
    :param filename: .smi filename
    :return:
    '''
    with open (filename, 'w', newline='') as f:
        writer = csv.writer(f, delimiter = ' ')
        for row in my_dict.items():
            writer.writerow(row)

    return filename
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--working_dir', '-dir', help='directory where the json files are')
    args = parser.parse_args()

    json_files = os.listdir(args.working_dir)

    for file in json_files:
        #print(file)
        json_file_name = file.split('.')

        json_file = os.path.join(args.working_dir, file)
        #print(json_file)
        activity_dict = read_json_file(json_file, choice='activity_dict')
        smiles_dict = read_json_file(json_file, choice='smiles_dict')

        lead_like_filtered_data = apply_lead_like_filters(smiles_dict)
        #print(lead_like_filtered_data)
        dbofscaffolds = getMurckoScaffold(lead_like_filtered_data)

        print(len(dbofscaffolds))
        if len(dbofscaffolds) >= 100:
            actives = mt100scaffolds(dbofscaffolds, activity_dict, lead_like_filtered_data)
            #print(actives)
        elif len(dbofscaffolds) < 100:
            actives = lt100scaffolds(dbofscaffolds, activity_dict, lead_like_filtered_data)
            #print(actives)
        write_csv(actives, '{}_subset_{}.smi'.format(json_file_name[0], len(dbofscaffolds)))
        #write out actives with their cluster ids
        cluster_id_dict = get_cluster_ids_for_actives(dbofscaffolds)
        write_csv(cluster_id_dict, '{}_cluster_ids.csv'.format(json_file_name[0]))

if __name__ == "__main__":
    main()