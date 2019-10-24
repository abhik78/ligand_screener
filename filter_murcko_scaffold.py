from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import csv

def read_smiles_file(filename):
    '''
    read a two columns csv file with smiles and returns a dictionary in this format
    'CHEMBL4076947': 'Cc1cc(Nc2nccc3nc([nH]c23)c4c(Cl)cc(cc4Cl)C#N)ncn1'
    :param filename: two columns csv file
    :return: dictionary
    '''
    smiles_dict = {}
    with open(filename) as csv_file:
        rows = csv.reader(csv_file)
        for row in rows:
            smiles_dict[row[1]] = row[0]
    return smiles_dict

def get_activity_dict(filename):
    '''
    read two columns csv file with activity and returns a dcitionary in this format
    'CHEMBL379218': '46.0'
    :param filename: two columns csv file
    :return: dictionary
    '''
    activity_dict = {}
    with open(filename) as csv_file:
        rows= csv.reader(csv_file)
        for row in rows:
            activity_dict[row[0]] = row[1]
    print(activity_dict)
    return activity_dict

def getMurckoScaffold(smiles_file):
    '''Reads two coulms csv file with smiles and chembl_id in it
    Returns a dictionary of Murcko scaffolds with the corresponding molecules
    'Cc1n[nH]c2ccc(cc12)c3cncc(OC[C@@H](N)Cc4ccccc4)c3': 'CHEMBL379218'
    :param smiles_file: two columns csv file
    :return: dictionary of scaffolds and chembl_id
    '''
    """ """
    smiles_list = list(read_smiles_file(smiles_file).values())
    chembl_id_list = list(read_smiles_file(smiles_file).keys())

    mols_list = [Chem.MolFromSmiles(x) for x in smiles_list]

    scaffolds = {}
    for mol, y in zip(mols_list, chembl_id_list):


        core = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold = Chem.MolToSmiles(core)

        if scaffold in scaffolds:
             scaffolds[scaffold].append(y)
        else:
             scaffolds[scaffold] = []
             scaffolds[scaffold].append(y)

    return scaffolds

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

dbofscaffolds = getMurckoScaffold('P24941_smiles.csv')
smiles = read_smiles_file('P24941_smiles.csv')
ki = get_activity_dict('P24941_activity.csv')
print(len(dbofscaffolds))
if len(dbofscaffolds) >= 100:
    actives = mt100scaffolds(dbofscaffolds, ki, smiles)
    print(actives)
elif len(dbofscaffolds) < 100:
    actives = lt100scaffolds(dbofscaffolds, ki, smiles)
    #print(actives)