from rdkit import Chem
from rdkit.Chem import AllChem
import csv
import os
import argparse

def read_charged_file(directory):
    '''
    check for ligand.charge file and return generator,
    :param directory: where ligand.charge file is
    :return: filename as generator
    '''

    for filename in os.listdir(directory):
        if filename.endswith('.charge'):
            yield filename

def open_file(filename):
    '''
    read ligand.charge file
    :param filename: ligand.charge
    :return: rows as genrator
    '''

    with open(filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            yield row

def create_3d_sdf_from_smiles(smiles_chemblid_dict, active_sdf):

    '''
    create 3D sdf files using rdkit
    :param smiles_chemblid_dict: dictionary of smiles and molecule name
    :param decoy_sdf: 3d sdf
    :return: combined sdf file
    '''
    writer = Chem.SDWriter(active_sdf)

    for k, v in smiles_chemblid_dict.items():
        #print(v)
        mol = Chem.MolFromSmiles(k)
        molH = Chem.AddHs(mol)

        AllChem.EmbedMolecule(molH)
        AllChem.UFFOptimizeMolecule(molH)
        molH.SetProp("_Name", v)
        writer.write(molH)
    writer.close()



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--active_dir', '-dir', help='directory where the decoy .charged files are')
    parser.add_argument('--sdf_file_dir', '-out', help='directory where to write out 3d rdkit sdf files')
    args = parser.parse_args()

    list_of_charged_file = [j for j in read_charged_file(args.active_dir)]
    smiles_chemblid_dict = {}
    for file in list_of_charged_file:
        charged_file = os.path.join(args.active_dir, file)

        for i in open_file(charged_file):
            smiles_chemblid_dict[i[0]] = i[1]
    active_sdf = os.path.join(args.sdf_file_dir, "active_3d_rdkit.sdf")
    create_3d_sdf_from_smiles(smiles_zincid_dict=smiles_chemblid_dict, active_sdf= active_sdf)

if __name__ == "__main__":
    main()