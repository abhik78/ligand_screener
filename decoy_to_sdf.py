from rdkit import Chem
from rdkit.Chem import AllChem
import csv
import os
import argparse

def read_picked_file(directory):
    '''
    read .picked file and create genrator
    :param directory: directory where .picked files are
    :return: generator
    '''

    for filename in os.listdir(directory):
        if filename.endswith('.picked'):
            yield filename

def open_decoy_picked(decoy_filename):
    '''
    read csv file and yiled rows
    :param decoy_filename: read .picked files
    :return: return rows as generator
    '''
    with open(decoy_filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        next(csv_reader)
        for row in csv_reader:
            yield row

def create_3d_sdf_from_smiles(smiles_zincid_dict, decoy_sdf):

    '''
    create 3DF sdf files using rdkit
    :param smiles_zincid_dict: dictionary of smiles and molecule name
    :param decoy_sdf: 3d sdf
    :return: combined sdf file
    '''
    writer = Chem.SDWriter(decoy_sdf)

    for k, v in smiles_zincid_dict.items():
        #print(v)
        mol = Chem.MolFromSmiles(k)
        molH = Chem.AddHs(mol)

        AllChem.EmbedMolecule(molH)
        AllChem.UFFOptimizeMolecule(molH)
        molH.SetProp("_Name", "ZIN"+v)
        writer.write(molH)
    writer.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--decoys_dir', '-dir', help='directory where the decoy .picked files are')
    parser.add_argument('--sdf_file_dir', '-out', help='directory where to write out 3d rdkit sdf files')
    args = parser.parse_args()

    list_0f_decoy_file = [j for j in read_picked_file(args.decoys_dir)]

    smiles_zincid_dict = {}
    for file in list_0f_decoy_file:
        decoy_file = os.path.join(args.decoys_dir, file)
        for i in open_decoy_picked(decoy_file):
            smiles_zincid_dict[i[0]] = i[1]
    #smiles_zincid_dict = {i[0]:i[1] for file in list_0f_decoy_file for i in open_decoy_picked(os.path.join(args.decoys_dir, file))}
    decoy_sdf = os.path.join(args.sdf_file_dir, "decoy_3d_rdkit.sdf")
    create_3d_sdf_from_smiles(smiles_zincid_dict=smiles_zincid_dict, decoy_sdf= decoy_sdf)

if __name__ == "__main__":
    main()