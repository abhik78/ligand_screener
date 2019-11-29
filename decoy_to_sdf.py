from rdkit import Chem
from rdkit.Chem import AllChem
import csv
import os
import argparse
import multiprocessing

def read_picked_file(directory):
    '''
    read .picked file and create genrator
    :param directory: directory where .picked files are
    :return: generator
    '''
    if os.path.exists(directory):
        for filename in os.listdir(directory):
            if filename.endswith('.picked'):
                yield filename

def open_text_file(filename):
    '''
    read csv file and yiled rows
    :param decoy_filename: read .picked files
    :return: return rows as generator
    '''
    with open(filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        next(csv_reader) #ignore the first line
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

        AllChem.EmbedMolecule(molH, useRandomCoords=True)
        AllChem.UFFOptimizeMolecule(molH)
        molH.SetProp("_Name", "ZIN"+v)
        writer.write(molH)
    writer.close()

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--target_ids', '-t', help='target ids as single column csv file')
    parser.add_argument('--data_dir', '-dir', help='directory where the decoy .picked files are')
    parser.add_argument('--sdf_file_dir', '-out', help='directory where to write out 3d rdkit sdf files')
    args = parser.parse_args()
    return args

class Processor:
    def __init__(self,args):
        self.args = args

    def process_test_system(self,unp_id):
        decoy_dir = os.path.join(self.args.data_dir, '{}/decoys'.format(unp_id))
        print(unp_id)
        if os.path.exists(decoy_dir):
            list_of_decoy_file = [j for j in read_picked_file(decoy_dir)]

            smiles_zincid_dict = {}

            for filename in list_of_decoy_file:

                print(filename)
                picked_file = os.path.join(decoy_dir, filename)
                for i in open_text_file(filename=picked_file):
                    smiles_zincid_dict[i[0]] = i[1]

            decoy_sdf = os.path.join(self.args.sdf_file_dir, "{}_decoy_3d_rdkit.sdf".format(unp_id))
            create_3d_sdf_from_smiles(smiles_zincid_dict=smiles_zincid_dict, decoy_sdf=decoy_sdf)


def main():

    p = multiprocessing.Pool(8)

    args = parse_arguments()
    unp_id_list = [row[0] for row in open_text_file(args.target_ids)]
    print(unp_id_list)
    proc = Processor(args)

    p.map(proc.process_test_system, unp_id_list)



if __name__ == "__main__":
    main()