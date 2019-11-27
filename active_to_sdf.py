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


def parse_arguments():
    """

    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', '-f', help='single column csv file with target uniprot ids, header uniprot_id')
    parser.add_argument('--dude_decoys_dir', '-dir', help='directory where the .charged files are inside dude-decoys')
    parser.add_argument('--sdf_file_dir', '-out', help='directory where to write out 3d rdkit sdf files')
    args = parser.parse_args()
    return args

class Processor:
    def __init__(self,args):
        self.args = args

    def process_test_system(self,unp_id):
        decoy_dir = os.path.join(self.args.data_dir, '{}/decoys'.format(unp_id))
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


#multiprocessing class, modify accordingly

def main():

    args = parse_arguments()
    list_of_targets= [id for id in open_file(args.input_file)]
    #needs to add steps to browse decoys directory and find the files
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