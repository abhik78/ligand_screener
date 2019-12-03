from rdkit import Chem
from rdkit.Chem import AllChem
import csv
import os
import argparse
import multiprocessing

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

def create_3d_sdf_from_smiles(smiles_chemblid_dict, active_sdf_file):

    '''
    create 3D sdf files using rdkit
    :param smiles_chemblid_dict: dictionary of smiles and molecule name
    :param decoy_sdf: 3d sdf
    :return: combined sdf file
    '''
    writer = Chem.SDWriter(active_sdf_file)

    for k, v in smiles_chemblid_dict.items():

        try:
            mol = Chem.MolFromSmiles(k)
            molH = Chem.AddHs(mol)

            AllChem.EmbedMolecule(molH)
            AllChem.UFFOptimizeMolecule(molH)
            molH.SetProp("_Name", v)
            writer.write(molH)
        except:
            print("{}".format(v))
    writer.close()


def parse_arguments():
    """

    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--target_ids', '-f', help='single column csv file with target uniprot ids, header uniprot_id')
    parser.add_argument('--data_dir', '-dir', help='directory where the .charged files are inside dude-decoys')
    parser.add_argument('--sdf_file_dir', '-out', help='directory where to write out 3d rdkit sdf files')
    args = parser.parse_args()
    return args

class Processor:
    def __init__(self,args):
        self.args = args

    def process_test_system(self,unp_id):
        decoy_dir = os.path.join(self.args.data_dir, '{}'.format(unp_id))
        if os.path.exists(decoy_dir):
            list_of_decoy_file = [j for j in read_charged_file(decoy_dir)]

            smiles_chemblid_dict = {}

            for filename in list_of_decoy_file:

                picked_file = os.path.join(decoy_dir, filename)
                for i in open_file(filename=picked_file):
                    molecule_name = i[1] + '_' + i[2]
                    print(molecule_name)
                    smiles_chemblid_dict[i[0]] = molecule_name
            active_sdf_file = os.path.join(self.args.sdf_file_dir, "{}_active_3d_rdkit.sdf".format(unp_id))
            create_3d_sdf_from_smiles(smiles_chemblid_dict=smiles_chemblid_dict, active_sdf_file=active_sdf_file)


def main():
    p = multiprocessing.Pool(8)

    args = parse_arguments()
    unp_id_list = [row[0] for row in open_file(args.target_ids)]
    print(unp_id_list)
    proc = Processor(args)

    p.map(proc.process_test_system, unp_id_list)


if __name__ == "__main__":
    main()