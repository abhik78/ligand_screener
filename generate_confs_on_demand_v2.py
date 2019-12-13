import os
import argparse
from ttictoc import TicToc
import glob
from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.conformer import ConformerGenerator, ConformerSettings
#from multiprocessing.pool import ThreadPool as Pool
import multiprocessing
from pathos.pools import ProcessPool

def read_sdf_file(directory):
    '''
    read .picked file and create genrator
    :param directory: directory where .picked files are
    :return: generator
    '''
    if os.path.exists(directory):
        for filename in os.listdir(directory):
            if filename.endswith('.sdf'):
                yield filename


def standardise(mol):
    '''
    standardised mols using ccdc api options
    :param mol: ccdc molecule
    :return: ccdc molecule
    '''
    mol.standardise_aromatic_bonds()
    mol.standardise_delocalised_bonds()
    return mol

def generate_confs(mol, nconformers):
    '''
    generate conformmers
    :param mol: ccdc molecules
    :param nconformers: number of conformers to generate
    :param nthreads: number of threads to use
    :return:
    '''
    settings = ConformerSettings()
    settings.max_conformers = nconformers
    gen = ConformerGenerator(settings, nthreads=1)
    standardised_mol = standardise(mol)
    confs = gen.generate(standardised_mol)

    return confs





def test_func(sdf_file):

        sdf_file_name = sdf_file.split('_')[0]
        print(sdf_file_name)
        #file_path = '/home/amukhopadhyay/ligand_screener_data/conf_test'
        #conf_dir = os.makedirs('/home/amukhopadhyay/ligand_screener_data/conf_test')
        sdf_dir = '/home/amukhopadhyay/ligand_screener_data/test_sdf'
        full_directory_path = os.path.join(sdf_dir, '{}'.format(sdf_file_name))
        os.makedirs(full_directory_path)
        os.chdir(full_directory_path)


        molecule_object_from_sdf_file = MoleculeReader(os.path.join(sdf_dir, sdf_file))
        #list_of_molecules = [m for m in molecule_object_from_sdf_file]
        for mol in molecule_object_from_sdf_file:
            conformers = generate_confs(mol, 2)
            with MoleculeWriter('%s_conformers.mol2' % mol.identifier) as mol_writer:
                for c in conformers:
                    mol_writer.write(c.molecule)






def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf_dir', '-dir', help='sdf file to generate conformers')
    #parser.add_argument('--conformers_file_dir', '-out', help='directory where to write out conformmers')
    #parser.add_argument('--number_of_conformers', '-n', help='number of conformers to generate')
    args = parser.parse_args()
    return args

def main():
    p = multiprocessing.Pool(8)

    args = parse_arguments()
    sdf_dir = os.path.join(args.sdf_dir)
    list_of_sdf_files = [filename for filename in read_sdf_file(sdf_dir)]
    #print((list_of_sdf_files))

    t = TicToc()
    if list_of_sdf_files:
        p.map(test_func, list_of_sdf_files)
    t.tic()





        #molecule_object_from_sdf_file = MoleculeReader(os.path.join(sdf_dir, file))
        #list_of_molecules = [m for m in molecule_object_from_sdf_file]

        # for mol in molecule_object_from_sdf_file:
        #     conformers = generate_confs(mol, int(args.number_of_conformers), 8)
        #
        #     with MoleculeWriter('%s_conformers.mol2' % mol.identifier) as mol_writer:
        #         for c in conformers:
        #             mol_writer.write(c.molecule)






    t.toc()
    print(t.elapsed)

if __name__=='__main__':
     main()