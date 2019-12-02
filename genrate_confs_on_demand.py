import os
import argparse
from ttictoc import TicToc

from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.conformer import ConformerGenerator, ConformerSettings

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

def generate_confs(mol, nconformers, nthreads):
    '''
    generate conformmers
    :param mol: ccdc molecules
    :param nconformers: number of conformers to generate
    :param nthreads: number of threads to use
    :return:
    '''
    settings = ConformerSettings()
    settings.max_conformers = nconformers
    gen = ConformerGenerator(settings, nthreads=nthreads)
    standardised_mol = standardise(mol)
    confs = gen.generate(standardised_mol)

    return confs
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf_dir', '-dir', help='sdf file to generate conformers')
    parser.add_argument('--conformers_file_dir', '-out', help='directory where to write out conformmers')
    parser.add_argument('--number_of_conformers', '-n', help='number of conformers to generate')
    args = parser.parse_args()
    return args

class Conformer_generator:
    def __init__(self, args):
        self.args = args

    def generate_conformer(self, mol):

        conformers = generate_confs(mol, int(self.args.number_of_conformers), 1)

        #full_file_path = os.path.join(directory, '%s_conformers.mol2' % mol.identifier)
        with MoleculeWriter('%s_conformers.mol2' % mol.identifier) as mol_writer:
            for c in conformers:
                mol_writer.write(c.molecule)

def work(args, mol):
    proc = Conformer_generator(args)
    proc.generate_conformer(mol)

def main():
    args = parse_arguments()
    sdf_dir = os.path.join(args.sdf_dir)
    list_of_sdf_files = [filename for filename in read_sdf_file(sdf_dir)]

    proc = Conformer_generator(args)
    t = TicToc()
    t.tic()

    for file in list_of_sdf_files:
        sdf_file_name = file.split('_')[0]
        full_directory_path = os.path.join(args.conformers_file_dir, '{}'.format(sdf_file_name))
        os.makedirs(full_directory_path)
        os.chdir(full_directory_path)
        list_of_molecules = []
        for m in MoleculeReader(os.path.join(sdf_dir, file)):
            list_of_molecules.append(m)

        print(type(list_of_molecules))
        for mol in list_of_molecules:
            print(mol.identifier)
            proc.generate_conformer(mol)
    t.toc()
    print(t.elapsed)

if __name__=='__main__':
     main()