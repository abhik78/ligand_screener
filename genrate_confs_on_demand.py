import os
import argparse
from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.conformer import ConformerGenerator, ConformerSettings


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf_mols', '-sdf', help='sdf file to generate conformers')
    parser.add_argument('--conformers_file_dir', '-out', help='directory where to write out conformmers')
    args = parser.parse_args()

    molecules = args.sdf_mols
    list_of_molecules = [m for m in MoleculeReader(molecules)]
    for mol in list_of_molecules:

        conformers = generate_confs(mol, 2, 1)
        full_directory_path = os.path.join(args.conformers_file_dir, '%s_conformers.mol2' % mol.identifier)
        with MoleculeWriter(full_directory_path) as mol_writer:
            for c in conformers:
                mol_writer.write(c.molecule)

if __name__=='__main__':
     main()