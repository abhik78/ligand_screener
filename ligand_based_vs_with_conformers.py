#!/usr/bin/env python
#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2016-08-15: created by the Cambridge Crystallographic Data Centre
#
######################################################################
"""
    ligand_based_VS.py - simple interface to the ligand screener.
"""
import os
import argparse
import csv
from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.conformer import ConformerGenerator, ConformerSettings
from ccdc.screening import Screener
import multiprocessing

os.environ['CSDHOME'] = '/home/amukhopadhyay/not-backed-up/CCDC/CSD_2020'


###############################################################################

def open_text_file(filename):
    '''
    read csv file and yiled rows
    :param decoy_filename: read .picked files
    :return: return rows as generator
    '''
    with open(filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader)  # ignore the first line
        for row in csv_reader:
            yield row

def read_mol2_file(directory):
    '''
    read .picked file and create genrator
    :param directory: directory where .picked files are
    :return: generator
    '''
    if os.path.exists(directory):
        for filename in os.listdir(directory):
            if filename.endswith('.mol2'):
                yield filename

def setup_screener():
    """Return the settings object for the ligand screener.

    :param output_dir: Directory for screener output files
    :return: Screener settings
    """
    settings = Screener.Settings()
    # settings.output_directory = os.path.join(os.getcwd(), output_dir)
    return settings


def standardise(mol):
    """Return standardised molecule."""
    mol.standardise_aromatic_bonds()
    mol.standardise_delocalised_bonds()
    return mol


def screen_molecules(screener, mols_to_screen, activity, conformers_dir, output_name):
    """Run the ligand screener and write out the screened conformations.
    Return sorted list of ranked scores.

    :param screener:
    :param mols_to_screen: Screening set
    :param activity: 1 if the molecule is active, 0 if it's a decoy
    :param nconformers: Number of conformers to screen for each molecule in screening set
    :param nthreads: Number of threads on which to run the conformer generation
    :param output_name: File name for the result molecules
    :return: sorted list of ranked scores
    """
    screen_set = [m for m in MoleculeReader(mols_to_screen)]  ### Read the molecules to screen
    scores = []

    molwriter = MoleculeWriter(output_name)
    for mol in screen_set:
        mol_id = mol.identifier
        list_of_conformers_files = read_mol2_file(conformers_dir)
        for conformers_file in list_of_conformers_files:
            if conformers_file.startswith(mol_id):
                print(conformers_file)
                conformers_file_path = os.path.join(conformers_dir, conformers_file)
                print(conformers_file_path)

                conformers = [[x for x in MoleculeReader(conformers_file_path)]]
                print(type(conformers))
                print("yeah!!!!!! start screening")
                res = screener.screen(conformers)  # Screening step
                scores.extend([(r.score, activity, r.identifier) for r in res])
                # Write results
                for r in res:
                    molwriter.write(r.molecule)
    molwriter.close()

    return sorted(scores)


def write_scores(results, output_file):
    """Write results to an output file.

    :param results: Scores and identifiers from a screen
    :param output_file: Filename for output
    """
    with open(output_file, 'w') as f:
        for r in results:
            score, activity, identifier = r
            f.write('%.3f, %d, %s\n' % (score, activity, identifier))


def parse_command_line_args():
    """Return the command line arguments.

    Parse the command line arguments and return them."""

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-q', '--query', help='Query file')
    parser.add_argument('-a', '--actives_dir', help='Actives sdf dir')
    parser.add_argument('-d', '--decoys_dir', help='Decoys sdf dir')
    parser.add_argument('-ac', '--conformers_dir_actives', help='pregenerated conformers directory of actives')
    parser.add_argument('-dc', '--conformers_dir_decoys', help='pregenerated conformers directory of decoys')
    parser.add_argument('--target_ids', '-f', help='single column csv file with target uniprot ids, header uniprot_id')
    parser.add_argument('-o', '--output_directory', help='Output directory')

    args = parser.parse_args()

    return args


class Processor:

    def __init__(self, args):
        self.args = args

    def test_func(self, unp_id):
        # Parse arguments
        overlay_dir = self.args.query
        active_dir = self.args.actives_dir
        decoy_dir = self.args.decoys_dir
        active_conformers_dir = os.path.join(self.args.conformers_dir_actives, '{}'.format(unp_id))
        decoy_conformers_dir = os.path.join(self.args.conformers_dir_decoys, '{}'.format(unp_id))

        actives_to_screen = os.path.join(active_dir, '{}_active_3d_rdkit.sdf'.format(unp_id))
        print(actives_to_screen)
        decoys_to_screen = os.path.join(decoy_dir, '{}_decoy_3d_rdkit.sdf'.format(unp_id))
        print(decoys_to_screen)



        output_dir = self.args.output_directory
        complete_output_dir = os.path.join(output_dir, '{}'.format(unp_id))
        overlay_mol = os.path.join(overlay_dir, '{}.sdf'.format(unp_id))
        print(overlay_mol)
        query = [m for m in MoleculeReader(overlay_mol)]  # Read the query mol or overlay of mols
        print("start screening")
        settings = setup_screener()
        os.makedirs(complete_output_dir)
        os.chdir(complete_output_dir)
        screener = Screener(query, settings=settings)  # Generate fields around the input query

        output_name_actives = os.path.join(complete_output_dir, "{}_actives_screened.mol2".format(unp_id))
        print(output_name_actives)
        output_name_decoys = os.path.join(complete_output_dir, "{}_decoys_screened.mol2".format(unp_id))
        print(output_name_decoys)
        actives_scores = screen_molecules(screener, actives_to_screen, 1, active_conformers_dir,
                                          output_name_actives)  ### Screen set of actives
        decoys_scores = screen_molecules(screener, decoys_to_screen, 0, decoy_conformers_dir,
                                         output_name_decoys)  ### Screen set of decoys
        print("writing scores")
        all_data = actives_scores
        all_data.extend(decoys_scores)
        screening_scores = sorted(all_data)

        output_name_scores = os.path.join(complete_output_dir, "{}_screening_scores.csv".format(unp_id))
        write_scores(screening_scores, output_name_scores)


def main():
    p = multiprocessing.Pool(8)
    args = parse_command_line_args()

    unp_id_list = [row[0] for row in open_text_file(args.target_ids)]
    proc = Processor(args)
    p.map(proc.test_func, unp_id_list)


if __name__ == '__main__':
    main()