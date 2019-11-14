import psycopg2
from psycopg2 import Error
import json
import csv
import sys
from collections import defaultdict
import argparse
import os

def merge_two_dict(d1, d2):
    '''
    d1 = {1: 2, 3: 4}
    d2 = {1: 6, 3: 7}
    returns defaultdict(<class 'list'>, {1: [2, 6], 3: [4, 7]}
    :return:
    '''

    dd = defaultdict(list)

    for d in (d1, d2): # you can list as many input dicts as you want here
        for key, value in d.items():
            dd[key].append(value)
    return dd

def write_json(json_dir, filename, my_dict):
    '''
    writes a json file in this format
    {"CHEMBL379218": ["Cc1n[nH]c2ccc(cc12)c3cncc(OC[C@@H](N)Cc4ccccc4)c3", 46.0]}

    :param json_dir: output directory
    :param filename: json filename
    :param my_dict: dictonary to write
    :return:
    '''
    complete_filepath = os.path.join(json_dir, filename)
    with open(complete_filepath, 'w') as f:
        json.dump(my_dict, f)
    return filename


def write_csv(my_dict, filename):
    '''
    function not used presently, may be needed in future, will keep it for now
    :param my_dict: dictionary to write
    :param filename:
    :return: csv file
    '''
    with open (filename, 'w', newline='') as f:
        writer = csv.writer(f)
        for row in my_dict.items():
            writer.writerow(row)

    return filename

def query_database(sql):
    '''
    :param sql: sql query
    :return: rows
    '''
    try:
        connection = psycopg2.connect(database="chembl_25"
                                      , user="postgres",
                                      password="abhik1234",
                                      host="127.0.0.1",
                                      port="5432")
        cursor = connection.cursor()
        cursor.execute(sql)
        rows = cursor.fetchall()

    except (Exception, psycopg2.Error) as error :
        print("Error while connecting to PostgreSQL", error)

    finally:
        if (connection):
            cursor.close()
            connection.close()

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', '-f', help='single column csv file with uniprot ids, header uniprot_id')
    parser.add_argument('--json_dir', '-d', help = 'path of output directory')

    args = parser.parse_args()

    with open(args.input_file) as csv_file:
        reader = csv.DictReader(csv_file)

        for line in reader:
            unp_id = line['uniprot_id']
            print(unp_id)

            sql1 = '''select a.canonical_smiles, b.chembl_id, d.standard_value, f.chembl_id as target_chembl_id from compound_structures a 
                        INNER JOIN molecule_dictionary b 
                        on a.molregno = b.molregno 
                        INNER JOIN compound_properties c 
                        on a.molregno = c.molregno  
                        INNER JOIN activities d 
                        on a.molregno = d.molregno  
                        INNER JOIN assays e                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                        on d.assay_id = e.assay_id
                        INNER JOIN target_dictionary f
                        on e.tid = f.tid
                        INNER JOIN target_components h
                        on f.tid = h.tid
                        INNER JOIN component_sequences i
                        on i.component_id = h.component_id
                        where c.num_ro5_violations <= 1
                        and d.standard_type = 'Ki' and d.standard_units = 'nM' and d.standard_value < 1000 and d.standard_relation = '='
                        and e.assay_type = 'B' and e.confidence_score = 9
                        and f.target_type = 'SINGLE PROTEIN'
                        and i.accession = '{}'           
                '''.format(unp_id)

            results = query_database(sql1)

            smiles_inner_dict = {}
            activity_inner_dict = {}
            highest_activity_dict = {}
            for row in results:


                smiles_inner_dict[row[1]] = row[0]
                set(smiles_inner_dict)

                activity_inner_dict.setdefault(row[1], []).append(float(row[2]))

                # find highest activity
                for k,v in activity_inner_dict.items():
                    highest_activity_dict[k] = max(v)

            new_dict = merge_two_dict(smiles_inner_dict, highest_activity_dict)
            print(new_dict)
            write_json(args.json_dir, '{}.json'.format(unp_id), new_dict)

if __name__ == "__main__":
    main()