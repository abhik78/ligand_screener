import psycopg2
from psycopg2 import Error
import json
import csv
import sys
from collections import defaultdict

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

def write_json(filename, my_dict):
    with open(filename, 'w') as f:
        json.dump(my_dict, f)
    return filename


def write_csv(my_dict, filename):
    with open (filename, 'w', newline='') as f:
        writer = csv.writer(f)
        for row in my_dict.items():
            writer.writerow(row)

    return filename

def query_database(sql):
    '''

    :param sql:
    :return:
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

with open(sys.argv[1]) as csv_file:
    reader = csv.DictReader(csv_file)

    for line in reader:
        unp_id = line['uniprot']

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

        smiles_full_dict = {}
        activity_full_dict = {}
        smiles_inner_dict = {}
        activity_inner_dict = {}
        d = {}
        for row in results:


            smiles_inner_dict[row[1]] = row[0]
            set(smiles_inner_dict)

            activity_inner_dict.setdefault(row[1], []).append(float(row[2]))
            activity_full_dict.setdefault(unp_id, activity_inner_dict)
            smiles_full_dict.setdefault(unp_id, smiles_inner_dict)

            # find highest activity
            for k,v in activity_inner_dict.items():
                d[k] = max(v)

        #print(d)
        #write_csv(d, '{}_activity.csv'.format(unp_id))
        #write_csv(smiles_inner_dict, '{}_smiles.csv'.format(unp_id))


        new_dict = merge_two_dict(smiles_inner_dict, d)
        print(new_dict)
        write_json('{}.json'.format(unp_id), new_dict)
