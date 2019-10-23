import psycopg2
from psycopg2 import Error
import json
import pandas as pd
import csv
import sys

def write_csv(dict, file_name_csv):
    '''

    :param list:
    :return:
    '''
    df = pd.DataFrame(dict)
    df.to_csv(file_name_csv)
    return file_name_csv



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
        for row in results:


            smiles_inner_dict[row[1]] = row[0]
            #smiles_inner_dict.setdefault(row[1], []).append(row[0])
            set(smiles_inner_dict)

            #activity_inner_dict[row[1]] = float(row[2])
            activity_inner_dict.setdefault(row[1], []).append(float(row[2]))

            #smiles_full_dict[row[3]] = (smiles_inner_dict)
            #activity_full_dict[row[3]] = activity_inner_dict
            activity_full_dict.setdefault(row[3], activity_inner_dict)
            smiles_full_dict.setdefault(row[3], smiles_inner_dict)
        write_csv(activity_full_dict, '{}_activity.csv'.format(unp_id))
        write_csv(smiles_full_dict, '{}_smiles.csv'.format(unp_id))



