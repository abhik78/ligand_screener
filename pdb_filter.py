from ccdc import io
import requests

class SearchAPI:
    def __init__(self, search_url):
        self.search_url = search_url
        self.search_options = '&wt=json&rows=100000'


    def url_response(self, url):
        """
        Getting JSON response from URL
        :param url: String
        :return: JSON
        """
        r = requests.get(url=url)
        # Status code 200 means 'OK'
        if r.status_code == 200:
            json_result = r.json()
            return json_result
        else:
            #print(r.status_code, r.reason)
            return None

    def run_search(self):

        full_query = self.search_url
        response = self.url_response(full_query)
        return response


def create_het_chembl_mapping_dict(filename= None):
    '''
    read sdf overlays file of PDB hetcodes for each targets and find out chembl id for each hetcodes
    :param filename: sdf file for overlays
    :return: dictionary of hetcode and chembl id
    '''

    sdf_reader = io.EntryReader(filename)
    id_list = []
    for sdf in sdf_reader:
        id_list.append((sdf.molecule.identifier))

    het_list = []

    for id in id_list:
        het_list.append((id.split("_")[2]))

    het_chembl_mapping_dict = {}
    for het in het_list:
        compound_url = 'https://www.ebi.ac.uk/pdbe/api/pdb/compound/mappings/{}'.format(het)
        pdbe_qeury = SearchAPI(search_url=compound_url)
        pdb_list_dict = pdbe_qeury.run_search()
        if pdb_list_dict:
            for k, v in pdb_list_dict.items():
                for i in v:
                    if i.get('chembl_id'):
                        het_chembl_mapping_dict[i.get('chembl_id')] = het
    return het_chembl_mapping_dict

if __name__ == '__main__':
    create_het_chembl_mapping_dict()