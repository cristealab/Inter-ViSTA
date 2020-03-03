
import requests
import pandas as pd
import numpy as np
import io
from Bio import SeqIO

##############################################################################

def list2string(l, sep=' '):
        
        query = ''

        for item in l:
            try:
                query = query + item + sep

            except TypeError: # exception to deal with NaNs
                pass
            
        return query
    
##############################################################################
    
class UniprotQuery:
    def __init__(self):
                
        self.url = 'https://www.uniprot.org/uploadlists/'

    def map_proteins(self, proteins, to_ = 'STRING_ID'):
        
        params = {'query': list2string(proteins, sep=' '),
                  'format': 'tab',
                  'from': 'ACC+ID',
                  'to': to_}
        
        r = requests.post(self.url, data = params)
        
        return pd.read_csv(io.StringIO(r.text), sep = '\t', header = 0, index_col = None)
    
    def retrieve_entries(self, proteins):
        
        params = {'query': list2string(proteins, sep=' '), 
                  'format': 'xml',
                  'from': 'ACC+ID',
                  'to': 'ACC'}

        r = requests.post(self.url, data = params)
             
        return r
        
    def parse_xml(self, response):
        
        data = {}

        for n, rec in enumerate(SeqIO.parse(io.StringIO(response.text), 'uniprot-xml'), 1):

            keys = rec.annotations.keys()

            if 'comment_subcellularlocation_location' in keys:
                locs = rec.annotations['comment_subcellularlocation_location']
            else:
                locs = []

            if 'gene_name_primary' in keys:
                genename = rec.annotations['gene_name_primary']
            else:
                genename = None

            if 'organism' in keys:
                organism = rec.annotations['organism']
            else:
                organism = None

            # it's now very easy to parse additional information if we want it :) 
            ## Let me know if there's anything else we should pull in from Uniprot
            data[rec.id] = {'localizations': locs, 
                            'GO IDs': [r.split('GO:', 1)[1] for r in rec.dbxrefs if 'GO:' in r], 
                            'genename': genename, 
                            'organism': organism}
            
            
            if n % 1000 == 0:
                print(n, ' out of ', len(self.proteins), ' xml entries parsed')
             
        # note it is very easy to turn this data output into a csv from here, 
        # we could even have a filepath as an input and write directly from this function instead of returning the data dictionary
        
        return data
    
    
def get_string_interactions(stringIds, species):
    
    url = 'https://string-db.org/api/tsv/network'
    params = {'species': species, 'caller_identity': 'Inter-ViSTA'}

    if len(stringIds) >= 2000:
        n_chunks = int(np.ceil(len(stringIds)/2000))
        dfs = []

        for chunk in range(n_chunks):
            ps = stringIds[2000*chunk:2000*(chunk+1)]

            p = list2string(ps, '%0D') #each protein on a new line
            params['identifiers'] = p
            r = requests.post(url, data = params)
            _df = pd.read_csv(io.StringIO(r.text), sep = '\t', header = 0, index_col = None)

            dfs.append(_df)
            time.sleep(1) # courtesy waiting period to not overload the server

        df = pd.concat(dfs, axis = 0, join = 'outer')
   
    else:
        params['identifiers'] = list2string(stringIds, sep='%0D')
        r = requests.post(url, data = params)
        
        df = pd.read_csv(io.StringIO(r.text), sep = '\t', header = 0, index_col = None)
        
    return df
