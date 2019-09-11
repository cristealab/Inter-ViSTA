
import requests
import pandas as pd
import numpy as np
import io
from Bio import SeqIO

class UniprotQuery:
    def __init__(self, proteins):
        
        query = ''

        for item in proteins:
            try:
                query = query + item + ' '

            except TypeError: # exception to deal with NaNs
                pass

        self.proteins = proteins
        self.query = query
        
        url = 'https://www.uniprot.org/uploadlists/'

        params = {'query': self.query, 
                  'format': 'xml',
                  'from': 'ACC+ID',
                  'to': 'ACC'}

        r = requests.post(url, data = params)
        
        self.response = r
        
    def parse_xml(self):
        
        data = {}

        for n, rec in enumerate(SeqIO.parse(io.StringIO(self.response.text), 'uniprot-xml'), 1):

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
                
# function would be used in this way:
# background = pd.read_csv(r'Z:\2 Programming\4 Inter-ViSTA\testing files\Background Gene Lists\background_gene_list_human_fibroblast.txt', sep='\t').accession.values
# data = UniprotQuery(background).parse_xml()

# then write the data to a csv or whatever