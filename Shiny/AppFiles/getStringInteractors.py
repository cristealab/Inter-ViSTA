# ADAPTED FROM 
# Michelle Kennedy
# 1.29.2019

# Below are the main functions needed to get the string interactions. 

import requests

# for outputting the response from STRING to a pandas data frame
import io 
import pandas as pd
import numpy as np
import time

class GetSTRINGInteractions:
    def __init__(self):
        pass
    
    def to_query_string(self, mylist, sep): #can also accept arrays
        '''convert a list to a string that can be used as a query string in an http post request'''
        
        l = ''

        for item in mylist:
            try:
                l = l + str(item) + sep
            
            except TypeError: # exception to deal with NaNs in mylist
                pass
        
        return l
    
    def map_identifiers_string(self, proteins, species):
        '''Use STRING's API to retrive the corresponding string identifiers for each protein.
        I highly reccommend using this function prior to querying for interactions. 
        It significantly decreases the request response time'''
        
        # STRING will only let you query 2000 proteins at a time, otherwise you get an error message back
        
        if len(proteins) >= 2000:
            n_chunks = int(np.ceil(len(proteins)/2000))
            dfs = []
            
            for chunk in range(n_chunks):
                ps = proteins[2000*chunk:2000*(chunk+1)]
                
                p = self.to_query_string(ps, '%0D') #each protein on a new line

                url = 'https://string-db.org/api/tsv/get_string_ids'
                params = {'identifiers': p, 'species':species, 'echo_query': 1, 'caller_identity': 'Princeton_University'}

                r = requests.post(url, data = params)
                _df = pd.read_csv(io.StringIO(r.text), sep = '\t', header = 0, index_col = None)
                
                dfs.append(_df)
                time.sleep(1)
                
            df = pd.concat(dfs, axis = 0, join = 'outer')
   
        else:
            ps = proteins
        
            p = self.to_query_string(ps, '%0D') #each protein on a new line

            url = 'https://string-db.org/api/tsv/get_string_ids'
            params = {'identifiers': p, 'species':species, 'echo_query': 1, 'caller_identity': 'Princeton_University'}

            r = requests.post(url, data = params)
            df = pd.read_csv(io.StringIO(r.text), sep = '\t', header = 0, index_col = None)
            
            
        df = df[['stringId', 'queryItem']].set_index('stringId')
        
        return df

    
    def get_interactions(self, IDs, species):
        
        # STRING will only let you query 2000 proteins at a time
        
        if len(IDs) > 2000:
            
            n_chunks = int(np.ceil(len(IDs)/2000))
            
            dfs = []
            
            for chunk in range(n_chunks):
                ID_list = IDs[2000*chunk:2000*(chunk+1)]
                
                p = self.to_query_string(ID_list, '%0D') #each ID on a new line

                url = 'https://string-db.org/api/tsv/network'
                params = {'identifiers': p, 'species':species, 'caller_identity': 'Princeton_University'}

                r = requests.post(url, data = params)
                _df = pd.read_csv(io.StringIO(r.text), sep = '\t', header = 0, index_col = None)
                dfs.append(_df)
                time.sleep(1)
            
            df = pd.concat(dfs, axis = 0, join = 'outer')
                
        else:
            ID_list = IDs
        
            p = self.to_query_string(ID_list, '%0D') #each ID on a new line

            url = 'https://string-db.org/api/tsv/network'
            params = {'identifiers': p, 'species':species, 'caller_identity': 'Princeton_University'}

            r = requests.post(url, data = params)
            df = pd.read_csv(io.StringIO(r.text), sep = '\t', header = 0, index_col = None)
        
        return df