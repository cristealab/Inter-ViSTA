def run_uniprot_query(proteins):
    uniprot = UniprotQuery()
    entries_xml = uniprot.retrieve_entries(proteins)
    entries = uniprot.parse_xml(entries_xml)
    
    return entries

def run_functions(proteins, species):
    uniprot = UniprotQuery()
    IDs_df = uniprot.map_proteins(list(proteins.values()), to_='STRING_ID').set_index('To')
  
    if IDs_df.empty: 
      return 
    
    known_interactions = get_string_interactions(IDs_df.index.values.tolist(), species = species)

    known_interactions['Original geneID_A'] = IDs_df.loc[known_interactions['ncbiTaxonId'].astype(str).str.cat(known_interactions['stringId_A'], sep='.').values, 'From'].values
    known_interactions['Original geneID_B'] = IDs_df.loc[known_interactions['ncbiTaxonId'].astype(str).str.cat(known_interactions['stringId_B'], sep='.').values, 'From'].values

    return known_interactions[['Original geneID_A', 'Original geneID_B', 'preferredName_A', 'preferredName_B', 'escore']]
