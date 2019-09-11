def run_functions(proteins, species):
	interactions = GetSTRINGInteractions()

	IDs_df = interactions.map_identifiers_string(list(proteins.values()), species)
	if IDs_df.empty: return 
	IDs = IDs_df.values.flatten()
	known_interactions = interactions.get_interactions(IDs, species)

	known_interactions['Original geneID_A'] = IDs_df.loc[known_interactions['ncbiTaxonId'].astype(str).str.cat(known_interactions['stringId_A'], sep='.').values, 'queryItem'].values
	known_interactions['Original geneID_B'] = IDs_df.loc[known_interactions['ncbiTaxonId'].astype(str).str.cat(known_interactions['stringId_B'], sep='.').values, 'queryItem'].values

	results = known_interactions.set_index(['Original geneID_A', 'Original geneID_B'])
	results = results.reset_index(level=['Original geneID_A', 'Original geneID_B'])
	results = results[results.columns[[0, 1, 4, 5, 12]]]

	return results

