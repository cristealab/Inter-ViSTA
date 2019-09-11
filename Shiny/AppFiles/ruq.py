def run_uniprot_query(proteins):
	uniprot = UniprotQuery(proteins).parse_xml()
	results = pd.DataFrame(uniprot).transpose()
	return uniprot