from src.pdb import (PDBFromUniprot, 
	CIFFromUniprot, LigandBindingSite, EditPDB, 
	SurfaceResidues, Rosetta, MutantListMaker,
	DDGMonomer, PDBPreProcessor, HeaderParser, CIFFParser)
from src.pfam import PfamFromUniprot
from src.blast import localBLAST
from src.genbank import GenBank
from src.hmm import HMM
from src.stats import Correlation
from src.seq import MutationListGenerator, SequenceGetter

def main():

#	# Get pdb_id for each structure in PDB with respective uniprot_id. Can be many pdb's per target. 
#	# THIS WORKS. DO NOT CHANGE. 

# Ignore this for now. Need to build tool that allows selection of a pdb based on header file information.
# Need to build capability for psiBLAST.	
#	pdb_from_uniprot = PDBFromUniprot()
#	pdb_getter = PDBFromUniprot()
#	uniprot_pdb_mapping = {}
#	for uniprot_id in steele_uniprot_codes:
#		uniprot_pdb_mapping[uniprot_id] = []
#		pdb_codes = pdb_from_uniprot.get_pdb_id(queryText)
#		if pdb_codes != []:
#			for pdb_id in pdb_codes:
#				uniprot_pdb_mapping[uniprot_id].append(pdb_id)
#		else:
#			print "%s returned no pdb codes" % uniprot_id
#	print "steele's mapping: ",  uniprot_pdb_mapping
#	uniprot_id = 'P0A9Q1'
#	queryText = ("<orgPdbQuery>" +
#				 	"<queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>" +
#						"<description>Simple query for a list of Uniprot Accession IDs: " + uniprot_id + " </description>" +
#						"<accessionIdList>" + uniprot_id + "</accessionIdList>" +
#					"</orgPdbQuery>")
#	pdb_from_uniprot = PDBFromUniprot()

	# proteins with no E. coli structures
	uniprot_ids = ['P07001', 'P0AB67', 'P27306', 'P0A7B3', 'P0ACS2', 'P0A9E2', 'P0A9B2', 'P0A6T1', 'P0AFG8', 'P06959', 'P0A9P0', 'P0AFG3', 'P0AFG6', 'P61889', 'P33940', 'P0AC53', 'P00350', 'P0AB67', 'P27306', 'P0A7B3', 'P0A9E2', 'P0AFG8', 'P0AFG3', 'P0AFG6', 'P33940', 'P0AC53', 'H6N162']
	genbank_ids = ['81171066', '11182439', '67470903', '71162387', '84027826', '84027822', '84027824', '2506692', '81175321', '170180374']
	ronming = ['3TCF', '3TCG', '3TCH']
	gene_names = ['pntA']
	pdb_codes = ['2BHS']
	queue = ['...', '3BP8', '1LB2', '1DPS', '1BG8', '2XUV', '3UCS']
	#pdb_codes = ['2ZHG', '2CGP', '2H27']
	rerun = ['4JDR', '3NBU', '2CMD', '2ZYA', '4TWZ', '1SRU', '1JTS']
	#pdb_codes = ['1X15', '4N72', '4JDR', '2CMD', '2ZYA', '3NBU', '4TWZ', '1SRU', '1W36']
	dna_bound = ['3JRG', '3JRH', '2ZHG', '2CGP', '2H27']
	key_error = ['1XGE'] # line 1433 in _edit_bfactor KeyErro 103
	list_index_out_of_range = ['4TWZ', '2BUI', '2BHS'] # list index out of range line 1395
	ligan_is_on_chain_B = ['1X15']
	homodimer_is_inverted = ['3PFL', '2ZYA', '1NYT', '2EFT', '1DFG'] # targeted residues are on wrong chain or both chains. need to get residues for both chains (copy B to A)
	hold = ['4OX6']
	two_chains = ['1NYE', '3UCS', '2O97']
	gaps_in_structure = ['2G67'] # need to fix write_csv_output to be robust to gaps in seq. currently just subtracting the first aln_diff
	pops_issue = []
	pocket_finder_issue = ['1IHF', '1OR7', '1SRU']
	ligand_issue = ['2ZHG', '1IHF']
	cleaned_file_mapping = {'1W36' : '1W36DBCY', '1W36': '1W36D', '1W36': '1W36B', '1W36': '1W36C', '4B2N' : '4B2NA', '1U60' : '1U60AB', '2XUV' : '2XUVABCD', '1Y00' : '1Y00AB'}
	done = ['1X15', '4N72', '4JDR', '2CMD', '2ZYA', '3NBU', '4TWZ', '1SRU', '1S7C', '4B2NA', '2J1N', '1U60AB', '1YAC', '2XUVABCD', '1Y00AB', '1W36D', '1W36B', '1W36C', '3NR7', '2L15', '1A04', '2GQQ', '1LB2', '3I2Z', '3TCH', '3TCH', '1H16', '3TCH']
	preliminary = []
	final = ['4N72', '3NBU', '3TCH', '1DPS', '1NYE', '3NR7', '2GQQ', '3I2Z', '1BG8', '1YAC', '2CMD', '1GS5', '3N1C', '2R97', '4LFU', '4JDR', '3PFL', '1S7C', '2GFY', '2GFW']
	# Download .pdb using pdb_id
	# THIS WORKS. DO NOT CHANGE
#	pdb_getter = PDBFromUniprot()
	print "===================================================================="
	print "Entering bioverse Design Pipeline"
	print "===================================================================="
#
	for pdb_code in pdb_codes:
	#	print "\n\nIdentifying residues of interest for %s" % pdb_code
		# Download .pdb using pdb_id
		# THIS WORKS. DO NOT CHANGE
		preprocessor = PDBPreProcessor(pdb_code)
		preprocessor.preprocess_check()
		#pdb_getter = PDBFromUniprot()
		#pdb_getter.fetch_pdb(pdb_code)
		#cif_getter = CIFFromUniprot(pdb_code)
		#cif_getter.fetch_mmCIF()
#
## Structural pipeline:
##	# Parse PDB header
##	# THIS WORKS. DO NOT CHANGE
		header_parser = HeaderParser(pdb_code)
		header_parser.get_header_dict()
#
##	# CIF parser
##	# THIS WORKS. DO NOT CHANGE
		cif_parser = CIFFParser(pdb_code)
		chains, genes, organisms, pdb_sequences = cif_parser.get_gene_annotations()
		sequence_annotations = cif_parser.collate_sequence_annotations()
		print sequence_annotations
		cif_parser.write_fasta(sequence_annotations)

# what is the gene_name on each chain and what is the sequence on each chain. then save each sequence to different fasta file
# ( chain, gene_name, sequence )

#
##	# Sequence Analysis
##	# THIS WORKS. DO NOT CHANGE
		seq_getter = SequenceGetter(genes, organisms)
		dna_sequences = seq_getter.get_DNA_sequence()
		seq_getter.get_protein_sequence()
#
##	# Preprocessing
##	# THIS WORKS. DO NOT CHANGE
		number_of_structures = preprocessor.process()
		preprocessor.get_complex_info()
#
##	# Find start site for structure
##	# THIS WORKS. DO NOT CHANGE
		local_blast = localBLAST(pdb_code, sequence_annotations)
		aln_diff = local_blast.align_pdb_seqs()
		print aln_diff
#
##	# Instantiate PDBEditor
##	# THIS WORKS. DO NOT CHANGE
		pdb_editor = EditPDB(pdb_code, server_mode=False)
#
##	# Get surface residues
##	# THIS WORKS. DO NOT CHANGE
	#	sr_getter = SurfaceResidues(pdb_code, server_mode=False)
	#	sr_getter.write_resi_sasa_output()
	#	sr_getter.write_frac_sasa_output()
	#	pdb_editor.edit_bfactor_sasa()
	#	sr_getter.write_surface_resi_output(0.3)
	#	pdb_editor.edit_bfactor_surface_residues()
#
#	# Get ligand
#	# THIS WORKS. DO NOT CHANGE
	#	ligand = LigandBindingSite(pdb_code, chains)
	#	ligand.get_residues_within_5A()
	#	ligand.write_residue_output()
	#	pdb_editor.edit_bfactor_ligand_binding_pocket()
#
#
## 	# Find pockets 
##	# THIS WORKS. DO NOT CHANGE
	#	rosetta = Rosetta(pdb_code)
	#	rosetta.find_pockets()
		pdb_editor.edit_bfactor_pocket_residues() 
#
#
##	# Create Rosetta minimized structure for non-redundant monomers in the structure
##	# This will be used for ddG monomer
		preprocessor = PDBPreProcessor(pdb_code, target_finder_mode=False)
		preprocessor.process()
		number_of_structures, chains, minimum_chains = preprocessor.count_structures_in_asymmetric_unit()
#
#
##	# Make 'mutants_list' file for ddg_monomer
##	# THIS WORKS. DO NOT CHANGE
		ListMaker1 = MutantListMaker(pdb_code, chains)
#		ListMaker.generate_mutant_list(lpocket=True, SurfRes=True)
		ListMaker1.generate_mutant_list(pocketres=True, lpocket=True, SurfRes=True)
#		ListMaker.generate_mutant_list(pocketres=True, lpocket=True)
		ListMaker2 = MutantListMaker(pdb_code, minimum_chains, asymmetric_unit=True)
		ListMaker2.filter_mutant_list(minimum_chains)
#
#
##	# Calculate ddG
##	# THIS WORKS. DO NOT CHANGE
		ddgMonomer = DDGMonomer(pdb_code, minimum_chains)
		ddgMonomer.get_targets(5.5)
#
##	# Write output csv
##	# THIS WORKS. DO NOT CHANGE
		mutation_list_generator = MutationListGenerator(sequence_annotations, dna_sequences, aln_diff, pdb_code=pdb_code)
		mutation_list_generator.write_csv_output(aln_diff)


	

			

main()