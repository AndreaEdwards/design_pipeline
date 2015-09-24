from src.pdb import (PDBFromUniprot, 
	CIFFromUniprot, LigandBindingSite, EditPDB, 
	SurfaceResidues, Rosetta, MutantListMaker,
	DDGMonomer)
from src.pfam import PfamFromUniprot
#from src.blast import BLAST
from src.genbank import GenBank
from src.hmm import HMM
from src.stats import Correlation

def main():

#	steele_uniprot_codes = ['P07001', 'P0AB67', 'P27306', 'P0A7B3', 'P0ACS2', 'P0A9E2', 'P0A9B2', 'P0A6T1', 'P0AFG8', 'P06959', 'P0A9P0', 'P0AFG3', 'P0AFG6', 'P61889', 'P33940', 'P0AC53', 'P00350']

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
	uniprot_ids = ['P0AB67', 'P27306', 'P0A7B3', 'P0A9E2', 'P0AFG8', 'P0AFG3', 'P0AFG6', 'P33940', 'P0AC53', 'H6N162']
	genbank_ids = ['81171066', '11182439', '67470903', '71162387', '84027826', '84027822', '84027824', '2506692', '81175321', '170180374']
	
	pdb_codes = ['2XUVABCD', '1Y00AB']
	pops_issue = []
	pocket_finder_issue = ['1W36DBCY', '2ZHG']
	ligand_issue = ['2ZHG']
	cleaned_file_mapping = {'1W36' : '1W36DBCY', '4B2N' : '4B2NA', '1U60' : '1U60AB', '2XUV' : '2XUVABCD', '1Y00' : '1Y00AB'}
	done = ['1X15', '4N72', '4JDR', '2CMD', '2ZYA', '3NBU', '4TWZ', '1SRU', '1S7C', '4B2NA', '2J1N', '1U60AB', '1YAC']
	# Download .pdb using pdb_id
	# THIS WORKS. DO NOT CHANGE
#	pdb_getter = PDBFromUniprot()
	print "===================================================================="
	print "Entering BioVerse Design Pipeline"
	print "===================================================================="
#
	for pdb_code in pdb_codes:
		print "\n\nIdentifying residues of interest for %s" % pdb_code
		# Download .pdb using pdb_id
		# THIS WORKS. DO NOT CHANGE
#		pdb_getter.fetch_pdb(pdb_code)
		pdb_editor = EditPDB(pdb_code)
#
## Structural pipeline:
#
##	# Get surface residues
##	# THIS WORKS. DO NOT CHANGE
		sr_getter = SurfaceResidues(pdb_code)
		sr_getter.write_resi_sasa_output()
		sr_getter.write_frac_sasa_output()
		pdb_editor.edit_bfactor_sasa()
		sr_getter.write_surface_resi_output(0.3)
		pdb_editor.edit_bfactor_surface_residues()
#
#
## How many pdbs are there for each uniprot id? What is the difference between these sequences
## What ligand is bound? Extract out all residues within 5
#
#	# Get ligand
#	# THIS WORKS. DO NOT CHANGE
		ligand = LigandBindingSite(pdb_code)
		ligand.get_residues_within_5A()
		ligand.write_residue_output()
		pdb_editor.edit_bfactor_ligand_binding_pocket()
#
#
## 	# Find pockets 
##	# THIS WORKS. DO NOT CHANGE
		rosetta = Rosetta(pdb_code)
		rosetta.find_pockets()
		pdb_editor.edit_bfactor_pocket_residues()
#
#
##	# Make 'mutants_list' file for ddg_monomer
##	# THIS WORKS. DO NOT CHANGE
		ListMaker = MutantListMaker(pdb_code)
#		ListMaker.generate_mutant_list(pocketres=True, lpocket=True, SurfRes=True)
		ListMaker.generate_mutant_list(pocketres=True, lpocket=True)
#
		ddgMonomer = DDGMonomer(pdb_code)
		ddgMonomer.get_targets(5.5)
#
##	# Linear regression
##	correlation = Correlation()
##	correlation.linregress('4PDJ_fracsasa.txt', 'test_pulled_bfactors.txt')
#
##	# Get pfam_id using uniprot_id. Can be many pfam_ids per target.
##	# THIS WORKS. DO NOT CHANGE
##	uniprot_id = 'P0A9Q1'
##	pfam_getter = PfamFromUniprot(uniprot_id)
##	id_list = pfam_getter.get_pfam_id()
#
##	# Fetch HMM
##	# THIS WORKS. DO NOT CHANGE
#	for genbank_id in genbank_ids:
##	genbank_id = 1788201
#		hmm = HMM()
#		hmm.fetch_hmm(genbank_id)
#
##	# Fetch genbank record
##	# THIS WORKS. DO NOT CHANGE
#		genbank_getter = GenBank()
##	genbank_id = 1788201
#		genbank_getter.fetch_record(genbank_id)
#
##	# wwwBLAST using genbank_id for each target.
##	# Buggy (sleeps if no response from server). 
##	#genbank_record_number = 1787583
##	for genbank_id in genbank_ids:
##		blast = BLAST(genbank_id)
##		blast.blast_record()
##		blast.save_blast_xml()
#
##	# Implement linear regression model. linear algebra... numpy?

			

main()