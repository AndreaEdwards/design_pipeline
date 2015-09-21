from src.pdb import (PDBFromUniprot, 
	CIFFromUniprot, LigandBindingSite, EditPDB, 
	SurfaceResidues, Rosetta, MutantListMaker,
	DDGMonomer)
from src.pfam import PfamFromUniprot
from src.blast import BLAST
from src.genbank import GenBank
from src.hmm import HMM
from src.stats import Correlation

def main():
	
#	# Get uniprot_id for target in database. Only one uniprot_id for each target.
#	# This needs to be changed for postgress
#	genbank_ids = []
#	uniprot_ids = []
#
#	Data = open('/Users/Andrea/repositories/design_pipeline/src/database/database.csv')
#	data = Data.read().split('\r')
#	Data.close
#	for line in data:
#		if line.startswith("#"):
#			pass
#		else:
#			genbank_id, uniprot_id = line.split(',')[1], line.split(',')[2]
#			genbank_ids.append(genbank_id)
#			uniprot_ids.append(uniprot_id)

#	# Get pdb_id for each structure in PDB with respective uniprot_id. Can be many pdb's per target. 
#	# THIS WORKS. DO NOT CHANGE. 
#	uniprot_id = 'P0A9Q1'
#	queryText = ("<orgPdbQuery>" +
#				 	"<queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>" +
#						"<description>Simple query for a list of Uniprot Accession IDs: " + uniprot_id + " </description>" +
#						"<accessionIdList>" + uniprot_id + "</accessionIdList>" +
#					"</orgPdbQuery>")
#	pdb_from_uniprot = PDBFromUniprot()
	
#	steele_uniprot_codes = ['P07001', 'P0AB67', 'P27306', 'P0A7B3', 'P0ACS2', 'P0A9E2', 'P0A9B2', 'P0A6T1', 'P0AFG8', 'P06959', 'P0A9P0', 'P0AFG3', 'P0AFG6', 'P61889', 'P33940', 'P0AC53', 'P00350']

#	steele_pdb_codes = ['1X15', '2ZHG', '1S7C', '3NBU', '4N72', '4JDR', '2CMD', '2ZYA']
#	pdb_codes = pdb_from_uniprot.get_pdb_id(queryText)
#	print pdb_codes

#	# Download .pdb using pdb_id
#	# THIS WORKS. DO NOT CHANGE
#	pdb_getter = PDBFromUniprot()
#	pdb_code = '4PDJ'
	#pdb_code = '2XGE'
#	pdb_getter.fetch_pdb(pdb_code)
#	cif_getter = CIFFromUniprot()
#	cif_getter.fetch_mmCIF(pdb_code, pdb_dir)

# Structural pipeline:
# How many pdbs are there for each uniprot id? What is the difference between these sequences
# What ligand is bound? Extract out all residues within 5

#	# Get ligand
#	# THIS WORKS. DO NOT CHANGE
	ligand = LigandBindingSite('4PDJ')
	ligand.get_residues_within_5A()
	ligand.write_residue_output()

#	# Get surface residues
#	# THIS WORKS. DO NOT CHANGE
	sr_getter = SurfaceResidues('4PDJ')
	sr_getter.write_resi_sasa_output()
	sr_getter.write_surface_resi_output(0.3)
#	sr_getter.write_frac_sasa_output()

# 	# Find pockets 
#	# THIS WORKS. DO NOT CHANGE
	rosetta = Rosetta('4PDJ')
	rosetta.find_pockets()

#	# Edit the B-factor column of a pdb
#	# THIS WORKS. Needs improvement with file handling
	pdb_editor = EditPDB('4PDJ')
	pdb_editor.edit_bfactor_sasa()
	pdb_editor.edit_bfactor_ligand_binding_pocket()
	pdb_editor.edit_bfactor_surface_residues()
#	pdb_editor.edit_bfactor_pockets()
	pdb_editor.edit_bfactor_pocket_residues()
#	pdb_editor.write_bfactor()

#	# Make 'mutants_list' file for ddg_monomer
#	# THIS WORKS. DO NOT CHANGE
	ListMaker = MutantListMaker('4PDJ')
	ListMaker.generate_mutant_list()

	ddgMonomer = DDGMonomer('4PDJ')
	ddgMonomer.get_targets(5.5)

#	# Linear regression
#	correlation = Correlation()
#	correlation.linregress('4PDJ_fracsasa.txt', 'test_pulled_bfactors.txt')

#	# Get pfam_id using uniprot_id. Can be many pfam_ids per target.
#	# THIS WORKS. DO NOT CHANGE
#	uniprot_id = 'P0A9Q1'
#	pfam_getter = PfamFromUniprot(uniprot_id)
#	id_list = pfam_getter.get_pfam_id()

#	# Fetch HMM
#	# THIS WORKS. DO NOT CHANGE
#	genbank_id = 1788201
#	hmm = HMM()
#	hmm.fetch_hmm(genbank_id)

#	# Fetch genbank record
#	# THIS WORKS. DO NOT CHANGE
#	genbank_getter = GenBank()
#	genbank_id = 1788201
#	genbank_getter.fetch_record(genbank_id)

#	# wwwBLAST using genbank_id for each target.
#	# Buggy (sleeps if no response from server). 
#	#genbank_record_number = 1787583
#	for genbank_id in genbank_ids:
#		blast = BLAST(genbank_id)
#		blast.blast_record()
#		blast.save_blast_xml()

#	# Implement linear regression model. linear algebra... numpy?

			

main()