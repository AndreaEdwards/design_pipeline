from src.pdb import PDBFromUniprot, CIFFromUniprot, LigandBindingSite, EditPDB
from src.pfam import PfamFromUniprot
from src.blast import BLAST
from src.genbank import GenBank
from src.hmm import HMM

def main():
	pass
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
#	pdb_codes = pdb_from_uniprot.get_pdb_id(queryText)
#	print pdb_codes

#	# Download .pdb using pdb_id
#	# THIS WORKS. DO NOT CHANGE
#	pdb_code = '3UUD'
#	pdb_dir = '/Users/Andrea/repositories/design_pipeline/pdbs/'
#	cif_getter = CIFFromUniprot()
#	cif_getter.fetch_mmCIF(pdb_code, pdb_dir)

#	# Edit the B-factor column of a pdb
#	# Currently this just prints the b-factors... needs work
#	pdb_editor = EditPDB('3UUD')
#	pdb_editor.b_factor()

# Structural pipeline:
# How many pdbs are there for each uniprot id? What is the difference between these sequences
# What ligand is bound? Extract out all residues within 5
# Get surface residues
# Find pockets

#	# Get pfam_id using uniprot_id. Can be many pfam_ids per target.
#	# THIS WORKS. DO NOT CHANGE
#	uniprot_id = 'P0A9Q1'
#	pfam_getter = PfamFromUniprot(uniprot_id)
#	id_list = pfam_getter.get_pfam_id()

#	# Fetch genbank record
#	# THIS WORKS. DO NOT CHANGE
#	genbank_getter = GenBank()
#	genbank_id = 1788201
#	genbank_getter.fetch_record(genbank_id)

#	# Fetch HMM
#	# THIS WORKS. DO NOT CHANGE
#	genbank_id = 1788201
#	hmm = HMM()
#	hmm.fetch_hmm(genbank_id)

#	# wwwBLAST using genbank_id for each target.
#	# Buggy (sleeps if no response from server). 
#	#genbank_record_number = 1787583
#	for genbank_id in genbank_ids:
#		blast = BLAST(genbank_id)
#		blast.blast_record()
#		blast.save_blast_xml()

			

main()