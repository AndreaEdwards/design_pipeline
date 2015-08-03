--Run with command??

CREATE TABLE protein_information (
	enzyme text NOT NULL,
	genbank text,
	uniprot text,
	pdb_id text,
);
COPY protein_information (
	enzyme, 
	genbank,
	uniprot,
	pdb_id)
FROM '/Users/andrea/repositories/design_pipeline/database.csv' DELIMITER ',' CSV;

--How do I make sure that pdb_id is an array? probably specific datatype.