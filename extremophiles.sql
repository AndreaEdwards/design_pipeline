CREATE TABLE protein_information (
	enzyme text NOT NULL,
	source_organism text,
	optimal_temperature_C int,
	optimal_ph float(8), --float(8) means that the values stored here have a minimum acceptable precision of 8 binary digits.
	optimal_salt_M float(8),
	genbank_accession_number text,
	pdb_id text,
	melting_temp_C float(8),
	reference text,
	temperature_type text,
	pH_type text,
	salt_type text,
	mesophilic_counterpart_pdb text,
	stoichiometry text
);
COPY protein_information (
	enzyme, 
	source_organism, 
	optimal_temperature_C, 
	optimal_ph, 
	optimal_salt_M,
	genbank_accession_number,
	pdb_id,
	melting_temp_C,
	reference,
	temperature_type,
	pH_type,
	salt_type,
	mesophilic_counterpart_pdb,
	stoichiometry)
FROM '/Users/andrea/repositories/protein_engineering/extremophiles.csv' DELIMITER ',' CSV;

--INSERT INTO protein_information (
--	enzyme, 
--	source_organism, 
--	optimal_temperature_C, 
--	optimal_ph, 
--	optimal_salt_M,
--	genbank_accession_number,
--	pdb_id,
--	melting_temp_C,
--	reference,
--	temperature_type,
--	pH_type,
--	salt_type,
--	mesophilic_counterpart_pdb,
--	stoichiometry)
--VALUES 
--	('L-Threonine Dehydrogenase','Pyrococcus horikoshii',70,10,NULL,'2DFV_C','2DFV',NULL,'','Thermophile','Alkaliphile','','4ILK_A','Monomer')
