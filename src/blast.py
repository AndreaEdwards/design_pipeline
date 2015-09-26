import os
from Bio.Blast import NCBIWWW, NCBIXML
from src.util import FileHandlers
import subprocess

PATH_TO_PSIBLAST = '/Users/andrea/ncbi-blast-2.2.31+/bin'

PATH_TO_BLAST_DATABASE = '/Users/andrea/ncbi-blast-2.2.31+/bin/nr'

PSIBLAST_OPTIONS = {
    'db' : PATH_TO_BLAST_DATABASE ,

#    'num_threads' : 40 ,    # messes up my PBS submission
    'num_iterations' : 2 ,
    'pseudocount' : 2 ,

    'evalue' : 1 ,
    'inclusion_ethresh' : .001 ,
    'comp_based_stats' : 1 ,

    'outfmt' : 7 ,
    'num_descriptions' : 3000 ,
    'num_alignments' : 300 ,

    # still keeping output for potential features later
    'out' : lambda x : x + '.pb' ,
    'out_ascii_pssm' : lambda x : x + '.pssm' ,
    'out_pssm' : lambda x : x + '.cp' ,
    'export_search_strategy' : lambda x : x + '.ss'
    }

PROTEIN_LETTERS = 'ACDEFGHIKLMNPQRSTVWY'

# helper for making bash scripts
def create_executable_str( executable , args = [] , options = {} , out_filename = '' , args_first = False , extra_options = '' , append = False ):
    """
    Returns a str for a commandline call to  <executable>  (input as a string)
    using a list og  <args>  and a dict of  <options>  (keys as options, values
    as option values, empty str for flags)
    
    Optionally write stout to  <out_filename>  (only if specified)
    Optionally provide the arguments to the method first  (<args_first>)
    
    Add additional text  <extra_options>  to the call
    """
    if not isinstance( executable , str ):
        raise IOError( 'variable  <executable>  must be a string!' )
    if not isinstance( args , list ) and not isinstance( args , tuple ):
        raise IOError( 'variable  <args>  must be a list! (or otherwise iterable)' )
    if not isinstance( options , dict ):
        raise IOError( 'variable  <executable>  must be a string!' )
        
    # setup the options and args
    options = ''.join( [' -' + str( i ) + ( ' ' + str( options[i] ) )*bool( options[i] ) for i in options.keys()] )
    args = ''.join( [' ' + str( i ) for i in args] )
    
    # choose general format
    if args_first:
        perform = executable + args + options + extra_options
    else:
        perform = executable + options + args + extra_options
    
    # optionally write stdout
    if out_filename:
        if append:
            perform += ' >> ' + out_filename
        else:
            perform += ' > ' + out_filename

    return perform

# runs a commandline, usually combined with create_executable_str above
def run_local_commandline( command , collect_stdout = False ):
    """
    Runs the explicit  <command>  by performing a system call with subprocess
    """
    # get the output
    print '\n'+ '='*80 + '\nPerforming system call:\n' + command + '\n' + '='*80 +'\n'
    if not collect_stdout:
        subprocess.call( command , shell = True )    # just the command, no output piping

    # older call that pipes the output into Python for manipulation
    else:
        #stdout = subprocess.Popen( command , shell = True , stdout = subprocess.PIPE , stdin = subprocess.PIPE , stderr = subprocess.STDOUT ).communicate()[0].strip()
        # shell = True, do NOT .split the command, = False, DO .split the command
        stdout = subprocess.Popen( command , shell = True , stdout = subprocess.PIPE , stdin = subprocess.PIPE , stderr = subprocess.STDOUT ).communicate()[0]#.strip()    # for consistency
        return stdout

# local
def run_psiblast( sequence_filename , run = True ):
    """
    Runs PSIBLAST on  <sequence_filename>  using the default options in
    PSIBLAST_OPTIONS and returns the relevant output file: "out_ascii_pssm"
    """
    root_filename = os.path.abspath( sequence_filename ).rstrip( '.fa' )
    
    # collect the options, set the input, derive the output filenames
    psiblast_options = {}
    psiblast_options.update( PSIBLAST_OPTIONS )
    psiblast_options['query'] = sequence_filename
    for i in psiblast_options.keys():
        if '__call__' in dir( psiblast_options[i] ):
            psiblast_options[i] = psiblast_options[i]( root_filename )

    for i in psiblast_options.keys():
        if isinstance( psiblast_options[i] , str ) and os.path.isfile( psiblast_options[i] ):
            psiblast_options[i] = os.path.abspath( psiblast_options[i] )
    
    command = create_executable_str( PATH_TO_PSIBLAST , args = [] , options = psiblast_options )

    if run:
        run_local_commandline( command )
    
        # the only output we need
        return psiblast_options['out_ascii_pssm']
    else:
        # just send the command
        return command , psiblast_options['out_ascii_pssm']


# simple method, scan for empty/non-existent pssm + if "no hits" were found
def check_psiblast_output( psiblast_pssm , psiblast_output = None , failed_output_str = 'No hits found' ):
    # well, if the pssm file is NOT empty, things are good
    # if it is empty, check the psiblast_output for "No hits found"
    not_empty = None
    success = True
    
    if os.path.isfile( psiblast_pssm ):
        f = open( psiblast_pssm , 'r' )
        not_empty = bool( f.read().strip() )
        f.close()
    if not not_empty:
        # pssm does not exist OR was empty
        if not os.path.isfile( psiblast_output ):
            # both files empty or missing
            success = False
        else:
            f = open( psiblast_output , 'r' )
            lines = f.read()
            f.close()
            
            if failed_output_str.lower() in lines.lower():
                not_empty = False
                # but successful, just empty
            else:
                # major problems, both failed to generate
                # OR empty but also not supposed to be, rerun
                success = False
            
    return success , not_empty

# modified by njc, hybrid method
# now robust to versions, based on separator rather than anticipated structure
def extract_pssm_from_psiblast_pssm( pssm_filename , aa_line_shift = -4 , columns = len( PROTEIN_LETTERS ) ):
    # most args are ignored---left for backward compatibility.

    # fields can be separated by spaces and perhaps a '-' sign, or in some cases just a '-' sign.
    #
    # note: this re captures, so we will be given each field's separator. we need this to restore
    # '-' signs.
    split_line = re.compile( '( +-?|-)' )

    f = open( pssm_filename , 'r' )
    lines = [i.rstrip( '\n' ) for i in f.xreadlines()]
    f.close()

    split_lines = []
    for line in lines:
        ff = split_line.split( line )
        if not ff[0]:
            ff.pop( 0 )    # re splits include a dummy empty field if the first field
                           # begins with a separator. (so drop it)
        ff = [sep[-1] + dat for sep , dat in zip( ff[::2] , ff[1::2] )]    # restore the separators to each field.
        # only care about "-" character if its there
        split_lines.append( ff )

    # HACK: assume lines with data have the largest number of columns.
    most_columns = max( [len( ff ) for ff in split_lines] )

    # line listing amino acids has four fewer columns: it is missing position, query, information,
    # and relative weight fields. there should be one and only one such line.
    aa_line = [ff for ff in split_lines if len( ff ) == most_columns + aa_line_shift]
    assert len( aa_line ) == 1
    aa_line = aa_line[0]
    
    # build a map from column index to aa, then run some sanity checks. note the line lists each aa twice.
    x2aa = dict( [(x , aa.strip()) for x , aa in enumerate( aa_line )] )
    num_aa = len( x2aa )/2
    assert num_aa*2 == len( x2aa )
    assert num_aa == columns    # sanity check
    assert not [x for x in xrange( num_aa ) if not x2aa[x] == x2aa[x + num_aa]]

    pssm_dict = {}
    for ff in split_lines:
        if not len( ff ) == most_columns: continue
        pos = int( ff[0] )
        pssm_dict[pos] = {
            'position':                pos ,
            'query identity':          ff[1] ,
            'log-likelihood':          dict( [(x2aa[x] , int(v)) for x , v in enumerate( ff[2:2 + num_aa] )] ) ,
            'approximate frequencies': dict( [(x2aa[x] , float(v)/100) for x , v in enumerate( ff[2+num_aa:2 + 2*num_aa] )] ) ,
            'information content':     float( ff[-2] ) ,
            '?':                       float( ff[-1] ) ,
        }
    #DEBUG pssm_dict['XTRA_AALINE'] = aaline
    return pssm_dict


class BLAST:
	def __init__(self, genbank_record_number):
		self.genbank_record_number = genbank_record_number
		self.dir_path = os.getcwd() + '/blast_results'
		self.xml_result = ''
		self._mkdir() 

	def blast_record(self):
		print "BLASTing record number %d ..." % int(self.genbank_record_number)
		result_handle = NCBIWWW.qblast("blastp", "nr", self.genbank_record_number)
		print "extracting result..."
		self.xml_result = result_handle.read()
		result_handle.close()
		return self.xml_result
	
	def save_blast_xml(self):
		print "Saving BLAST result (.xml) for record number %d" % int(self.genbank_record_number)
		file_name = str(self.genbank_record_number) + "_blast.xml"
		file_path = self.dir_path + '/' + file_name
		save_file = open(file_path, "w")
		save_file.write(self.xml_result)
		save_file.close()

	def parse_xml(self, file_path):
		result_handle = open(file_path)
		blast_record = NCBIXML.read(result_handle)
		result_handle.close()
		return blast_record

	def _mkdir(self):
		file_handlers = FileHandlers()
		file_handlers.make_results_folder(self.dir_path.split('/')[-1])

class localPSIBLAST:
	def __init__(self):
		pass

def main():
##	Example of how to use the xml output:
#	file_paths = file_handlers.search_directory()
#	xml_files = file_handlers.find_files(file_paths, 'xml')
#	for xml_file in xml_files:
#		blast_record = blast.parse_xml(xml_file)
#		for description in blast_record.description:
#			print description.title
#		for alignment in blast_record.alignments:
#			#print alignment
#			for hsp in alignment.hsps:
#				#print('**** Alignment ****')
#				#print('sequence:', alignment.title)
#				#print('length:', alignment.length)
#				#print('e value', hsp.expect)
#				#print(hsp.query[0:75] + '...')
#				#print(hsp.match[0:75] + '...')
#				#print(hsp.sbjct[0:75] + '...')
#				print(hsp.identities)
	run_psiblast( '/Users/andrea/repositories/design_pipeline/src/4B2N.fa' )

main()
