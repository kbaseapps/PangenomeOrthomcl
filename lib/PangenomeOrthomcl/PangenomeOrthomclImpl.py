#BEGIN_HEADER
import subprocess
import MySQLdb
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from biokbase.workspace.client import Workspace as workspaceService
#END_HEADER


class PangenomeOrthomcl:
    '''
    Module Name:
    PangenomeOrthomcl

    Module Description:
    A KBase module: PangenomeOrthomcl
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.scratch = config['scratch']
        self.workspaceURL = config['workspace-url']
        #END_CONSTRUCTOR
        pass

    def build_pangenome_with_orthomcl(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN build_pangenome_with_orthomcl
        process = subprocess.Popen(["service", "mysql", "start"], stdout=subprocess.PIPE)
        output = process.communicate()[0]
        print("mysql: " + output)
        db = MySQLdb.connect(host="localhost", user="root", passwd="12345");
        cur = db.cursor()
        cur.execute("DROP DATABASE IF EXISTS orthomcl")
        cur.execute("CREATE DATABASE orthomcl")
        cur.close()
        db.close()
        orthomcl_cfg = self.scratch + '/orthomcl.cfg'
        f = open(orthomcl_cfg, 'w')
        f.write("dbConnectString=dbi:mysql:orthomcl:mysql_local_infile=1:localhost:3306\n")
        f.write("dbLogin=root\n")
        f.write("dbPassword=12345\n")
        f.write("similarSequencesTable=SimilarSequences\n")
        f.write("orthologTable=Ortholog\n")
        f.write("inParalogTable=InParalog\n")
        f.write("coOrthologTable=CoOrtholog\n")
        f.write("interTaxonMatchView=InterTaxonMatch\n")
        f.write("percentMatchCutoff=50\n")
        f.write("evalueExponentCutoff=-5\n")
        f.write("oracleIndexTblSpc=NONE\n")
        f.close()
        plbin = "/kb/deployment/plbin"
        process = subprocess.Popen(["perl", plbin + "/orthomclInstallSchema", orthomcl_cfg], stdout=subprocess.PIPE)
        output = process.communicate()[0]
        token = ctx['token']
        ws = workspaceService(self.workspaceURL, token=token)
        genomeset = ws.get_objects([{'ref':params["intput_genomeset_ref"]}])[0]['data']
        feature_info = {}
        records = []
        feature_index = 0
        for param_key in genomeset["elements"]:
            genome_ref = genomeset["elements"][param_key]["ref"]
            genome = ws.get_objects([{'ref' : genome_ref}])[0]['data']
            for feature in genome['features']:
                if feature['type'] == 'CDS' and 'protein_translation' in feature:
                    sequence = feature['protein_translation']
                    feature_index += 1
                    id = str(feature_index)
                    feature_info[id] = [genome_ref, feature["id"]]
                    record = SeqRecord(Seq(sequence), id=id, description="#" + genome_ref + ", " + feature["id"])
                    records.append(record)
        fasta_file = self.scratch + "/blast_temp_input.fa"
        SeqIO.write(records, fasta_file, "fasta")
        returnVal = {}
        #END build_pangenome_with_orthomcl

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method build_pangenome_with_orthomcl return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
