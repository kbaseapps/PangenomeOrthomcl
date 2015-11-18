#BEGIN_HEADER
import subprocess
import MySQLdb
import os
import shutil
import traceback
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
    def add_lines(self, log, lines):
        for line in lines:
            if len(line) > 0:
                log += "|" + line + "\n"
        return log
    
    def add(self, log, process):
        process_out = process.communicate()
        output = process_out[0]
        if output is not None and len(output) > 0:
            log += "Output:\n"
            log = self.add_lines(log, output.split("\n"))
        errors = process_out[1]
        if errors is not None and len(errors) > 0:
            log += "Errors:\n"
            log = self.add_lines(log, errors.split("\n"))
        return log
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.scratch = config['scratch']
        if os.path.exists(self.scratch):
            shutil.rmtree(self.scratch)
        os.makedirs(self.scratch)
        self.workspaceURL = config['workspace-url']
        #END_CONSTRUCTOR
        pass

    def build_pangenome_with_orthomcl(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN build_pangenome_with_orthomcl
        log = ""
        try:
            log += "Starting mysql service\n"
            log = self.add(log, subprocess.Popen(["service", "mysql", "start"], cwd=self.scratch, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
            #######################################################
            log += "Preparing database\n"
            db = MySQLdb.connect(host="localhost", user="root", passwd="12345");
            cur = db.cursor()
            cur.execute("DROP DATABASE IF EXISTS orthomcl")
            cur.execute("CREATE DATABASE orthomcl")
            cur.close()
            db.close()
            #######################################################
            log += "Preparing orthomcl config file\n"
            orthomcl_cfg = self.scratch + '/orthomcl.cfg'
            f = open(orthomcl_cfg, 'w')
            f.write("dbVendor=mysql\n");
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
            #######################################################
            log += "Running orthomclInstallSchema\n"
            log = self.add(log, subprocess.Popen(["perl", plbin + "/orthomclInstallSchema", orthomcl_cfg], cwd=self.scratch, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
            #######################################################
            log += "Loading GenomeSet object from workspace\n"
            token = ctx['token']
            ws = workspaceService(self.workspaceURL, token=token)
            genomeset = ws.get_objects([{'ref':params["intput_genomeset_ref"]}])[0]['data']
            #######################################################
            log += "Preparing genome refs\n"
            genome_refs = []
            for param_key in genomeset["elements"]:
                genome_refs.append(genomeset["elements"][param_key]["ref"])
            if len(genome_refs) > 10:
                raise ValueError('Number of genomes exceeds 10, which is too many for all-against-all blastp')
            feature_info = {}
            compliant_fasta_dir = self.scratch + "/compliantFasta"
            os.makedirs(compliant_fasta_dir)
            for genome_pos, genome_ref in enumerate(genome_refs):
                #######################################################
                log += "Loading Genome object from workspace for ref [" + genome_ref + "]\n"
                genome = ws.get_objects([{'ref' : genome_ref}])[0]['data']
                #######################################################
                log += "Preparing fasta file for ref [" + genome_ref + "]\n"
                genome_id = str(genome_pos + 1)
                records = []
                for feature_pos, feature in enumerate(genome['features']):
                    if feature['type'] == 'CDS' and 'protein_translation' in feature:
                        sequence = feature['protein_translation']
                        id = str(feature_pos + 1)
                        record = SeqRecord(Seq(sequence), id=id, description="")
                        records.append(record)
                        feature_info[genome_id + "|" + id] = {"fid": feature["id"], "gid": genome_id, "gref": genome_ref}
                fasta_file = self.scratch + "/" + genome_id + ".fasta"
                SeqIO.write(records, fasta_file, "fasta")
                #######################################################
                log += "Running orthomclAdjustFasta for ref [" + genome_ref + "]\n"
                log = self.add(log, subprocess.Popen(["perl", plbin + "/orthomclAdjustFasta", genome_id, fasta_file, "1"], cwd=compliant_fasta_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
            #######################################################
            log += "Running orthomclFilterFasta\n"
            log = self.add(log, subprocess.Popen(["perl", plbin + "/orthomclFilterFasta", compliant_fasta_dir, "50", "10"], cwd=self.scratch, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
            #######################################################
            log += "Running formatdb\n"
            protdb = "goodProteins.fasta"  # created by orthomclFilterFasta
            log = self.add(log, subprocess.Popen(["formatdb", "-i", protdb], cwd=self.scratch, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
            #######################################################
            log += "Running blastp\n"
            blastp_args = ["blastall", "-p", "blastp", "-d", protdb, "-i", protdb, "-F", "m S", "-v", "100000", "-b", "100000", "-e", "1e-5", "-m", "8", "-a", "1"]
            blast_output = self.scratch + "/blastres.txt"
            with open(blast_output, "w") as outfile:
                log = self.add(log, subprocess.Popen(blastp_args, cwd=self.scratch, stdout=outfile, stderr=subprocess.PIPE))
            #######################################################
            log += "Running orthomclBlastParser\n"
            sim_seq_file = self.scratch + "/similarSequences.txt"
            with open(sim_seq_file, "w") as outfile:
                log = self.add(log, subprocess.Popen(["perl", plbin + "/orthomclBlastParser", blast_output, compliant_fasta_dir], cwd=self.scratch, stdout=outfile, stderr=subprocess.PIPE))
            #######################################################
            log += "Running orthomclLoadBlast\n"
            log = self.add(log, subprocess.Popen(["perl", plbin + "/orthomclLoadBlast", orthomcl_cfg, sim_seq_file], cwd=self.scratch, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
            #######################################################
            log += "Running orthomclPairs\n"
            orthomcl_pairs_file = self.scratch + "/orthomcl_pairs.log"
            log = self.add(log, subprocess.Popen(["perl", plbin + "/orthomclPairs", orthomcl_cfg, orthomcl_pairs_file, "cleanup=no"], cwd=self.scratch, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
            #######################################################
            log += "Running orthomclDumpPairsFiles\n"
            log = self.add(log, subprocess.Popen(["perl", plbin + "/orthomclDumpPairsFiles", orthomcl_cfg], cwd=self.scratch, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
            #######################################################
            log += "Running mcl\n"
            mcl_output_file = self.scratch + "/mclOutput"
            log = self.add(log, subprocess.Popen(["mcl", "mclInput", "--abc", "-I", "1.5", "-o", mcl_output_file], cwd=self.scratch, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
            #######################################################
            log += "Running orthomclMclToGroups\n"
            groups_file = self.scratch + "/groups.txt"
            with open(groups_file, "w") as outfile, open(mcl_output_file, "r") as infile:
                log = self.add(log, subprocess.Popen(["perl", plbin + "/orthomclMclToGroups", "grp", "1000"], cwd=self.scratch, stdin=infile, stdout=outfile, stderr=subprocess.PIPE))
            #######################################################
            log += "Parsing groups file\n"
            prots = "";
            with open(groups_file, "r") as infile:
                for line in infile.readlines():
                    words = line.rstrip().split(" ")
                    for id in words[1:]:
                        feature = feature_info[id]
                        prots += id + "(" + feature["fid"] + ") "
            raise ValueError(prots)
        except Exception, err:
            log += traceback.format_exc() + "\n"
            raise ValueError(log)
        returnVal = {}
        #END build_pangenome_with_orthomcl

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method build_pangenome_with_orthomcl return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
