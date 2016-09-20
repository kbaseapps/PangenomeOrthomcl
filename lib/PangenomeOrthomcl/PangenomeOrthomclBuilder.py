import subprocess
import MySQLdb
import os
import shutil
import json
import uuid
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from biokbase.workspace.client import Workspace as workspaceService
from GenomeAnnotationAPI.GenomeAnnotationAPIClient import GenomeAnnotationAPI

class PangenomeOrthomclBuilder:
    '''
    Module Name:
    PangenomeOrthomclBuilder
    '''

    def __init__(self, scratch, workspaceURL, params, token, provenance):
        self.scratch = scratch
        self.workspaceURL = workspaceURL
        self.params = params
        self.token = token
        self.provenance = provenance
        self.plbin = "/kb/deployment/plbin"
        self.log = ""
        self.ws = workspaceService(self.workspaceURL, token=self.token)


    def run(self):
        self.log_line("Input parameters: " + json.dumps(self.params))
        if os.path.exists(self.scratch):
            shutil.rmtree(self.scratch)
        os.makedirs(self.scratch)
        self.startup_mysql()
        self.prepare_mysql_db()
        orthomcl_cfg = self.prepare_othomcl_config()
        self.orthomcl_install_schema(orthomcl_cfg)
        genomeset = self.load_genomeset_object()
        genome_refs = self.prepare_genome_refs(genomeset)
        compliant_fasta_dir = self.scratch + "/compliantFasta"
        feature_info = self.load_genome_features_prepare_fasta(genome_refs, compliant_fasta_dir)
        self.orthomcl_filter_fasta(compliant_fasta_dir)
        protdb = "goodProteins.fasta"  # created by orthomclFilterFasta
        blast_output = self.run_blast(protdb)
        sim_seq_file = self.orthomcl_blast_parser(compliant_fasta_dir, blast_output)
        self.load_blast_output_to_db(orthomcl_cfg, sim_seq_file)
        self.orthomcl_pairs(orthomcl_cfg)
        self.prepare_mcl_input(orthomcl_cfg)
        mcl_output_file = self.run_mcl()
        groups_file = self.orthomcl_group_mcl_output(mcl_output_file)
        orthologs = [];
        ids_in_orths = {};
        cluster_ind = self.parse_orthomcl_groups(groups_file, feature_info, orthologs, 
                                                 ids_in_orths)
        self.add_single_gene_families(feature_info, orthologs, ids_in_orths, cluster_ind)
        return self.save_pangenome_and_report(genome_refs, orthologs)


    def startup_mysql(self):
        self.log_line("Starting mysql service")
        self.log_process(subprocess.Popen(["service", "mysql", "start"], 
                cwd=self.scratch, stdout=subprocess.PIPE, stderr=subprocess.PIPE))

    def prepare_mysql_db(self):
        self.log_line("Preparing database")
        db = MySQLdb.connect(host="localhost", user="root", passwd="12345");
        cur = db.cursor()
        cur.execute("DROP DATABASE IF EXISTS orthomcl")
        cur.execute("CREATE DATABASE orthomcl")
        cur.close()
        db.close()

    def prepare_othomcl_config(self):
        self.log_line("Preparing orthomcl config file")
        orthomcl_cfg = self.scratch + "/orthomcl.cfg"
        f = open(orthomcl_cfg, "w")
        f.write("dbVendor=mysql\n");
        f.write("dbConnectString=dbi:mysql:orthomcl:mysql_local_infile=1:localhost:" + 
                "3306\n")
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
        return orthomcl_cfg

    def orthomcl_install_schema(self, orthomcl_cfg):
        self.log_line("Running orthomclInstallSchema")
        self.log_process(subprocess.Popen(["perl", self.plbin + "/orthomclInstallSchema", 
                orthomcl_cfg], cwd=self.scratch, stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE))

    def load_genomeset_object(self):
        genomeset = None
        if "input_genomeset_ref" in self.params and self.params["input_genomeset_ref"] is not None:
            self.log_line("Loading GenomeSet object from workspace")
            genomeset = self.ws.get_objects([{"ref": self.params["input_genomeset_ref"]}])[0]["data"]
        return genomeset

    def prepare_genome_refs(self, genomeset):
        self.log_line("Preparing genome refs")
        genome_refs = []
        if genomeset is not None:
            for param_key in genomeset["elements"]:
                genome_refs.append(genomeset["elements"][param_key]["ref"])
            self.log_line("Genome references from genome set: " + ", ".join(genome_refs))
        if "input_genome_refs" in self.params and self.params["input_genome_refs"] is not None:
            for genome_ref in self.params["input_genome_refs"]:
                if genome_ref is not None:
                    genome_refs.append(genome_ref)
            self.log_line("Final list of genome references: " + ", ".join(genome_refs))
        if len(genome_refs) < 2:
            raise ValueError("Number of genomes should be more than 1")
        if len(genome_refs) > 20:
            self.log_line("WARNING! Number of genomes exceeds 20, which can make " +
            "all-against-all blastp working unexpectedly long.")
        return genome_refs

    def load_genome_features_prepare_fasta(self, genome_refs, compliant_fasta_dir):
        feature_info = {}
        os.makedirs(compliant_fasta_dir)
        for genome_pos, genome_ref in enumerate(genome_refs):
            ############################# Genome loading ##########################
            self.log_line("Loading Genome object from workspace for ref [" + 
                          genome_ref + "]")
            info = self.ws.get_object_info_new({"objects": [{"ref": genome_ref}]})[0]
            genome_ref = str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])
            gaapi = GenomeAnnotationAPI(os.environ['SDK_CALLBACK_URL'], token=self.token)
            genome_combined = gaapi.get_combined_data({"ref": genome_ref, "exclude_genes": 1, 
                                                       "exclude_summary": 1})
            cds_map = genome_combined["feature_by_id_by_type"][genome_combined["cds_type"]]
            protein_map = genome_combined["protein_by_cds_id"]
            cds_ids = list(cds_map.keys())
            ############################# Features + Fasta ##########################
            self.log_line("Preparing fasta file for ref [" + genome_ref + "]")
            genome_id = str(genome_pos + 1)
            records = []
            for feature_pos, feature_id in enumerate(cds_ids):
                cds = cds_map[feature_id]
                if feature_id not in protein_map:
                    continue
                protein = protein_map[feature_id]
                if "protein_amino_acid_sequence" in protein:
                    sequence = protein["protein_amino_acid_sequence"]
                    id = str(feature_pos + 1)
                    record = SeqRecord(Seq(sequence), id=id, description="")
                    records.append(record)
                    func = None
                    if "protein_function" in protein:
                        func = protein["protein_function"]
                    if ((not func) or len(func) == 0) and "feature_function" in cds:
                        func = cds["feature_function"]
                    feature_info[genome_id + "|" + id] = {"fid": feature_id, "fpos": 
                            feature_pos, "gref": genome_ref, "func": func}
            fasta_file = self.scratch + "/" + genome_id + ".fasta"
            SeqIO.write(records, fasta_file, "fasta")
            ############################# Adjusting Fasta by Orthomcl ##########################
            self.log_line("Running orthomclAdjustFasta for ref [" + genome_ref + "]")
            self.log_process(subprocess.Popen(["perl", self.plbin + "/orthomclAdjustFasta", 
                    genome_id, fasta_file, "1"], cwd=compliant_fasta_dir, 
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE))
        return feature_info

    def orthomcl_filter_fasta(self, compliant_fasta_dir):
        self.log_line("Running orthomclFilterFasta")
        self.log_process(subprocess.Popen(["perl", self.plbin + "/orthomclFilterFasta", 
                compliant_fasta_dir, "50", "10"], cwd=self.scratch, 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE))

    def run_blast(self, protdb):
        ############################# Formatdb ##########################
        self.log_line("Running formatdb")
        self.log_process(subprocess.Popen(["formatdb", "-i", protdb],
                cwd=self.scratch, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
        ############################# BLAST ##########################
        self.log_line("Running blastp")
        blastp_args = ["blastall", "-p", "blastp", "-d", protdb, "-i", protdb, "-F", "m S", 
                       "-v", self.get_param(self.params, "num_descriptions", "100000"), 
                       "-b", self.get_param(self.params, "num_alignments", "100000"), 
                       "-e", self.get_param(self.params, "evalue", "1e-5"), 
                       "-m", "8",   # Alignment view is tabular (for orthomclBlastParser)
                       "-a", "1"]   # Number of processors is always 1
        self.add_param(self.params, "word_size", "-W", blastp_args)
        self.add_param(self.params, "gapopen", "-G", blastp_args)
        self.add_param(self.params, "gapextend", "-E", blastp_args)
        self.add_param(self.params, "matrix", "-M", blastp_args)
        self.add_param(self.params, "threshold", "-f", blastp_args)
        self.add_param(self.params, "comp_based_stats", "-C", blastp_args)
        self.add_param(self.params, "seg", "-F", blastp_args)
        self.add_param(self.params, "lcase_masking", "-U", blastp_args, True)
        self.add_param(self.params, "xdrop_gap_final", "-Z", blastp_args)
        self.add_param(self.params, "window_size", "-A", blastp_args)
        self.add_param(self.params, "use_sw_tback", "-s", blastp_args, True)
        self.log_line("Blastp command line: " + " ".join(blastp_args))
        blast_output = self.scratch + "/blastres.txt"
        with open(blast_output, "w") as outfile:
            self.log_process(subprocess.Popen(blastp_args, cwd=self.scratch, 
                    stdout=outfile, stderr=subprocess.PIPE))
        return blast_output

    def orthomcl_blast_parser(self, compliant_fasta_dir, blast_output):
        self.log_line("Running orthomclBlastParser")
        sim_seq_file = self.scratch + "/similarSequences.txt"
        with open(sim_seq_file, "w") as outfile:
            self.log_process(subprocess.Popen(["perl", self.plbin + "/orthomclBlastParser", 
                    blast_output, compliant_fasta_dir], cwd=self.scratch, stdout=outfile, 
                    stderr=subprocess.PIPE))
        return sim_seq_file

    def load_blast_output_to_db(self, orthomcl_cfg, sim_seq_file):
        self.log_line("Running orthomclLoadBlast")
        self.log_process(subprocess.Popen(["perl", self.plbin + "/orthomclLoadBlast", 
                orthomcl_cfg, sim_seq_file], cwd=self.scratch, stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE))

    def orthomcl_pairs(self, orthomcl_cfg):
        self.log_line("Running orthomclPairs")
        orthomcl_pairs_file = self.scratch + "/orthomcl_pairs.log"
        self.log_process(subprocess.Popen(["perl", self.plbin + "/orthomclPairs", 
                orthomcl_cfg, orthomcl_pairs_file, "cleanup=no"], cwd=self.scratch, 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE))
        return orthomcl_pairs_file

    def prepare_mcl_input(self, orthomcl_cfg):
        self.log_line("Running orthomclDumpPairsFiles")
        self.log_process(subprocess.Popen(["perl", self.plbin + "/orthomclDumpPairsFiles", 
                orthomcl_cfg], cwd=self.scratch, stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE))

    def run_mcl(self):
        self.log_line("Running mcl")
        mcl_output_file = self.scratch + "/mclOutput"
        mcl_args = ["mcl", "mclInput", "--abc", 
                    "-I", self.get_param(self.params, "mcl_main_i", "1.5"), 
                    "-o", mcl_output_file]
        self.add_param(self.params, "mcl_p", "-P", mcl_args)
        self.add_param(self.params, "mcl_s", "-S", mcl_args)
        self.add_param(self.params, "mcl_r", "-R", mcl_args)
        self.add_param(self.params, "mcl_pct", "-pct", mcl_args)
        self.add_param(self.params, "mcl_warn_p", "-warn-pct", mcl_args)
        self.add_param(self.params, "mcl_warn_factor", "-warn-factor", mcl_args)
        self.add_param(self.params, "mcl_init_l", "-l", mcl_args)
        self.add_param(self.params, "mcl_main_l", "-L", mcl_args)
        self.add_param(self.params, "mcl_init_i", "-i", mcl_args)
        self.log_line("Mcl command line: " + " ".join(mcl_args))
        self.log_process(subprocess.Popen(mcl_args, cwd=self.scratch, 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE))
        return mcl_output_file

    def orthomcl_group_mcl_output(self, mcl_output_file):
        self.log_line("Running orthomclMclToGroups")
        groups_file = self.scratch + "/groups.txt"
        with open(groups_file, "w") as outfile, open(mcl_output_file, "r") as infile:
            self.log_process(subprocess.Popen(["perl", self.plbin + "/orthomclMclToGroups", 
                    "grp", "1000"], cwd=self.scratch, stdin=infile, stdout=outfile, 
                    stderr=subprocess.PIPE))
        return groups_file

    def parse_orthomcl_groups(self, groups_file, feature_info, orthologs, ids_in_orths):
        self.log_line("Parsing groups file")
        cluster_ind = 0
        with open(groups_file, "r") as infile:
            for line_pos, line in enumerate(infile.readlines()):
                cluster_ind = line_pos + 1
                cluster_id = "cluster" + str(cluster_ind)
                function = ""
                items = []
                words = line.rstrip().split(" ")
                for id in words[1:]:
                    feature = feature_info[id]
                    items.append([feature["fid"], feature["fpos"], feature["gref"]])
                    func = feature["func"]
                    if func is not None and len(func) > len(function):
                        function = func
                    ids_in_orths[id] = True
                orthologs.append({"function": function, "id": cluster_id, 
                                  "orthologs": items})
        return cluster_ind

    def add_single_gene_families(self, feature_info, orthologs, ids_in_orths, cluster_ind):
        self.log_line("Adding single-gene families (they're not reported by OrthoMCL)")
        singles = 0
        for id in feature_info:
            if id in ids_in_orths:
                continue
            cluster_ind += 1
            singles += 1
            cluster_id = "cluster" + str(cluster_ind)
            feature = feature_info[id]
            function = feature["func"]
            items = [[feature["fid"], feature["fpos"], feature["gref"]]]
            orthologs.append({"function": function, "id": cluster_id, 
                              "orthologs": items})
        self.log_line(str(singles) + " single-gene families were added")

    def save_pangenome_and_report(self, genome_refs, orthologs):
        self.log_line("Saving pangenome object")
        output_obj_name = self.params["output_pangenome_id"]
        pangenome = {"genome_refs": genome_refs, "id": output_obj_name, "name": 
                     output_obj_name, "orthologs": orthologs, "type": "orthomcl"}
        input_ws_objects = []
        if "input_genomeset_ref" in self.params and self.params["input_genomeset_ref"] is not None:
            input_ws_objects.append(self.params["input_genomeset_ref"])
        if "input_genome_refs" in self.params and self.params["input_genome_refs"] is not None:
            for genome_ref in self.params["input_genome_refs"]:
                if genome_ref is not None:
                    input_ws_objects.append(genome_ref)
        self.provenance[0]["input_ws_objects"] = input_ws_objects
        self.provenance[0]["description"] = "Orthologous groups construction using OrthoMCL tool"
        info = self.ws.save_objects({"workspace": self.params["output_workspace"], "objects":
            [{"type": "KBaseGenomes.Pangenome", "name": output_obj_name, 
              "data": pangenome, "provenance": self.provenance}]})[0]
        pangenome_ref = str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])
        report = "Input genomes: " + str(len(genome_refs)) + "\n" + \
            "Output orthologs: " + str(len(orthologs)) + "\n"
        report_obj = {"objects_created": [{"ref": pangenome_ref, 
            "description": "Pangenome object"}], "text_message": report}
        report_name = "orthomcl_report_" + str(hex(uuid.getnode()))
        report_info = self.ws.save_objects({"workspace": self.params["output_workspace"],
            "objects": [{"type": "KBaseReport.Report", "data": report_obj, 
                "name": report_name, "meta": {}, "hidden": 1, "provenance": self.provenance}]})[0]
        return {"pangenome_ref": pangenome_ref, 
                "report_name": report_name, "report_ref": str(report_info[6]) + "/" + 
                str(report_info[0]) + "/" + str(report_info[4])}

    def log_line(self, line):
        self.log += line + "\n"
        print(line)

    def log_lines(self, lines):
        for line in lines:
            if len(line) > 0:
                self.log_line("|" + line)

    def log_process(self, process):
        process_out = process.communicate()
        output = process_out[0]
        if output is not None and len(output) > 0:
            self.log_line("Output:")
            self.log_lines(output.split("\n"))
        errors = process_out[1]
        if errors is not None and len(errors) > 0:
            self.log_line("Errors:")
            self.log_lines(errors.split("\n"))

    def get_param(self, params, param_name, def_value):
        ret = None
        if param_name in params and params[param_name] is not None and \
                len(str(params[param_name])) > 0:
            ret = str(params[param_name])
        else:
            ret = str(def_value)
        return ret

    def add_param(self, params, param_name, cli_arg, target_args, bool=False):
        if param_name in params and params[param_name] is not None and \
                len(str(params[param_name])) > 0:
            value = params[param_name]
            if bool:
                if value == 1:
                    target_args.append(cli_arg)
            else:
                target_args.append(cli_arg)
                target_args.append(str(value))
