import unittest
import os
import json
import time

from os import environ
from ConfigParser import ConfigParser
from pprint import pprint

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

from biokbase.workspace.client import Workspace as workspaceService
from PangenomeOrthomcl.PangenomeOrthomclImpl import PangenomeOrthomcl


class PangenomeOrthomclTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        cls.ctx = {'token': token}
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('PangenomeOrthomcl'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = PangenomeOrthomcl(cls.cfg)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_PangenomeOrthomcl_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_build_pangenome_with_orthomcl(self):
        contig_obj_name = "contigset.1"
        contig = {'id': '1', 'length': 10, 'md5': 'md5', 'sequence': 'agcttttcat'}
        obj = {'contigs': [contig], 'id': 'id', 'md5': 'md5', 'name': 'name', 
                'source': 'source', 'source_id': 'source_id', 'type': 'type'}
        self.getWsClient().save_objects({'workspace': self.getWsName(), 'objects':
            [{'type': 'KBaseGenomes.ContigSet', 'name': contig_obj_name, 'data': obj}]})
        genome_fasta_files = ["Escherichia_coli_042_uid161985.faa", "Escherichia_coli_BW2952_uid59391.faa", 
                              "Escherichia_coli_K12_MG1655_uid57779.faa"]
        genomeset_obj = {"description" : "", "elements" : {}}
        for genome_index, genome_file_name in enumerate(genome_fasta_files):
            test_dir = os.path.dirname(os.path.realpath(__file__))
            file_path = test_dir + "/data/" + genome_file_name
            features = []
            for record in SeqIO.parse(file_path, "fasta"):
                id = record.id
                sequence = str(record.seq)
                if len(sequence) <= 100:
                    features.append({"id" : id, "location" : [["1", 0, "+", 0]], "type" : "CDS", 
                                     "protein_translation" : sequence, "aliases" : [], 
                                     "annotations" : [], "function" : ""})
            genome_obj = {"complete" : 0, "contig_ids" : ["1"], "contig_lengths" : [10],
                          "contigset_ref" : self.getWsName() + "/" + contig_obj_name, 
                          "dna_size" : 10, "domain" : "Bacteria", "gc_content" : 0.5,
                          "genetic_code" : 11, "id" : genome_file_name, "md5" : "md5",
                          "num_contigs" : 1, "scientific_name" : genome_file_name,
                          "source" : "test folder", "source_id" : "noid", "features" : features}
            genome_obj_name = "genome." + str(genome_index)
            self.getWsClient().save_objects({'workspace': self.getWsName(), 'objects':
                    [{'type': 'KBaseGenomes.Genome', 'name': genome_obj_name, 'data': genome_obj}]})
            genomeset_obj["elements"]["param" + str(genome_index)] = {"ref" :
                    self.getWsName() + "/" + genome_obj_name}
        genomeset_obj_name = "genomeset.1"
        self.getWsClient().save_objects({'workspace': self.getWsName(), 'objects':
                [{'type': 'KBaseSearch.GenomeSet', 'name': genomeset_obj_name, 'data': genomeset_obj}]})
        output_pangenome_name = "pangenome.1"
        ret = self.getImpl().build_pangenome_with_orthomcl(self.getContext(), {
                "intput_genomeset_ref" : self.getWsName() + "/" + genomeset_obj_name,
                "output_workspace" : self.getWsName(), "output_pangenome_id" : output_pangenome_name})[0]
        print(ret["output_log"])
        pangenome = self.getWsClient().get_objects([{'ref': ret["pangenome_ref"]}])[0]['data']
        self.assertEqual(len(pangenome["orthologs"]), 350)
        pass
        