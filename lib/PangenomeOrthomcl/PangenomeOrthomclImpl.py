# -*- coding: utf-8 -*-
#BEGIN_HEADER
import traceback
from PangenomeOrthomcl.PangenomeOrthomclBuilder import PangenomeOrthomclBuilder
#END_HEADER


class PangenomeOrthomcl:
    '''
    Module Name:
    PangenomeOrthomcl

    Module Description:
    A KBase module: PangenomeOrthomcl
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.6"
    GIT_URL = "https://github.com/kbaseapps/PangenomeOrthomcl"
    GIT_COMMIT_HASH = "3c1358e049e76c7581484cc8fa358ae8a57b0fd7"

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
        """
        :param params: instance of type "BuildPangenomeWithOrthmclParams"
           (Input parameters of build_pangenome_with_orthomcl method.
           input_genomeset_ref - optional input reference to genome set
           object (alternative way is input_genome_refs); input_genome_refs -
           optional input list of references to genome objects (alternative
           way is input_genomeset_ref); output_workspace - workspace for
           saving resulting pangenome; output_pangenome_id - name of
           resulting pangenome object; num_descriptions - [blastp, -v] Store
           one-line descriptions for this number of database sequences.
           Default value is 100000. num_alignments - [blastp, -b] Store
           alignments for this number of database sequences. Default value is
           100000. evalue - [blastp, -e] Expect value (E) for saving hits.
           Default value is 1e-5. word_size - [blastp, -W] Word size of
           initial match. Valid word sizes are 2-7. Default value is 3.
           gapopen - [blastp, -G] Cost to open a gap. Default value is 11.
           gapextend - [blastp, -E] Cost to extend a gap. Default value is 1.
           matrix - [blastp, -M] Scoring matrix name. Default value is
           BLOSUM62. threshold - [blastp, -f] Minimum score to add a word to
           the BLAST lookup table. Default value is 11. comp_based_stats -
           [blastp, -C] Use composition-based statistics (0: no
           composition-based statistics; 1: Composition-based statistics as
           in NAR 29:2994-3005, 2001; 2: Composition-based score adjustments
           as in Bioinformatics 21:902-911, 2005, conditioned on sequence
           properties; 3 - Composition-based score adjustment as in
           Bioinformatics 21:902-911, 2005, unconditionally). Default value
           is 2. seg - [blastp, -F] Filter query sequence with SEG (yes/no).
           Default value is yes. lcase_masking - [blastp, -U] Use lower case
           filtering in query and subject sequence(s). Default value is
           false(0). xdrop_gap_final - [blastp, -Z] Heuristic value (in bits)
           for final gapped alignment. Default value is 25. window_size -
           [blastp, -A] Multiple hits window size, use 0 to specify 1-hit
           algorithm. Default value is 40. use_sw_tback - [blastp, -s]
           Compute locally optimal Smith-Waterman alignments. Default value
           is false(0). mcl_p - [mcl, -P] Prune number. Default value is
           10000. mcl_s - [mcl, -S] Selection number. Default value is 1100.
           mcl_r - [mcl, -R] Recovery number. Default value is 1400. mcl_pct
           - [mcl, -pct] Recovery percentage. Default value is 90. mcl_warn_p
           - [mcl, -warn-p] Warn if pruning reduces mass to this weight.
           Default value is 10. mcl_warn_factor - [mcl, -warn-factor] Warn if
           pruning reduces entry count by this value. Default value is 1000.
           mcl_init_l - [mcl, -l] Initial loop length. Default value is 0.
           mcl_main_l - [mcl, -L] Main loop length. Default value is 10000.
           mcl_init_i - [mcl, -i] Initial inflation. Default value is 2.0.
           mcl_main_i - [mcl, -I] Main inflation. Default value is 1.5.) ->
           structure: parameter "input_genomeset_ref" of type
           "ws_genomeset_id" (The workspace ID for a GenomeSet data object.
           @id ws KBaseSearch.GenomeSet), parameter "input_genome_refs" of
           list of type "ws_genome_id" (The workspace ID for a GenomeSet data
           object. @id ws KBaseGenome.Genome), parameter "output_workspace"
           of type "workspace" (Name of workspace.), parameter
           "output_pangenome_id" of type "ws_pangenome_id" (The workspace ID
           for a Pangenome data object. @id ws KBaseGenomes.Pangenome),
           parameter "num_descriptions" of Long, parameter "num_alignments"
           of Long, parameter "evalue" of String, parameter "word_size" of
           Long, parameter "gapopen" of Long, parameter "gapextend" of Long,
           parameter "matrix" of String, parameter "threshold" of Long,
           parameter "comp_based_stats" of String, parameter "seg" of String,
           parameter "lcase_masking" of type "boolean" (Indicates true or
           false values, false = 0, true = 1 @range [0,1]), parameter
           "xdrop_gap_final" of Double, parameter "window_size" of Long,
           parameter "use_sw_tback" of type "boolean" (Indicates true or
           false values, false = 0, true = 1 @range [0,1]), parameter "mcl_p"
           of Long, parameter "mcl_s" of Long, parameter "mcl_r" of Long,
           parameter "mcl_pct" of Long, parameter "mcl_warn_p" of Long,
           parameter "mcl_warn_factor" of Long, parameter "mcl_init_l" of
           Long, parameter "mcl_main_l" of Long, parameter "mcl_init_i" of
           Double, parameter "mcl_main_i" of Double
        :returns: instance of type "BuildPangenomeWithOrthmclResult" (Output
           results of build_pangenome_with_orthomcl method. One of
           'pangenome_ref' and 'error' fields should be defined.) ->
           structure: parameter "pangenome_ref" of type "ws_pangenome_id"
           (The workspace ID for a Pangenome data object. @id ws
           KBaseGenomes.Pangenome), parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN build_pangenome_with_orthomcl
        try:
            runner = PangenomeOrthomclBuilder(self.scratch, self.workspaceURL,
                                              params, ctx["token"], ctx["provenance"])
            returnVal = runner.run()
        except Exception, err:
            raise ValueError(traceback.format_exc())
        #END build_pangenome_with_orthomcl

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method build_pangenome_with_orthomcl return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
