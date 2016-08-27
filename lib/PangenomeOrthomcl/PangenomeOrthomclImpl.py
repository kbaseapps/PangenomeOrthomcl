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
