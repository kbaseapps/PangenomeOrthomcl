/*
A KBase module: PangenomeOrthomcl
*/
module PangenomeOrthomcl {

    typedef structure {
        string intput_genomeset_ref;
        string output_workspace;
        string output_pangenome_id;
    } BuildPangenomeWithOrthmclParams;

    /*
        One of 'pangenome_ref' and 'error' fields should be defined.
    */
    typedef structure {
        string output_log;
        string pangenome_ref;
    } BuildPangenomeWithOrthmclResult;

    funcdef build_pangenome_with_orthomcl(BuildPangenomeWithOrthmclParams params)
        returns (BuildPangenomeWithOrthmclResult) authentication required;
};