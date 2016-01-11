/*
A KBase module: PangenomeOrthomcl
*/
module PangenomeOrthomcl {

    /*
        Indicates true or false values, false = 0, true = 1
        @range [0,1]
    */
    typedef int boolean;

    /*
        Name of workspace.
    */
    typedef string workspace;

    /* 
        The workspace ID for a GenomeSet data object.
        @id ws KBaseGenome.Genome
    */
    typedef string ws_genome_id;

    /* 
        The workspace ID for a GenomeSet data object.
        @id ws KBaseSearch.GenomeSet
    */
    typedef string ws_genomeset_id;

    /* 
        The workspace ID for a Pangenome data object.
        @id ws KBaseGenomes.Pangenome
    */
    typedef string ws_pangenome_id;

    /*
        Input parameters of build_pangenome_with_orthomcl method.
        input_genomeset_ref - optional input reference to genome set 
            object (alternative way is input_genome_refs);
        input_genome_refs - optional input list of references to
            genome objects (alternative way is input_genomeset_ref);
        output_workspace - workspace for saving resulting pangenome;
        output_pangenome_id - name of resulting pangenome object;
        num_descriptions - [blastp, -v] Store one-line descriptions for 
            this number of database sequences. Default value is 100000.
        num_alignments - [blastp, -b] Store alignments for this number of 
            database sequences. Default value is 100000.
        evalue - [blastp, -e] Expect value (E) for saving hits. Default
            value is 1e-5.
        word_size - [blastp, -W] Word size of initial match. Valid word 
            sizes are 2-7. Default value is 3.
        gapopen - [blastp, -G] Cost to open a gap. Default value is 11.
        gapextend - [blastp, -E] Cost to extend a gap. Default value is 1.
        matrix - [blastp, -M] Scoring matrix name. Default value is BLOSUM62.
        threshold - [blastp, -f] Minimum score to add a word to the BLAST 
            lookup table. Default value is 11.
        comp_based_stats - [blastp, -C] Use composition-based statistics 
            (0: no composition-based statistics; 1: Composition-based 
            statistics as in NAR 29:2994-3005, 2001; 2: Composition-based 
            score adjustments as in Bioinformatics 21:902-911, 2005, 
            conditioned on sequence properties; 3 - Composition-based 
            score adjustment as in Bioinformatics 21:902-911, 2005, 
            unconditionally). Default value is 2.
        seg - [blastp, -F] Filter query sequence with SEG (yes/no). Default
            value is yes.
        lcase_masking - [blastp, -U] Use lower case filtering in query and 
            subject sequence(s). Default value is false(0).
        xdrop_gap_final - [blastp, -Z] Heuristic value (in bits) for final 
            gapped alignment. Default value is 25.
        window_size - [blastp, -A] Multiple hits window size, use 0 to 
            specify 1-hit algorithm. Default value is 40.
        use_sw_tback - [blastp, -s] Compute locally optimal Smith-Waterman 
            alignments. Default value is false(0).
        mcl_p - [mcl, -P] Prune number. Default value is 10000.
        mcl_s - [mcl, -S] Selection number. Default value is 1100.
        mcl_r - [mcl, -R] Recovery number. Default value is 1400.
        mcl_pct - [mcl, -pct] Recovery percentage. Default value is 90.
        mcl_warn_p - [mcl, -warn-p] Warn if pruning reduces mass to this 
            weight. Default value is 10.
        mcl_warn_factor - [mcl, -warn-factor] Warn if pruning reduces entry 
            count by this value. Default value is 1000.
        mcl_init_l - [mcl, -l] Initial loop length. Default value is 0.
        mcl_main_l - [mcl, -L] Main loop length. Default value is 10000.
        mcl_init_i - [mcl, -i] Initial inflation. Default value is 2.0.
        mcl_main_i - [mcl, -I] Main inflation. Default value is 1.5.
    */
    typedef structure {
        ws_genomeset_id input_genomeset_ref;
        list<ws_genome_id> input_genome_refs;
        workspace output_workspace;
        ws_pangenome_id output_pangenome_id;
        int num_descriptions;
        int num_alignments;
        string evalue;
        int word_size;
        int gapopen;
        int gapextend;
        string matrix;
        int threshold;
        string comp_based_stats;
        string seg;
        boolean lcase_masking;
        float xdrop_gap_final;
        int window_size;
        boolean use_sw_tback;
        int mcl_p;
        int mcl_s;
        int mcl_r;
        int mcl_pct;
        int mcl_warn_p;
        int mcl_warn_factor;
        int mcl_init_l;
        int mcl_main_l;
        float mcl_init_i;
        float mcl_main_i;
    } BuildPangenomeWithOrthmclParams;

    /*
        Output results of build_pangenome_with_orthomcl method.
        One of 'pangenome_ref' and 'error' fields should be defined.
    */
    typedef structure {
        ws_pangenome_id pangenome_ref;
        string report_name;
        string report_ref;
    } BuildPangenomeWithOrthmclResult;

    funcdef build_pangenome_with_orthomcl(BuildPangenomeWithOrthmclParams params)
        returns (BuildPangenomeWithOrthmclResult) authentication required;
};