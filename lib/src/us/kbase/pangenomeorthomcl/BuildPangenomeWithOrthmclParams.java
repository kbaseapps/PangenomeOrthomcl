
package us.kbase.pangenomeorthomcl;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: BuildPangenomeWithOrthmclParams</p>
 * <pre>
 * Input parameters of build_pangenome_with_orthomcl method.
 * input_genomeset_ref - optional input reference to genome set 
 *     object (alternative way is input_genome_refs);
 * input_genome_refs - optional input list of references to
 *     genome objects (alternative way is input_genomeset_ref);
 * output_workspace - workspace for saving resulting pangenome;
 * output_pangenome_id - name of resulting pangenome object;
 * num_descriptions - [blastp, -v] Store one-line descriptions for 
 *     this number of database sequences. Default value is 100000.
 * num_alignments - [blastp, -b] Store alignments for this number of 
 *     database sequences. Default value is 100000.
 * evalue - [blastp, -e] Expect value (E) for saving hits. Default
 *     value is 1e-5.
 * word_size - [blastp, -W] Word size of initial match. Valid word 
 *     sizes are 2-7. Default value is 3.
 * gapopen - [blastp, -G] Cost to open a gap. Default value is 11.
 * gapextend - [blastp, -E] Cost to extend a gap. Default value is 1.
 * matrix - [blastp, -M] Scoring matrix name. Default value is BLOSUM62.
 * threshold - [blastp, -f] Minimum score to add a word to the BLAST 
 *     lookup table. Default value is 11.
 * comp_based_stats - [blastp, -C] Use composition-based statistics 
 *     (0: no composition-based statistics; 1: Composition-based 
 *     statistics as in NAR 29:2994-3005, 2001; 2: Composition-based 
 *     score adjustments as in Bioinformatics 21:902-911, 2005, 
 *     conditioned on sequence properties; 3 - Composition-based 
 *     score adjustment as in Bioinformatics 21:902-911, 2005, 
 *     unconditionally). Default value is 2.
 * seg - [blastp, -F] Filter query sequence with SEG (yes/no). Default
 *     value is yes.
 * lcase_masking - [blastp, -U] Use lower case filtering in query and 
 *     subject sequence(s). Default value is false(0).
 * xdrop_gap_final - [blastp, -Z] Heuristic value (in bits) for final 
 *     gapped alignment. Default value is 25.
 * window_size - [blastp, -A] Multiple hits window size, use 0 to 
 *     specify 1-hit algorithm. Default value is 40.
 * use_sw_tback - [blastp, -s] Compute locally optimal Smith-Waterman 
 *     alignments. Default value is false(0).
 * mcl_p - [mcl, -P] Prune number. Default value is 10000.
 * mcl_s - [mcl, -S] Selection number. Default value is 1100.
 * mcl_r - [mcl, -R] Recovery number. Default value is 1400.
 * mcl_pct - [mcl, -pct] Recovery percentage. Default value is 90.
 * mcl_warn_p - [mcl, -warn-p] Warn if pruning reduces mass to this 
 *     weight. Default value is 10.
 * mcl_warn_factor - [mcl, -warn-factor] Warn if pruning reduces entry 
 *     count by this value. Default value is 1000.
 * mcl_init_l - [mcl, -l] Initial loop length. Default value is 0.
 * mcl_main_l - [mcl, -L] Main loop length. Default value is 10000.
 * mcl_init_i - [mcl, -i] Initial inflation. Default value is 2.0.
 * mcl_main_i - [mcl, -I] Main inflation. Default value is 1.5.
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "input_genomeset_ref",
    "input_genome_refs",
    "output_workspace",
    "output_pangenome_id",
    "num_descriptions",
    "num_alignments",
    "evalue",
    "word_size",
    "gapopen",
    "gapextend",
    "matrix",
    "threshold",
    "comp_based_stats",
    "seg",
    "lcase_masking",
    "xdrop_gap_final",
    "window_size",
    "use_sw_tback",
    "mcl_p",
    "mcl_s",
    "mcl_r",
    "mcl_pct",
    "mcl_warn_p",
    "mcl_warn_factor",
    "mcl_init_l",
    "mcl_main_l",
    "mcl_init_i",
    "mcl_main_i"
})
public class BuildPangenomeWithOrthmclParams {

    @JsonProperty("input_genomeset_ref")
    private java.lang.String inputGenomesetRef;
    @JsonProperty("input_genome_refs")
    private List<String> inputGenomeRefs;
    @JsonProperty("output_workspace")
    private java.lang.String outputWorkspace;
    @JsonProperty("output_pangenome_id")
    private java.lang.String outputPangenomeId;
    @JsonProperty("num_descriptions")
    private Long numDescriptions;
    @JsonProperty("num_alignments")
    private Long numAlignments;
    @JsonProperty("evalue")
    private java.lang.String evalue;
    @JsonProperty("word_size")
    private Long wordSize;
    @JsonProperty("gapopen")
    private Long gapopen;
    @JsonProperty("gapextend")
    private Long gapextend;
    @JsonProperty("matrix")
    private java.lang.String matrix;
    @JsonProperty("threshold")
    private Long threshold;
    @JsonProperty("comp_based_stats")
    private java.lang.String compBasedStats;
    @JsonProperty("seg")
    private java.lang.String seg;
    @JsonProperty("lcase_masking")
    private Long lcaseMasking;
    @JsonProperty("xdrop_gap_final")
    private Double xdropGapFinal;
    @JsonProperty("window_size")
    private Long windowSize;
    @JsonProperty("use_sw_tback")
    private Long useSwTback;
    @JsonProperty("mcl_p")
    private Long mclP;
    @JsonProperty("mcl_s")
    private Long mclS;
    @JsonProperty("mcl_r")
    private Long mclR;
    @JsonProperty("mcl_pct")
    private Long mclPct;
    @JsonProperty("mcl_warn_p")
    private Long mclWarnP;
    @JsonProperty("mcl_warn_factor")
    private Long mclWarnFactor;
    @JsonProperty("mcl_init_l")
    private Long mclInitL;
    @JsonProperty("mcl_main_l")
    private Long mclMainL;
    @JsonProperty("mcl_init_i")
    private Double mclInitI;
    @JsonProperty("mcl_main_i")
    private Double mclMainI;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("input_genomeset_ref")
    public java.lang.String getInputGenomesetRef() {
        return inputGenomesetRef;
    }

    @JsonProperty("input_genomeset_ref")
    public void setInputGenomesetRef(java.lang.String inputGenomesetRef) {
        this.inputGenomesetRef = inputGenomesetRef;
    }

    public BuildPangenomeWithOrthmclParams withInputGenomesetRef(java.lang.String inputGenomesetRef) {
        this.inputGenomesetRef = inputGenomesetRef;
        return this;
    }

    @JsonProperty("input_genome_refs")
    public List<String> getInputGenomeRefs() {
        return inputGenomeRefs;
    }

    @JsonProperty("input_genome_refs")
    public void setInputGenomeRefs(List<String> inputGenomeRefs) {
        this.inputGenomeRefs = inputGenomeRefs;
    }

    public BuildPangenomeWithOrthmclParams withInputGenomeRefs(List<String> inputGenomeRefs) {
        this.inputGenomeRefs = inputGenomeRefs;
        return this;
    }

    @JsonProperty("output_workspace")
    public java.lang.String getOutputWorkspace() {
        return outputWorkspace;
    }

    @JsonProperty("output_workspace")
    public void setOutputWorkspace(java.lang.String outputWorkspace) {
        this.outputWorkspace = outputWorkspace;
    }

    public BuildPangenomeWithOrthmclParams withOutputWorkspace(java.lang.String outputWorkspace) {
        this.outputWorkspace = outputWorkspace;
        return this;
    }

    @JsonProperty("output_pangenome_id")
    public java.lang.String getOutputPangenomeId() {
        return outputPangenomeId;
    }

    @JsonProperty("output_pangenome_id")
    public void setOutputPangenomeId(java.lang.String outputPangenomeId) {
        this.outputPangenomeId = outputPangenomeId;
    }

    public BuildPangenomeWithOrthmclParams withOutputPangenomeId(java.lang.String outputPangenomeId) {
        this.outputPangenomeId = outputPangenomeId;
        return this;
    }

    @JsonProperty("num_descriptions")
    public Long getNumDescriptions() {
        return numDescriptions;
    }

    @JsonProperty("num_descriptions")
    public void setNumDescriptions(Long numDescriptions) {
        this.numDescriptions = numDescriptions;
    }

    public BuildPangenomeWithOrthmclParams withNumDescriptions(Long numDescriptions) {
        this.numDescriptions = numDescriptions;
        return this;
    }

    @JsonProperty("num_alignments")
    public Long getNumAlignments() {
        return numAlignments;
    }

    @JsonProperty("num_alignments")
    public void setNumAlignments(Long numAlignments) {
        this.numAlignments = numAlignments;
    }

    public BuildPangenomeWithOrthmclParams withNumAlignments(Long numAlignments) {
        this.numAlignments = numAlignments;
        return this;
    }

    @JsonProperty("evalue")
    public java.lang.String getEvalue() {
        return evalue;
    }

    @JsonProperty("evalue")
    public void setEvalue(java.lang.String evalue) {
        this.evalue = evalue;
    }

    public BuildPangenomeWithOrthmclParams withEvalue(java.lang.String evalue) {
        this.evalue = evalue;
        return this;
    }

    @JsonProperty("word_size")
    public Long getWordSize() {
        return wordSize;
    }

    @JsonProperty("word_size")
    public void setWordSize(Long wordSize) {
        this.wordSize = wordSize;
    }

    public BuildPangenomeWithOrthmclParams withWordSize(Long wordSize) {
        this.wordSize = wordSize;
        return this;
    }

    @JsonProperty("gapopen")
    public Long getGapopen() {
        return gapopen;
    }

    @JsonProperty("gapopen")
    public void setGapopen(Long gapopen) {
        this.gapopen = gapopen;
    }

    public BuildPangenomeWithOrthmclParams withGapopen(Long gapopen) {
        this.gapopen = gapopen;
        return this;
    }

    @JsonProperty("gapextend")
    public Long getGapextend() {
        return gapextend;
    }

    @JsonProperty("gapextend")
    public void setGapextend(Long gapextend) {
        this.gapextend = gapextend;
    }

    public BuildPangenomeWithOrthmclParams withGapextend(Long gapextend) {
        this.gapextend = gapextend;
        return this;
    }

    @JsonProperty("matrix")
    public java.lang.String getMatrix() {
        return matrix;
    }

    @JsonProperty("matrix")
    public void setMatrix(java.lang.String matrix) {
        this.matrix = matrix;
    }

    public BuildPangenomeWithOrthmclParams withMatrix(java.lang.String matrix) {
        this.matrix = matrix;
        return this;
    }

    @JsonProperty("threshold")
    public Long getThreshold() {
        return threshold;
    }

    @JsonProperty("threshold")
    public void setThreshold(Long threshold) {
        this.threshold = threshold;
    }

    public BuildPangenomeWithOrthmclParams withThreshold(Long threshold) {
        this.threshold = threshold;
        return this;
    }

    @JsonProperty("comp_based_stats")
    public java.lang.String getCompBasedStats() {
        return compBasedStats;
    }

    @JsonProperty("comp_based_stats")
    public void setCompBasedStats(java.lang.String compBasedStats) {
        this.compBasedStats = compBasedStats;
    }

    public BuildPangenomeWithOrthmclParams withCompBasedStats(java.lang.String compBasedStats) {
        this.compBasedStats = compBasedStats;
        return this;
    }

    @JsonProperty("seg")
    public java.lang.String getSeg() {
        return seg;
    }

    @JsonProperty("seg")
    public void setSeg(java.lang.String seg) {
        this.seg = seg;
    }

    public BuildPangenomeWithOrthmclParams withSeg(java.lang.String seg) {
        this.seg = seg;
        return this;
    }

    @JsonProperty("lcase_masking")
    public Long getLcaseMasking() {
        return lcaseMasking;
    }

    @JsonProperty("lcase_masking")
    public void setLcaseMasking(Long lcaseMasking) {
        this.lcaseMasking = lcaseMasking;
    }

    public BuildPangenomeWithOrthmclParams withLcaseMasking(Long lcaseMasking) {
        this.lcaseMasking = lcaseMasking;
        return this;
    }

    @JsonProperty("xdrop_gap_final")
    public Double getXdropGapFinal() {
        return xdropGapFinal;
    }

    @JsonProperty("xdrop_gap_final")
    public void setXdropGapFinal(Double xdropGapFinal) {
        this.xdropGapFinal = xdropGapFinal;
    }

    public BuildPangenomeWithOrthmclParams withXdropGapFinal(Double xdropGapFinal) {
        this.xdropGapFinal = xdropGapFinal;
        return this;
    }

    @JsonProperty("window_size")
    public Long getWindowSize() {
        return windowSize;
    }

    @JsonProperty("window_size")
    public void setWindowSize(Long windowSize) {
        this.windowSize = windowSize;
    }

    public BuildPangenomeWithOrthmclParams withWindowSize(Long windowSize) {
        this.windowSize = windowSize;
        return this;
    }

    @JsonProperty("use_sw_tback")
    public Long getUseSwTback() {
        return useSwTback;
    }

    @JsonProperty("use_sw_tback")
    public void setUseSwTback(Long useSwTback) {
        this.useSwTback = useSwTback;
    }

    public BuildPangenomeWithOrthmclParams withUseSwTback(Long useSwTback) {
        this.useSwTback = useSwTback;
        return this;
    }

    @JsonProperty("mcl_p")
    public Long getMclP() {
        return mclP;
    }

    @JsonProperty("mcl_p")
    public void setMclP(Long mclP) {
        this.mclP = mclP;
    }

    public BuildPangenomeWithOrthmclParams withMclP(Long mclP) {
        this.mclP = mclP;
        return this;
    }

    @JsonProperty("mcl_s")
    public Long getMclS() {
        return mclS;
    }

    @JsonProperty("mcl_s")
    public void setMclS(Long mclS) {
        this.mclS = mclS;
    }

    public BuildPangenomeWithOrthmclParams withMclS(Long mclS) {
        this.mclS = mclS;
        return this;
    }

    @JsonProperty("mcl_r")
    public Long getMclR() {
        return mclR;
    }

    @JsonProperty("mcl_r")
    public void setMclR(Long mclR) {
        this.mclR = mclR;
    }

    public BuildPangenomeWithOrthmclParams withMclR(Long mclR) {
        this.mclR = mclR;
        return this;
    }

    @JsonProperty("mcl_pct")
    public Long getMclPct() {
        return mclPct;
    }

    @JsonProperty("mcl_pct")
    public void setMclPct(Long mclPct) {
        this.mclPct = mclPct;
    }

    public BuildPangenomeWithOrthmclParams withMclPct(Long mclPct) {
        this.mclPct = mclPct;
        return this;
    }

    @JsonProperty("mcl_warn_p")
    public Long getMclWarnP() {
        return mclWarnP;
    }

    @JsonProperty("mcl_warn_p")
    public void setMclWarnP(Long mclWarnP) {
        this.mclWarnP = mclWarnP;
    }

    public BuildPangenomeWithOrthmclParams withMclWarnP(Long mclWarnP) {
        this.mclWarnP = mclWarnP;
        return this;
    }

    @JsonProperty("mcl_warn_factor")
    public Long getMclWarnFactor() {
        return mclWarnFactor;
    }

    @JsonProperty("mcl_warn_factor")
    public void setMclWarnFactor(Long mclWarnFactor) {
        this.mclWarnFactor = mclWarnFactor;
    }

    public BuildPangenomeWithOrthmclParams withMclWarnFactor(Long mclWarnFactor) {
        this.mclWarnFactor = mclWarnFactor;
        return this;
    }

    @JsonProperty("mcl_init_l")
    public Long getMclInitL() {
        return mclInitL;
    }

    @JsonProperty("mcl_init_l")
    public void setMclInitL(Long mclInitL) {
        this.mclInitL = mclInitL;
    }

    public BuildPangenomeWithOrthmclParams withMclInitL(Long mclInitL) {
        this.mclInitL = mclInitL;
        return this;
    }

    @JsonProperty("mcl_main_l")
    public Long getMclMainL() {
        return mclMainL;
    }

    @JsonProperty("mcl_main_l")
    public void setMclMainL(Long mclMainL) {
        this.mclMainL = mclMainL;
    }

    public BuildPangenomeWithOrthmclParams withMclMainL(Long mclMainL) {
        this.mclMainL = mclMainL;
        return this;
    }

    @JsonProperty("mcl_init_i")
    public Double getMclInitI() {
        return mclInitI;
    }

    @JsonProperty("mcl_init_i")
    public void setMclInitI(Double mclInitI) {
        this.mclInitI = mclInitI;
    }

    public BuildPangenomeWithOrthmclParams withMclInitI(Double mclInitI) {
        this.mclInitI = mclInitI;
        return this;
    }

    @JsonProperty("mcl_main_i")
    public Double getMclMainI() {
        return mclMainI;
    }

    @JsonProperty("mcl_main_i")
    public void setMclMainI(Double mclMainI) {
        this.mclMainI = mclMainI;
    }

    public BuildPangenomeWithOrthmclParams withMclMainI(Double mclMainI) {
        this.mclMainI = mclMainI;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("BuildPangenomeWithOrthmclParams"+" [inputGenomesetRef=")+ inputGenomesetRef)+", inputGenomeRefs=")+ inputGenomeRefs)+", outputWorkspace=")+ outputWorkspace)+", outputPangenomeId=")+ outputPangenomeId)+", numDescriptions=")+ numDescriptions)+", numAlignments=")+ numAlignments)+", evalue=")+ evalue)+", wordSize=")+ wordSize)+", gapopen=")+ gapopen)+", gapextend=")+ gapextend)+", matrix=")+ matrix)+", threshold=")+ threshold)+", compBasedStats=")+ compBasedStats)+", seg=")+ seg)+", lcaseMasking=")+ lcaseMasking)+", xdropGapFinal=")+ xdropGapFinal)+", windowSize=")+ windowSize)+", useSwTback=")+ useSwTback)+", mclP=")+ mclP)+", mclS=")+ mclS)+", mclR=")+ mclR)+", mclPct=")+ mclPct)+", mclWarnP=")+ mclWarnP)+", mclWarnFactor=")+ mclWarnFactor)+", mclInitL=")+ mclInitL)+", mclMainL=")+ mclMainL)+", mclInitI=")+ mclInitI)+", mclMainI=")+ mclMainI)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
