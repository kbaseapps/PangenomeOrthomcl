
package us.kbase.pangenomeorthomcl;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: BuildPangenomeWithOrthmclResult</p>
 * <pre>
 * Output results of build_pangenome_with_orthomcl method.
 * One of 'pangenome_ref' and 'error' fields should be defined.
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "output_log",
    "pangenome_ref"
})
public class BuildPangenomeWithOrthmclResult {

    @JsonProperty("output_log")
    private String outputLog;
    @JsonProperty("pangenome_ref")
    private String pangenomeRef;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("output_log")
    public String getOutputLog() {
        return outputLog;
    }

    @JsonProperty("output_log")
    public void setOutputLog(String outputLog) {
        this.outputLog = outputLog;
    }

    public BuildPangenomeWithOrthmclResult withOutputLog(String outputLog) {
        this.outputLog = outputLog;
        return this;
    }

    @JsonProperty("pangenome_ref")
    public String getPangenomeRef() {
        return pangenomeRef;
    }

    @JsonProperty("pangenome_ref")
    public void setPangenomeRef(String pangenomeRef) {
        this.pangenomeRef = pangenomeRef;
    }

    public BuildPangenomeWithOrthmclResult withPangenomeRef(String pangenomeRef) {
        this.pangenomeRef = pangenomeRef;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((("BuildPangenomeWithOrthmclResult"+" [outputLog=")+ outputLog)+", pangenomeRef=")+ pangenomeRef)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
