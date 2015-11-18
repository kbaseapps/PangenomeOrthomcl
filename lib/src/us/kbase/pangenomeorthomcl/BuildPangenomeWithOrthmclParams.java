
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
 * <p>Original spec-file type: BuildPangenomeWithOrthmclParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "intput_genomeset_ref",
    "output_workspace",
    "output_pangenome_id"
})
public class BuildPangenomeWithOrthmclParams {

    @JsonProperty("intput_genomeset_ref")
    private String intputGenomesetRef;
    @JsonProperty("output_workspace")
    private String outputWorkspace;
    @JsonProperty("output_pangenome_id")
    private String outputPangenomeId;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("intput_genomeset_ref")
    public String getIntputGenomesetRef() {
        return intputGenomesetRef;
    }

    @JsonProperty("intput_genomeset_ref")
    public void setIntputGenomesetRef(String intputGenomesetRef) {
        this.intputGenomesetRef = intputGenomesetRef;
    }

    public BuildPangenomeWithOrthmclParams withIntputGenomesetRef(String intputGenomesetRef) {
        this.intputGenomesetRef = intputGenomesetRef;
        return this;
    }

    @JsonProperty("output_workspace")
    public String getOutputWorkspace() {
        return outputWorkspace;
    }

    @JsonProperty("output_workspace")
    public void setOutputWorkspace(String outputWorkspace) {
        this.outputWorkspace = outputWorkspace;
    }

    public BuildPangenomeWithOrthmclParams withOutputWorkspace(String outputWorkspace) {
        this.outputWorkspace = outputWorkspace;
        return this;
    }

    @JsonProperty("output_pangenome_id")
    public String getOutputPangenomeId() {
        return outputPangenomeId;
    }

    @JsonProperty("output_pangenome_id")
    public void setOutputPangenomeId(String outputPangenomeId) {
        this.outputPangenomeId = outputPangenomeId;
    }

    public BuildPangenomeWithOrthmclParams withOutputPangenomeId(String outputPangenomeId) {
        this.outputPangenomeId = outputPangenomeId;
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
        return ((((((((("BuildPangenomeWithOrthmclParams"+" [intputGenomesetRef=")+ intputGenomesetRef)+", outputWorkspace=")+ outputWorkspace)+", outputPangenomeId=")+ outputPangenomeId)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
