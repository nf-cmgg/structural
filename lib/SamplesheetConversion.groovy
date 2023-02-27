import groovyx.gpars.dataflow.DataflowBroadcast
import groovy.json.JsonSlurper
import org.yaml.snakeyaml.Yaml
import java.nio.file.Path
import nextflow.Nextflow
import nextflow.Channel


class SamplesheetConversion {
    public static DataflowBroadcast convert(
        Path samplesheetFile,
        Path schemaFile
    ) {

        def Map schema = (Map) new JsonSlurper().parseText(schemaFile.text)
        def Map schemaFields = schema.get("properties")
        def ArrayList allFields = schemaFields.keySet().collect()
        def ArrayList requiredFields = schema.get("required")

        def String fileType = getFileType(samplesheetFile)
        def String delimiter = fileType == "csv" ? "," : fileType == "tsv" ? "\t" : null
        def DataflowBroadcast samplesheet

        if(fileType == "yaml"){
            samplesheet = Channel.fromList(new Yaml().load((samplesheetFile.text)))
        }
        else {
            samplesheet = Channel.value(samplesheetFile).splitCsv(header:true, strip:true, sep:delimiter)
        }

        // Field checks + returning the channels
        def Map uniques = [:]
        def Boolean headerCheck = true
        def Integer sampleCount = 0

        return samplesheet.map({ entry ->

            sampleCount++

            // Check the header once for CSV/TSV and for every sample for YAML
            if(headerCheck) {
                def ArrayList entryKeys = entry.keySet().collect()
                def ArrayList differences = allFields.plus(entryKeys)
                differences.removeAll(allFields.intersect(entryKeys))

                def String yamlInfo = fileType == "yaml" ? " for sample ${sampleCount}." : ""

                def ArrayList samplesheetDifferences = entryKeys.intersect(differences)
                if(samplesheetDifferences.size > 0) {
                    throw new Exception("[Samplesheet Error] The samplesheet contains following unwanted field(s): ${samplesheetDifferences}${yamlInfo}")
                }

                def ArrayList requiredDifferences = requiredFields.intersect(differences)
                if(requiredDifferences.size > 0) {
                    throw new Exception("[Samplesheet Error] The samplesheet requires '${requiredFields.join(",")}' as header field(s), but is missing these: ${requiredDifferences}${yamlInfo}")
                }

                if(fileType in ["csv", "tsv"]) {
                    headerCheck = false
                }
            }

            // Check required dependencies
            def Map dependencies = schema.get("dependentRequired")
            if(dependencies) {
                for( dependency in dependencies ){
                    if(entry[dependency.key] != "" && entry[dependency.key]) {
                        def ArrayList missingValues = []
                        for( value in dependency.value ){
                            if(entry[value] == "" || !entry[value]) {
                                missingValues.add(value)
                            }
                        }
                        if (missingValues) {
                            throw new Exception("[Samplesheet Error] ${dependency.value} field(s) should be defined when '${dependency.key}' is specified, but  the field(s) ${missingValues} are/is not defined.")
                        }
                    }
                }
            }

            def Map meta = [:]
            def ArrayList output = []

            for( field : schemaFields ){
                def String key = field.key
                def String regexPattern = field.value.pattern && field.value.pattern != '' ? field.value.pattern : '^.*$'
                def String metaNames = field.value.meta

                def String input = entry[key]

                if((input == null || input == "") && key in requiredFields){
                    throw new Exception("[Samplesheet Error] Sample ${sampleCount} does not contain an input for required field '${key}'.")
                }
                else if(!(input ==~ regexPattern) && input != '' && input) {
                    throw new Exception("[Samplesheet Error] The '${key}' value for sample ${sampleCount} does not match the pattern '${regexPattern}'.")
                }
                else if(field.value.unique){
                    if(!(key in uniques)){uniques[key] = []}
                    if(input in uniques[key] && input){
                        throw new Exception("[Samplesheet Error] The '${key}' value needs to be unique. '${input}' was found twice in the samplesheet.")
                    }
                    uniques[key].add(input)
                }

                if(metaNames) {
                    for(name : metaNames.tokenize(',')) {
                        meta[name] = (input != '' && input) ? checkAndTransform(input, field, sampleCount) : field.value.default ?: null
                    }
                }
                else {
                    def inputFile = (input != '' && input) ? checkAndTransform(input, field, sampleCount) : field.value.default ? checkAndTransform(field.value.default, field, sampleCount) : []
                    output.add(inputFile)
                }
            }
            output.add(0, meta)
            return output
        })

    }

    // Function to infer the file type of the samplesheet
    private static String getFileType(
        Path samplesheetFile
    ) {
        def String extension = samplesheetFile.getExtension()
        if (extension in ["csv", "tsv", "yml", "yaml"]) {
            return extension == "yml" ? "yaml" : extension
        }

        def String header = getHeader(samplesheetFile)

        def Integer commaCount = header.count(",")
        def Integer tabCount = header.count("\t")

        if ( commaCount == tabCount ){
            throw new Exception("[Samplesheet Error] Could not derive file type from ${samplesheetFile}. Please specify the file extension (CSV, TSV, YML and YAML are supported).")
        }
        if ( commaCount > tabCount ){
            return "csv"
        }
        else {
            return "tsv"
        }
    }

    // Function to get the header from a CSV or TSV file
    private static String getHeader(
        Path samplesheetFile
    ) {
        def String header
        samplesheetFile.withReader { header = it.readLine() }
        return header
    }

    // Function to check and transform an input field from the samplesheet
    private static checkAndTransform(
        input,
        LinkedHashMap$Entry field,
        Integer sampleCount
    ) {
        def String type = field.value.type
        def String format = field.value.format
        def String key = field.key

        def ArrayList supportedTypes = ["string", "integer", "boolean"]
        if(!(type in supportedTypes)) {
            throw new Exception("[Samplesheet Schema Error] The type '${type}' specified for ${key} is not supported. Please specify one of these instead: ${supportedTypes}")
        }

        if(type == "string" || !type) {
            def ArrayList supportedFormats = ["file-path", "directory-path"]
            if(!(format in supportedFormats) && format) {
                throw new Exception("[Samplesheet Schema Error] The string format '${format}' specified for ${key} is not supported. Please specify one of these instead: ${supportedFormats} or don't supply a format for a simple string.")
            }
            if(format == "file-path" || format =="directory-path") {
                def inputFile = Nextflow.file(input)
                if(!inputFile.exists()){
                    throw new Exception("[Samplesheet Error] The '${key}' file or directory (${input}) for sample ${sampleCount} does not exist.")
                }
                return inputFile
            }
            else {
                return input as String
            }
        }
        else if(type == "integer") {
            try {
                return input as Integer
            } catch(java.lang.NumberFormatException e) {
                throw new Exception("[Samplesheet Error] The '${key}' value (${input}) for sample ${sampleCount} is not a valid integer.")
            }
        }
        else if(type == "boolean") {
            if(input.toLowerCase() == "true") {
                return true
            }
            else if(input.toLowerCase() == "false") {
                return false
            }
            else {
                throw new Exception("[Samplesheet Error] The '${key}' value (${input}) for sample ${sampleCount} is not a valid boolean.")
            }
        }
    }
}
