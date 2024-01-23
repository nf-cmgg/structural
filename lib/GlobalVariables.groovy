class GlobalVariables {
    // When making changes here, make sure to also update the following files: conf/modules.config
    public static List<String> svCallers = ["delly", "manta", "smoove"] //, "gridss"
    public static List<String> cnvCallers = ["qdnaseq", "wisecondorx"]
    public static List<String> repeatsCallers = ["expansionhunter"]

    public static List<String> allCallers = svCallers + cnvCallers + repeatsCallers

    // Callers that need the sex
    public static List<String> sexCallers = ["expansionhunter", "qdnaseq"]

}
