class DynamicContainers {
    static String apiUrl = "https://api.biocontainers.pro/ga4gh/trs/v2"

    public static String findLatest(String containerName, Double version, String type="docker") {
        this.findLatest(containerName, version as String, type)
    }

    public static String findLatest(String containerName, Float version, String type="docker") {
        this.findLatest(containerName, version as String, type)
    }

    public static String findLatest(String containerName, Integer version, String type="docker") {
        this.findLatest(containerName, version as String, type)
    }

    public static String findLatest(String containerName, String version, String type="docker") {
        type = type ? type.toLowerCase() : "docker"
        def String url = "${this.apiUrl}/tools/${containerName}/versions/${containerName}-${version}"
        def biocontainersGet = new URL(url).openConnection()
        biocontainersGet.setRequestProperty("accept", "application/json")
        def getRC = biocontainersGet.getResponseCode()

        if (getRC != 200) { return null }

        def parser = new groovy.json.JsonSlurper()
        def jsonResponse = parser.parseText(biocontainersGet.getInputStream().getText())

        def ArrayList images = jsonResponse["images"].findAll { it["image_type"].toLowerCase() == type }
        if(images.size() > 1 && type != "conda") {
            return this.getLatestImage(images)
        } 
        else if(images.size() == 0 && type == "singularity") {
            images = jsonResponse["images"].findAll { it["image_type"].toLowerCase() == "docker" }
            return this.getLatestImage(images)
        }
        else {
            return images[0]["image_name"]
        }

    }

    static String getLatestImage(ArrayList images) {
        def Date latestTime
        def String latestImage
        for(image : images) {
            def Integer year = image["updated"].split("-")[0] as Integer
            def Integer month = image["updated"].split("-")[1] as Integer
            def Integer day = image["updated"].split("-")[2].split("T")[0] as Integer
            def Integer hour = image["updated"].split("T")[1].split(":")[0] as Integer
            def Integer minutes = image["updated"].split("T")[1].split(":")[1] as Integer
            def Integer seconds = image["updated"].split("T")[1].split(":")[2].split("Z")[0] as Integer
            def Date time = new Date(year, month, day, hour, minutes, seconds)
            if (time > latestTime) {
                latestTime = time
                latestImage = image["image_name"]
            }
        }
        return latestImage
    }
}
