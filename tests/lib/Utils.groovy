import java.nio.file.Path
import groovy.transform.CompileStatic
import nextflow.Nextflow

/**
 * Utility functions for testing.
 */
@CompileStatic
class Utils {

    static String getRecursiveFileNames(Path fileOrDir, String outputDir) {
        if (Nextflow.file(fileOrDir.toString()).directory) {
            return fileOrDir.list().collect { file -> getRecursiveFileNames(file, outputDir) }
        }
        return fileOrDir.toString().replace("${outputDir}/", '')
    }

}
