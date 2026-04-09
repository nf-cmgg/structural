import java.nio.file.Path
import groovy.transform.CompileDynamic

/**
 * Utility functions for testing.
 */
@CompileDynamic
class Utils {

    static String getRecursiveFileNames(Path fileOrDir, String outputDir) {
        if (Path.of(fileOrDir.toString()).toFile().directory()) {
            return fileOrDir.list().collect { file -> getRecursiveFileNames(file, outputDir) }
        }
        return fileOrDir.toString().replace("${outputDir}/", '')
    }

}
