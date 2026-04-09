import java.nio.file.Path
import groovy.transform.CompileDynamic

/**
 * Utility functions for testing.
 */
@CompileDynamic
class Utils {

    static Object getRecursiveFileNames(Path fileOrDir, String outputDir) {
        /* groovylint-disable-next-line UnnecessaryGetter */
        if (Path.of(fileOrDir.toString()).toFile().isDirectory()) {
            return fileOrDir.list().collect { file -> getRecursiveFileNames(file, outputDir) }
        }
        return fileOrDir.toString().replace("${outputDir}/", '')
    }

}
