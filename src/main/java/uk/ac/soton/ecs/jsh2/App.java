package uk.ac.soton.ecs.jsh2;

import java.io.File;
import java.io.IOException;

/**
 * OpenIMAJ Hello world!
 */
public class App {

    public static final String WINDOW_DIR = "D:/detect/input/AZdoc/new";
    public static final String WINDOW_OUT_DIR = "D:/detect/input/AZdoc";
    public static final String LINUX_DIR_IN = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/AZdoc/new";
    public static final String LINUX_DIR_OUT = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/AZdoc/";

    public static void main(String[] args) throws IOException {
        testDetectBox();
    }

    private static void testDetectBox() throws IOException {
        File fin = new File(WINDOW_DIR);
        File fout = new File(WINDOW_OUT_DIR);
        if (fin.exists() && fin.isDirectory())
            for (final File file : fin.listFiles()) {
                if (file.isFile())
                    BoundDetector.detectBound(file, fout, true);
            }
    }

}
