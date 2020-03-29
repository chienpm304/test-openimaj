package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;

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
//        BoundDetector.detectBound(new File(WINDOW_DIR+"/20191227_084121.jpg"), new File(WINDOW_OUT_DIR), true);
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


    private static void testFilter(File fin, File fout) throws IOException{
        MBFImage frame = ImageUtilities.readMBF(fin);
        // Test convertAbs
        ConvertAbsolute convert = new ConvertAbsolute(1.9f, -90f/255f);
        for(int i = 0; i < frame.numBands(); i++){
            frame.getBand(i).processInplace(convert);
        }
        ImageUtilities.write(frame, new File(fout.getAbsolutePath() + "/filter/" + fin.getName()));
    }
}
