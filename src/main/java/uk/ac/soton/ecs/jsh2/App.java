package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.ColourSpace;

import java.io.File;
import java.io.FileWriter;
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
//        testDetectBox();
//        BoundDetector.detectBound(new File(WINDOW_DIR+"/20191227_084126.jpg"), new File(WINDOW_OUT_DIR), true);
//        generateColorImage();
//        writeHex();
        writeRGB();
    }

    private static void writeRGB() throws IOException {

        FileWriter csvWriter = new FileWriter("D:/magic_color/data_rgb.csv");
        csvWriter.append("rin,gin,bin,rout,gout,bout,\n");


        for (int k = 0; k <= 3; k++) {
            MBFImage frameIn = ImageUtilities.readMBF(new File("D:/magic_color/in/" + k + ".jpg"));
            MBFImage frameOut = ImageUtilities.readMBF(new File("D:/magic_color/out/" + k + ".jpg"));

            float[][] rin = frameIn.getBand(0).pixels;
            float[][] gin = frameIn.getBand(1).pixels;
            float[][] bin = frameIn.getBand(2).pixels;


            float[][] rout = frameOut.getBand(0).pixels;
            float[][] gout = frameOut.getBand(1).pixels;
            float[][] bout = frameOut.getBand(2).pixels;


            int height = frameIn.getHeight();
            int width = frameIn.getWidth();

            int redo, greeno, blueo, redi, greeni, bluei;
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    redi = (int) (rin[y][x] * 255);
                    greeni = (int) (gin[y][x] * 255);
                    bluei = (int) (bin[y][x] * 255);

                    redo = (int) (rout[y][x] * 255);
                    greeno = (int) (gout[y][x] * 255);
                    blueo = (int) (bout[y][x] * 255);
                    csvWriter.append(String.join(",",
                            String.valueOf(redi), String.valueOf(greeni), String.valueOf(bluei),
                            String.valueOf(redo), String.valueOf(greeno), String.valueOf(blueo), "\n"));

                }
            }
        }


        csvWriter.flush();
        csvWriter.close();
        System.out.println("done");
    }


    private static void writeHex() throws IOException {

        FileWriter csvWriter = new FileWriter("D:/magic_color/data.csv");
        csvWriter.append("in, out,\n");


        for (int k = 0; k <= 3; k++) {
            MBFImage frameIn = ImageUtilities.readMBF(new File("D:/magic_color/in/" + k + ".jpg"));
            MBFImage frameOut = ImageUtilities.readMBF(new File("D:/magic_color/out/" + k + ".jpg"));

            float[][] rin = frameIn.getBand(0).pixels;
            float[][] gin = frameIn.getBand(1).pixels;
            float[][] bin = frameIn.getBand(2).pixels;


            float[][] rout = frameOut.getBand(0).pixels;
            float[][] gout = frameOut.getBand(1).pixels;
            float[][] bout = frameOut.getBand(2).pixels;


            int height = frameIn.getHeight();
            int width = frameIn.getWidth();

            int red, green, blue, hexIn, hexOut;
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    red = (int) (rin[y][x] * 255);
                    green = (int) (gin[y][x] * 255);
                    blue = (int) (bin[y][x] * 255);
                    hexIn = (red * 65536) + (green * 256) + blue;

                    red = (int) (rout[y][x] * 255);
                    green = (int) (gout[y][x] * 255);
                    blue = (int) (bout[y][x] * 255);
                    hexOut = (red * 65536) + (green * 256) + blue;


                    csvWriter.append(String.join(",", String.valueOf(hexIn), String.valueOf(hexOut), "\n"));// String.valueOf(green), String.valueOf(blue), String.valueOf(hexIn),"\n"));//, green, blue, hex);
                }
            }
        }


        csvWriter.flush();
        csvWriter.close();
        System.out.println("done");
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

    private static void testFilter(File fin, File fout) throws IOException {
        MBFImage frame = ImageUtilities.readMBF(fin);
        // Test convertAbs
        ConvertAbsolute convert = new ConvertAbsolute(1.9f, -90f / 255f);
        for (int i = 0; i < frame.numBands(); i++) {
            frame.getBand(i).processInplace(convert);
        }
        ImageUtilities.write(frame, new File(fout.getAbsolutePath() + "/filter/" + fin.getName()));
    }

    private static void generateColorImage() throws IOException {
        int color = 0x000000;
        int height = 4096;
        int width = 4096;

        int red, green, blue;
        MBFImage image = new MBFImage(width, height, ColourSpace.RGB);
        float[][] rband = image.getBand(0).pixels;
        float[][] gband = image.getBand(1).pixels;
        float[][] bband = image.getBand(2).pixels;

        int index = 0, x, y;

        for (int b = 0; b < 256; b++) {

            for (int r = 0; r < 256; r++) {
                for (int g = 0; g < 256; g++) {
                    y = index / height;
                    x = index % width;
                    rband[y][x] = ImageUtilities.BYTE_TO_FLOAT_LUT[r];
                    gband[y][x] = ImageUtilities.BYTE_TO_FLOAT_LUT[g];
                    bband[y][x] = ImageUtilities.BYTE_TO_FLOAT_LUT[b];
                    index++;
                }
            }
        }

//        for(int y = 0; y < height; ++y)
//            for (int x = 0; x < width; ++x) {
//                red = color >> 16 & 255;
//                green = color >> 8 & 255;
//                blue = color & 255;
//                rband[y][x] = ImageUtilities.BYTE_TO_FLOAT_LUT[red];
//                gband[y][x] = ImageUtilities.BYTE_TO_FLOAT_LUT[green];
//                bband[y][x] = ImageUtilities.BYTE_TO_FLOAT_LUT[blue];
//                color+=1;
//            }
//        }


        ImageUtilities.write(image, new File("D:/detect/input/AZdoc/1.jpg"));
        System.out.printf("done");
    }


}
