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
//        generateColorImage();
//        writeHex();
//        writeRGB();
//        compare(new File("D:/magic_color/cropped/2.jpg"),new File("D:/magic_color/cropped/3.jpg"));
        genLUT();
    }

    private static void genLUT() throws IOException {
        // Additional Info:
        // Lookup texture is organised as 8x8 quads of 64x64 pixels representing all possible RGB colors:

        int[] data = new int[512*512];
        int color, red, green, blue;
        for (int by = 0; by < 8; by++) {
            for (int bx = 0; bx < 8; bx++) {
                for (int g = 0; g < 64; g++) {
                    for (int r = 0; r < 64; r++) {
//                        image.setPixel(r + bx * 64, g + by * 64, qRgb((int)(r * 255.0 / 63.0 + 0.5),
//                                                                            (int)(g * 255.0 / 63.0 + 0.5),
//                                                                            (int)((bx + by * 8.0) * 255.0 / 63.0 + 0.5)));

                        red = (int) (r * 255.0 / 63.0 + 0.5);
                        green = (int) (g * 255.0 / 63.0 + 0.5);
                        blue = (int) ((bx + by * 8.0) * 255.0 / 63.0 + 0.5);
//                        System.out.println(red+", "+green+", "+blue);
                        color = packRGB(red, green, blue);
                        data[(r + bx * 64) + (g + by * 64) * 512] = color;
                    }
                }
            }
        }
        MBFImage image = new MBFImage(data,512, 512);
        ImageUtilities.write(image, new File("D:/magic_color/lut_512x512.jpg"));
    }

    private static int packRGB(int R, int G, int B) {
        return (R*65536)+(G*256)+B;
    }

    private static void writeRGB() throws IOException {

        FileWriter csvWriter = new FileWriter("D:/magic_color/data_rgb.csv");
        csvWriter.append("rin,gin,bin,rout,gout,bout,\n");


//        for (int k = 0; k <= 3; k++) {
        MBFImage frameIn = ImageUtilities.readMBF(new File("D:/magic_color/in/full.jpg"));
        MBFImage frameOut = ImageUtilities.readMBF(new File("D:/magic_color/out/full.jpg"));

//            DisplayUtilities.display(frameIn);
//            DisplayUtilities.display(frameOut);

        float[][] rin = frameIn.getBand(0).pixels;
        float[][] gin = frameIn.getBand(1).pixels;
        float[][] bin = frameIn.getBand(2).pixels;


        float[][] rout = frameOut.getBand(0).pixels;
        float[][] gout = frameOut.getBand(1).pixels;
        float[][] bout = frameOut.getBand(2).pixels;

        if (frameIn.getWidth() != frameOut.getWidth() || frameIn.getHeight() != frameOut.getHeight())
            throw new RuntimeException("input and output not match!");

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
//        }


        csvWriter.flush();
        csvWriter.close();
        System.out.println("done");
    }


    private static void writeHex() throws IOException {

        FileWriter csvWriter = new FileWriter("D:/magic_color/data_hex_full.csv");
        csvWriter.append("in, out,\n");


//        for (int k = 0; k <= 0; k++) {
        MBFImage frameIn = ImageUtilities.readMBF(new File("D:/magic_color/in/full.jpg"));
        MBFImage frameOut = ImageUtilities.readMBF(new File("D:/magic_color/out/full.jpg"));

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
//        }


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
        int height = 4096;
        int width = 4096;

        int red, green, blue;
        MBFImage image = new MBFImage(width, height, ColourSpace.RGB);
        float[][] rband = image.getBand(0).pixels;
        float[][] gband = image.getBand(1).pixels;
        float[][] bband = image.getBand(2).pixels;

        int index = 0, x, y;

        for (int r = 0; r < 256; r++) {
            for (int g = 0; g < 256; g++) {
                for (int b = 0; b < 256; b++) {
                    y = index / height;
                    x = index % width;
                    rband[y][x] = ImageUtilities.BYTE_TO_FLOAT_LUT[r];
                    gband[y][x] = ImageUtilities.BYTE_TO_FLOAT_LUT[g];
                    bband[y][x] = ImageUtilities.BYTE_TO_FLOAT_LUT[b];
                    index++;
                }
            }
        }


        ImageUtilities.write(image, new File("D:/magic_color/2.jpg"));
        System.out.printf("done");
    }


    static boolean compare(File file1, File file2) throws IOException {
        MBFImage image1 = ImageUtilities.readMBF(file1);
        MBFImage image2 = ImageUtilities.readMBF(file2);

        FileWriter csvWriter = new FileWriter("D:/magic_color/compare_" + file1.getName() + "_" + file2.getName() + ".csv");
        csvWriter.append("rin,gin,bin,rout,gout,bout,\n");


        if (image1.getWidth() != image2.getWidth()
                || image1.getHeight() != image2.getHeight())
            return false;

        int height = image1.getHeight();
        int width = image1.getWidth();
        float[][] rin = image1.getBand(0).pixels;
        float[][] gin = image1.getBand(1).pixels;
        float[][] bin = image1.getBand(2).pixels;


        float[][] rout = image2.getBand(0).pixels;
        float[][] gout = image2.getBand(1).pixels;
        float[][] bout = image2.getBand(2).pixels;

        int red, green, blue, hexIn, hexOut;
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
        return true;
    }

}
