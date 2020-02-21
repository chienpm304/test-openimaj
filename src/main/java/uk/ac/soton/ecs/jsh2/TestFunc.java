package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.colour.Transforms;
import org.openimaj.image.processing.algorithm.GammaCorrection;
import org.openimaj.image.processing.resize.ResizeProcessor;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2dImpl;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static uk.ac.soton.ecs.jsh2.App.*;

public class TestFunc {

    public static void testIntersect() {
        Line2d.IntersectionResult res =
                MyHelper.findIntersection(
                        new Line2d(6, 5, 8, 5),
                        new Line2d(1, 5, 3, 5));
        System.out.println(res.type);
        System.out.println(res.intersectionPoint);
    }

    private static void bruteForceCanny() throws IOException {
        File fin = new File(Constants.WINDOW_DIR);
        File fout = new File(Constants.LINUX_DIR_OUT);
        if (fin.exists() && fin.isDirectory())
            for (final File file : fin.listFiles()) {
                if (file.isFile()) {

                    MBFImage frame = ImageUtilities.readMBF(file);
                    scaleFactor = Constants.STANDARD_WIDTH / (float) frame.getWidth();
                    scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;

                    frame = frame.process(new ResizeProcessor(scaleFactor));
                    width = frame.getWidth();
                    height = frame.getHeight();

                    MBFImage hsv = Transforms.RGB_TO_HSV(frame);

                    for (Constants.CANNY_SIGMA = 3f; Constants.CANNY_SIGMA < 6; Constants.CANNY_SIGMA += 2f) {
                        for (Constants.CANNY_LOW_THRESH = 0.01f; Constants.CANNY_LOW_THRESH < 0.2; Constants.CANNY_LOW_THRESH += 0.03) {
                            for (Constants.CANNY_HIGH_THRESH = Constants.CANNY_LOW_THRESH + 0.05f; Constants.CANNY_HIGH_THRESH < 0.3; Constants.CANNY_HIGH_THRESH += 0.03) {
                                Constants.LINUX_DIR_OUT = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/AZdoc/java/edges/" + Constants.CANNY_SIGMA + "_" + Constants.CANNY_LOW_THRESH + "_" + Constants.CANNY_HIGH_THRESH;
                                fout = new File(Constants.LINUX_DIR_OUT);
                                if (!fout.exists()) fout.mkdirs();
                                FImage edges = applyCannyDetector(hsv.getBand(Constants.S_CHANNEL_ID));
                                ImageUtilities.write(edges, new File(fout.getAbsolutePath() + "/" + file.getName()));

                            }
                        }
                    }
                }
            }

    }

    private static void testAngleGap() {
        System.out.println(calcAngleDiffInDegree(-89, 89));
    }

    private static void testMergeLines() {
        Line2d line1 = new Line2d(1, 1, 1, 6);
        Line2d line2 = new Line2d(1.2f, 4.8f, 1.2f, 6.5f);
        System.out.println("before merge: ");
        System.out.println("keep: " + line1.toString());

        mergeLine(line1, line2);
        System.out.println("after merge: " + line1.toString());
    }

    private static Tetragram detectLSD(File fin, File fout) throws IOException {

        MBFImage frame = ImageUtilities.readMBF(fin);
        scaleFactor = Constants.STANDARD_WIDTH / (float) frame.getWidth();
        scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;

        frame = frame.process(new ResizeProcessor(scaleFactor));//.process(new FGaussianConvolve(2f));

        width = frame.getWidth();
        height = frame.getHeight();

        System.out.println("processing: " + fin.getName() + "...");
        GammaCorrection gc = new GammaCorrection(Constants.GAMMA);
        for (int i = 0; i < frame.numBands(); i++) {
            frame.getBand(i).processInplace(gc);
        }

        Point2dImpl center = new Point2dImpl(width / 2, height / 2);

        List<Line2d> lines = getLinesUsingLineSegmentDetector(frame);

        removeSimilarAndNoiseLines(lines);

        System.out.println("After merge: " + lines.size());

        drawLines(frame, center, lines, RGBColour.GRAY,false);

        List<LineHolder> results = findBounds(lines);


        if (!results.isEmpty())
            drawBound(frame, center, results.get(0).lines, RGBColour.GREEN, RGBColour.YELLOW);


        ImageUtilities.write(frame, new File(fout.getAbsolutePath() + "/out/" + fin.getName()));
        return null;
    }

    private static List<Line2d> getLinesUsingLineSegmentDetector(MBFImage image) {
        LSD lsd = new LSD(image);
        return lsd.getLines();
    }


}
