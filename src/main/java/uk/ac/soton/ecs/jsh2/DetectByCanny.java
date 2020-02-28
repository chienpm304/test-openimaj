package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.processing.algorithm.EqualisationProcessor;
import org.openimaj.image.processing.algorithm.GammaCorrection;
import org.openimaj.image.processing.convolution.FGaussianConvolve;
import org.openimaj.image.processing.edges.CannyEdgeDetector;
import org.openimaj.image.processing.resize.ResizeProcessor;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static uk.ac.soton.ecs.jsh2.App.*;

public class DetectByCanny {

    public static float _CANNY_LOW_THRESH = 0.05f;
    public static float _CANNY_HIGH_THRESH = 0.25f;
    public static float _CANNY_SIGMA = 5f;

    public static final double _HOUGH_LINE_RHO = 1;
    public static final double _HOUGH_LINE_THETA = Math.PI / 180d;;

    public static int _HOUGH_LINE_MAX_LINE_GAP = 30;
    public static int _HOUGH_LINE_THRESHOLD = 30;
    public static int _HOUGH_MIN_LINE_LENGTH = 20;
    public static int _HOUGH_MAX_NUM_LINE = 1000;

    public static final double _GAMMA = 2.2d;
    public static final float _GAUSSIAN_BLUR_SIGMA = 2.5f;


    private static float scaleFactor;
    private static int width, height;

    public static Tetragram detectWithCanny(File fin, File fout) throws IOException {

        MBFImage frame = ImageUtilities.readMBF(fin);

        initSizeAndResizeImage(frame);

        System.out.println("processing: " + fin.getName() + "...");

        enhanceInputImage(frame);
        ImageUtilities.write(frame, new File(fout.getAbsolutePath() + "/enhanced/" + fin.getName()));


        Point2dImpl center = new Point2dImpl(width / 2, height / 2);

        FImage edges = applyCannyDetector(frame.flattenMax());
        ImageUtilities.write(edges, new File(fout.getAbsolutePath() + "/edged/" + fin.getName()));

        List<Line2d> lines = getLinesUsingHoughTransformP(edges);

        MBFImage tmp = frame.clone();
        drawLines(tmp, center, lines, RGBColour.GREEN, false);
        ImageUtilities.write(tmp, new File(fout.getAbsolutePath() + "/raw_detected/" + fin.getName()));
//
//
//        System.out.println("Before merge: " + lines.size());
//
//        removeNoiseLines(lines);
//
//        System.out.println("After merge: " + lines.size());
//
//        MBFImage merged = frame.clone();
//        drawLines(merged, center, lines, RGBColour.GREEN, true);
//        ImageUtilities.write(merged, new File(fout.getAbsolutePath() + "/merged/" + fin.getName()));

        Tetragram bound = findBounds2(lines);

        if (bound!=null) {
            drawBound(frame, center, bound.toLineList(), RGBColour.GREEN, RGBColour.YELLOW);
        }

        ImageUtilities.write(frame, new File(fout.getAbsolutePath() + "/out/" + fin.getName()));
        return null;
    }

    private static void removeNoiseLines(List<Line2d> lines) {
        // do nothing
    }

    private static void enhanceInputImage(MBFImage frame) {
        // process gamma correction
        GammaCorrection gc = new GammaCorrection(_GAMMA);
        for (int i = 0; i < frame.numBands(); i++) {
            frame.getBand(i).processInplace(gc);

        }
        frame.processInplace(new FGaussianConvolve(_GAUSSIAN_BLUR_SIGMA));
    }

    private static void initSizeAndResizeImage(MBFImage frame) {
        if(frame.getWidth() > frame.getHeight())
            scaleFactor = Constants.STANDARD_WIDTH / (float) frame.getHeight();
        else
            scaleFactor = Constants.STANDARD_WIDTH / (float) frame.getWidth();

        scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;

        if(scaleFactor!=1.0f) {
            frame.processInplace(new ResizeProcessor(scaleFactor));
        }

        width = frame.getWidth();
        height = frame.getHeight();
    }

    static FImage applyCannyDetector(FImage grey) {
        CannyEdgeDetector canny = new CannyEdgeDetector(_CANNY_LOW_THRESH, _CANNY_HIGH_THRESH, _CANNY_SIGMA);
        canny.processImage(grey);
        return grey;
    }

    private static List<Line2d> getLinesUsingHoughTransformP(FImage image) {
        List<Point2d> points = new ArrayList<>();
        for (int i = 0; i < image.width; i++)
            for (int j = 0; j < image.height; j++) {
                if (image.getPixel(i, j) > 0.5f)
                    points.add(new Point2dImpl(i, j));
            }
        HoughLinesP ht = new HoughLinesP(
                points,
                width,
                height,
                _HOUGH_LINE_RHO,
                _HOUGH_LINE_THETA,
                _HOUGH_LINE_THRESHOLD,
                _HOUGH_LINE_MAX_LINE_GAP,
                _HOUGH_MIN_LINE_LENGTH,
                _HOUGH_MAX_NUM_LINE);
        return ht.getLines();
    }

}
