package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.colour.Transforms;
import org.openimaj.image.processing.algorithm.GammaCorrection;
import org.openimaj.image.processing.convolution.FGaussianConvolve;
import org.openimaj.image.processing.edges.CannyEdgeDetector;
import org.openimaj.image.processing.resize.ResizeProcessor;
import org.openimaj.image.typography.hershey.HersheyFont;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;

import static java.lang.Math.PI;

/**
 * OpenIMAJ Hello world!
 */
public class App {
    public static final float THRESHOLD_BIN_INV = 0.07133f;
    public static final float STANDARD_WIDTH = 720;
    public static final int H_CHANNEL_ID = 0;
    public static final int S_CHANNEL_ID = 1;
    public static final int V_CHANNEL_ID = 2;

    public static final String WINDOW_DIR = "D:/detect/input/AZdoc/s";
    public static final double GAMMA = 2d;


    public static String WINDOW_OUT_DIR = "D:/detect/input/AZdoc";

    private static final String LINUX_DIR_IN = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/AZdoc/in";
    private static String LINUX_DIR_OUT = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/AZdoc/java";
    public static final int THRESHOLD_STEP = 50;
    public static final float NEW_THRESHOLD = 0.1f;


    public static float scaleFactor = 1.0f;
    private static int height;
    private static int width;

    public static final float GAUSSIAN_BLUR_SIGMA = 3f;

    public static float CANNY_LOW_THRESH = 0.01f;
    public static float CANNY_HIGH_THRESH = 0.05f;
    public static float CANNY_SIGMA = 5f;

    public static final double HOUGH_LINE_RHO = 1;
    public static final double HOUGH_LINE_THETA = Math.PI / 180d;

    public static int HOUGHLINE_THRESHOLD = 80;
    public static int HOUGH_LINE_LENGTH = 150;

    private static final int LINE_GAP_REMOVAL = 20;
    private static final int BOUNDING_GAP_REMOVAL = 3;

    public static final int HISTOGRAM_NBINS = 64;


    public static void testIntersect() {
        Line2d.IntersectionResult res =
                MyHelper.findIntersection(
                        new Line2d(6, 5, 8, 5),
                        new Line2d(1, 5, 3, 5));
        System.out.println(res.type);
        System.out.println(res.intersectionPoint);
    }

    public static void main(String[] args) throws IOException {
        testDetectBox();
//        testIntersect();
//        testMergeLines();
//        testAngleGap();
//        detectWithCanny(new File(LINUX_DIR_IN+"/8.jpg"), new File(LINUX_DIR_OUT));
//        bruteForceCanny();
    }

    private static void bruteForceCanny() throws IOException {
        File fin = new File(WINDOW_DIR);
        File fout = new File(LINUX_DIR_OUT);
        if (fin.exists() && fin.isDirectory())
            for (final File file : fin.listFiles()) {
                if (file.isFile()) {

                    MBFImage frame = ImageUtilities.readMBF(file);
                    scaleFactor = STANDARD_WIDTH / (float) frame.getWidth();
                    scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;

                    frame = frame.process(new ResizeProcessor(scaleFactor));
                    width = frame.getWidth();
                    height = frame.getHeight();

                    MBFImage hsv = Transforms.RGB_TO_HSV(frame);

                    for (CANNY_SIGMA = 3f; CANNY_SIGMA < 6; CANNY_SIGMA += 2f) {
                        for (CANNY_LOW_THRESH = 0.01f; CANNY_LOW_THRESH < 0.2; CANNY_LOW_THRESH += 0.03) {
                            for (CANNY_HIGH_THRESH = CANNY_LOW_THRESH + 0.05f; CANNY_HIGH_THRESH < 0.3; CANNY_HIGH_THRESH += 0.03) {
                                LINUX_DIR_OUT = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/AZdoc/java/edges/" + CANNY_SIGMA + "_" + CANNY_LOW_THRESH + "_" + CANNY_HIGH_THRESH;
                                fout = new File(LINUX_DIR_OUT);
                                if (!fout.exists()) fout.mkdirs();
                                FImage edges = applyCannyDetector(hsv.getBand(S_CHANNEL_ID), file, fout);
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

    private static void testDetectBox() throws IOException {
        File fin = new File(WINDOW_DIR);
        File fout = new File(WINDOW_OUT_DIR);
        if (fin.exists() && fin.isDirectory())
            for (final File file : fin.listFiles()) {
                if (file.isFile())
                    detectLSD(file, fout);
            }
    }
    private static Tetragram detectLSD(File fin, File fout) throws IOException {

        MBFImage frame = ImageUtilities.readMBF(fin);
        scaleFactor = STANDARD_WIDTH / (float) frame.getWidth();
        scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;

        frame = frame.process(new ResizeProcessor(scaleFactor));//.process(new FGaussianConvolve(2f));

        width = frame.getWidth();
        height = frame.getHeight();

        System.out.println("processing: " + fin.getName() + "...");
        GammaCorrection gc = new GammaCorrection(GAMMA);
        for (int i = 0; i < frame.numBands(); i++) {
            frame.getBand(i).processInplace(gc);
        }

        Point2dImpl center = new Point2dImpl(width / 2, height / 2);

        List<Line2d> lines = getLinesUsingLineSegmentDetector(frame);

        removeSimilarAndNoiseLines(lines);

        System.out.println("After merge: " + lines.size());

        drawLines(frame, center, lines, RGBColour.GRAY);

        List<LineHolder> results = findBounds(lines);


        if (!results.isEmpty())
            drawBound(frame, center, results.get(0).lines, RGBColour.GREEN, RGBColour.YELLOW);


        ImageUtilities.write(frame, new File(fout.getAbsolutePath() + "/out/" + fin.getName()));
        return null;
    }

    private static Tetragram detectWithCanny(File fin, File fout) throws IOException {

        MBFImage frame = ImageUtilities.readMBF(fin);
        scaleFactor = STANDARD_WIDTH / (float) frame.getWidth();
        scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;

        frame = frame.process(new ResizeProcessor(scaleFactor));
//        ImageUtilities.write(frame, new File(fout.getAbsolutePath()+"/gc/"+fin.getName()));

        width = frame.getWidth();
        height = frame.getHeight();

        System.out.println("processing: " + fin.getName() + "...");
        GammaCorrection gc = new GammaCorrection(GAMMA);
        for (int i = 0; i < frame.numBands(); i++) {
            frame.getBand(i).processInplace(gc);
        }

        MBFImage hsv = Transforms.RGB_TO_HSV(frame);
//        hsv.getBand(0).fill(0);
//        hsv.getBand(2).fill(0);

//
//        ImageUtilities.write(hsv.getBand(0), new File(fout.getAbsolutePath()+"/h/"+fin.getName()));
//        ImageUtilities.write(hsv.getBand(1), new File(fout.getAbsolutePath()+"/s/"+fin.getName()));
//        ImageUtilities.write(hsv.getBand(2), new File(fout.getAbsolutePath()+"/v/"+fin.getName()));
//        if(1==1) return null;
        Point2dImpl center = new Point2dImpl(width / 2, height / 2);


        FImage edges = applyCannyDetector(hsv.getBand(S_CHANNEL_ID), fin, fout);
        ImageUtilities.write(edges, new File(fout.getAbsolutePath() + "/edged/" + fin.getName()));

        List<Line2d> lines = getLinesUsingHoughTransformP(edges);

        System.out.println("Before merge: "+lines.size());
        removeSimilarAndNoiseLines(lines);

        System.out.println("After merge: " + lines.size());

        drawLines(frame, center, lines, RGBColour.GRAY);

        List<LineHolder> results = findBounds(lines);


//        drawLines(frame,center, lines, RGBColour.BLUE);

        if (!results.isEmpty())
            drawBound(frame, center, results.get(0).lines, RGBColour.GREEN, RGBColour.YELLOW);


        ImageUtilities.write(frame, new File(fout.getAbsolutePath() + "/out/" + fin.getName()));
        return null;
    }


    private static List<LineHolder> findBounds(List<Line2d> lines) {

        ArrayList<LineHolder> res = new ArrayList<>();

        Line2d base1, base2, fit1, fit2;
        double b1, b2, f1, f2;
        int n = lines.size();
        boolean addDummyLine = false;
        int angle_step = 5;
        boolean detected = false;

        while (res.isEmpty() && angle_step <= 45) {
            for (int i = 0; i < n - 1; i++) {
                base1 = lines.get(i);
                b1 = getAngleInDegree(base1);
                for (int j = i + 1; j < n; j++) {
                    base2 = lines.get(j);
                    b2 = getAngleInDegree(base2);
                    detected = false;
                    if (calcAngleDiffInDegree(b1, b2) <= angle_step) {

                        for (int k = 0; k < lines.size() - 1; k++) {
                            if (k == i || k == j) continue;

                            fit1 = lines.get(k);

                            for (int l = k + 1; l < lines.size(); l++) {
                                if (l == i || l == j) continue;
                                fit2 = lines.get(l);

                                LineHolder lh = considerToAddBound(base1, base2, fit1, fit2, angle_step);
                                if (lh != null) {
                                    res.add(lh);
                                    detected = true;
                                }
                            }
                        }

                        if (!detected) {
                            Line2d tmp1 = base1.clone();
                            Line2d tmp2 = base2.clone();
                            // added 2 dummy lines
                            if (b1 < 45) {
                                sortByXAxis(tmp1);
                                sortByXAxis(tmp2);
                            } else {
                                sortByYAxis(tmp1);
                                sortByYAxis(tmp2);
                            }

                            fit1 = new Line2d(base1.begin, base2.begin);
                            fit2 = new Line2d(base1.end, base2.end);

                            LineHolder lh = considerToAddBound(base1, base2, fit1, fit2, angle_step);
                            if (lh != null) {
                                res.add(lh);
                                detected = true;
                            }
                        }
                    }
                }
            }
            angle_step += 5;
        }
        Collections.sort(res);
        System.out.println("Detected " + res.size() + " bounds");
//        for (LineHolder lh : res) {
//            System.out.println(lh.toString());
//        }


        return res;
    }

    private static LineHolder considerToAddBound(Line2d base1, Line2d base2, Line2d fit1, Line2d fit2, double angle_step) {
        double f1 = getAngleInDegree(fit1);
        double f2 = getAngleInDegree(fit2);

        double b1 = getAngleInDegree(base1);
        double b2 = getAngleInDegree(base2);

        double bases_distance = (base1.distanceToLine(base2.begin) + base1.distanceToLine(base2.end)) / 2d;

        if (calcAngleDiffInDegree(b1, f1) >= 90 - angle_step * 2 && calcAngleDiffInDegree(b2, f1) >= 90 - angle_step * 2) {
            if (calcAngleDiffInDegree(b1, f2) >= 90 - angle_step * 2 && calcAngleDiffInDegree(b2, f2) >= 90 - angle_step * 2) {
                LineHolder lh = new LineHolder(base1, fit1, base2, fit2);
                double base1_fit1 = Math.min(fit1.distanceToLine(base1.begin), fit1.distanceToLine(base1.end));
                double base1_fit2 = Math.min(fit2.distanceToLine(base1.begin), fit2.distanceToLine(base1.end));
                double base2_fit1 = Math.min(fit1.distanceToLine(base2.begin), fit1.distanceToLine(base2.end));
                double base2_fit2 = Math.min(fit2.distanceToLine(base2.begin), fit2.distanceToLine(base2.end));

                double fit1_base1 = Math.min(base1.distanceToLine(fit1.begin), base1.distanceToLine(fit1.end));
                double fit1_base2 = Math.min(base2.distanceToLine(fit1.begin), base2.distanceToLine(fit1.end));
                double fit2_base1 = Math.min(base1.distanceToLine(fit2.begin), base1.distanceToLine(fit2.end));
                double fit2_base2 = Math.min(base2.distanceToLine(fit2.begin), base2.distanceToLine(fit2.end));

                double fits_distance = (fit1.distanceToLine(fit2.begin) + fit1.distanceToLine(fit2.end)) / 2d;


                lh.rank =
                        base1_fit1 + base1_fit2 + base2_fit1 + base2_fit2
                                - fits_distance - bases_distance +
                                fit1_base1 + fit1_base2 + fit2_base1 + fit2_base2;

//                                    lh.rank = Math.min(base1_fit1, fit1_base1) +
//                                            Math.min(base1_fit2, fit1_base2) +
//                                            Math.min(base2_fit1, fit2_base1) +
//                                            Math.min(base2_fit2, fit2_base2);

                return lh;
            }
        }
        return null;
    }

    private static double getAngleInDegree(Line2d line) {
        return line.calculateHorizontalAngle() * 180 / PI;
    }

    private static double calcAngleDiffInDegree(double a1, double a2) {
        double gap = 0;
        if (a1 * a2 >= 0)
            return Math.abs(a1 - a2);

        if (a1 < 0) gap = a2 - a1;
        else gap = a1 - a2;

        if (gap >= 90)
            gap = 180 - gap;

        return gap;
    }

    private static List<Line2d> getLinesUsingHoughTransformP(FImage image) {
        List<Point2d> points = new ArrayList<>();
        for (int i = 0; i < image.width; i++)
            for (int j = 0; j < image.height; j++) {
                if (image.getPixel(i, j) > 0.5f)
                    points.add(new Point2dImpl(i, j));
            }
        HoughLinesP ht = new HoughLinesP(points, width, height, HOUGH_LINE_RHO, HOUGH_LINE_THETA, HOUGHLINE_THRESHOLD, HOUGH_LINE_LENGTH, 200);
        return ht.getLines();
    }

    private static List<Line2d> getLinesUsingLineSegmentDetector(MBFImage image){
        LSD lsd = new LSD(image);
        return lsd.getLines();
    }

    private static List<Line2d> removeSimilarAndNoiseLines(List<Line2d> lines) {
        for (int i = 0; i < lines.size(); i++) {
            Line2d l = lines.get(i);

            if (l.begin.getX() < BOUNDING_GAP_REMOVAL && l.end.getX() < BOUNDING_GAP_REMOVAL //LEFT

                    || l.begin.getY() < BOUNDING_GAP_REMOVAL && l.end.getY() < BOUNDING_GAP_REMOVAL //TOP

                    || l.begin.getX() > width - BOUNDING_GAP_REMOVAL && l.end.getX() > width - BOUNDING_GAP_REMOVAL // RIGHT

                    || l.begin.getY() > height - BOUNDING_GAP_REMOVAL && l.end.getY() > height - BOUNDING_GAP_REMOVAL) {

                lines.remove(i);
                i--;
            }
        }

        Line2d l1, l2;
        Line2d kpLine, rmLine;
        for (int i = 0; i < lines.size() - 1; i++) {
            l1 = lines.get(i);
            for (int j = i + 1; j < lines.size(); j++) {
                l2 = lines.get(j);
                if (l1 != l2
                        && l2.distanceToLine(l1.begin) < LINE_GAP_REMOVAL
                        && l2.distanceToLine(l1.end) < LINE_GAP_REMOVAL
                        && l1.distanceToLine(l2.begin) < LINE_GAP_REMOVAL
                        && l1.distanceToLine(l2.end) < LINE_GAP_REMOVAL
                ) {
                    if (l1.calculateLength() < l2.calculateLength()) { // keep line i
                        Collections.swap(lines, i, j);
                        rmLine = lines.remove(j);
                        kpLine = lines.get(i);
                        j = i;

                    } else {
                        rmLine = lines.remove(j);
                        kpLine = lines.get(i);
                        j--;
                    }

                    //merge rmLine to kpLine
                    mergeLine(kpLine, rmLine);
                    lines.set(i, kpLine);
                    l1 = lines.get(i);
                }
            }
        }
        for (int i = 0; i < lines.size(); i++) {
            Line2d l = lines.get(i);

            if (l.calculateLength() < HOUGH_LINE_LENGTH) {

                lines.remove(i);
                i--;
            }
        }
        return lines;
    }

    private static void mergeLine(Line2d keep, Line2d remove) {

        //determinate the merge axis
        double hAngle = keep.calculateHorizontalAngle();
//        System.out.println("hAngle = "+hAngle);

        // Given the equation for a line as ax - by + c = 0
        float a = keep.end.getY() - keep.begin.getY();
        float b = keep.end.getX() - keep.begin.getX();
        float c = keep.end.getX() * keep.begin.getY() - keep.begin.getX() * keep.end.getY();

        float newX, newY;

        Point2d expandedPoint;
        if (hAngle < PI / 4) {//x-axis

            // begin.x < end.x
            sortByXAxis(keep);
            sortByXAxis(remove);

            if (remove.begin.getX() < keep.begin.getX()) {
                // expand left
                newX = remove.begin.getX();
                newY = (a * newX + c) / (b);
                if (newY < 0) {
                    newY = 0;
                    newX = (b * newY - c) / (a);
                } else if (newY > height) {
                    newY = height - 1;
                    newX = (b * newY - c) / (a);
                }
                keep.begin.setX(newX);
                keep.begin.setY(newY);
            } else if (remove.end.getX() > keep.end.getX()) {
                // expand right
                newX = remove.end.getX();
                newY = (a * newX + c) / (b);

                if (newY < 0) {
                    newY = 0;
                    newX = (b * newY - c) / (a);
                } else if (newY > height) {
                    newY = height - 1;
                    newX = (b * newY - c) / (a);
                }

                keep.end.setX(newX);
                keep.end.setY(newY);


            }
        } else { // y-axis
            // begin.x < end.x
            sortByYAxis(keep);
            sortByYAxis(remove);

            if (remove.begin.getY() < keep.begin.getY()) {
                // expand top
                newY = remove.begin.getY();
                newX = (b * newY - c) / (a);

                if (newX < 0) {
                    newX = 0;
                    newY = (a * newX + c) / (b);
                } else if (newX > width) {
                    newX = width - 1;
                    newY = (b * newY - c) / (a);
                }

                keep.begin.setX(newX);
                keep.begin.setY(newY);
            } else if (remove.end.getY() > keep.end.getY()) {
                // expand bottom
                newY = remove.end.getY();
                newX = (b * newY - c) / (a);

                if (newX < 0) {
                    newX = 0;
                    newY = (a * newX + c) / (b);
                } else if (newX > width) {
                    newX = width - 1;
                    newY = (b * newY - c) / (a);
                }

                keep.end.setX(newX);
                keep.end.setY(newY);
            }
        }
//        return keep;
    }

    private static void sortByYAxis(Line2d keep) {
        if (keep.begin.getY() > keep.end.getY()) {
            Point2d tmp = keep.begin;
            keep.setBeginPoint(keep.end);
            keep.setEndPoint(tmp);
        }
    }

    private static void sortByXAxis(Line2d keep) {
        if (keep.begin.getX() > keep.end.getX()) {
            Point2d tmp = keep.begin;
            keep.setBeginPoint(keep.end);
            keep.setEndPoint(tmp);
        }
    }





    private static FImage applyCannyDetector(FImage grey, File fin, File fout) throws IOException {
        CannyEdgeDetector canny = new CannyEdgeDetector(CANNY_LOW_THRESH, CANNY_HIGH_THRESH, CANNY_SIGMA);
//        FImage grey = frame.getBand(S_CHANNEL_ID);


//        grey.processInplace(new FGaussianConvolve(GAUSSIAN_BLUR_SIGMA));
//        ImageUtilities.write(grey, new File(fout.getAbsolutePath()+"/blur/"+fin.getName()));


        canny.processImage(grey);
        return grey;
    }


    static Random rnd = new Random();

    private static Float[] getRandomColor() {
        return RGBColour.RGB(rnd.nextInt(256), rnd.nextInt(256), rnd.nextInt(256));
    }

    static double distance(Point2d p1, Point2d p2) {

        float x1 = p1.getX(), x2 = p2.getX(), y1 = p1.getY(), y2 = p2.getY();

        return Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }

    private static void drawBound(MBFImage frame, Point2d center, List<Line2d> lines, Float[] baseColor, Float[] fitColor) {
        Float[] lineColor;
        for (Line2d line : lines) {
            if (lines.indexOf(line) % 2 == 0)
                lineColor = baseColor;
            else
                lineColor = fitColor;

            frame.drawLine(line, 2, lineColor);
            frame.drawPoint(line.begin, RGBColour.RED, 5);
            frame.drawPoint(line.end, RGBColour.YELLOW, 5);
        }
    }

    private static void drawLines(MBFImage frame, Point2d center, List<Line2d> lines, Float[] lineColor) {
        System.out.println(lines.size() + " lines");
        frame.drawText(lines.size() + "", (int) center.getX(),
                (int) center.getY(),
                HersheyFont.ROMAN_DUPLEX,
                40,
                RGBColour.RED);

        Float[] orgColor = lineColor;

        for (Line2d line : lines) {
            frame.drawLine(line, 2, lineColor);
            frame.drawPoint(line.begin, RGBColour.RED, 5);
            frame.drawPoint(line.end, RGBColour.YELLOW, 5);
            frame.drawText(Math.round(line.calculateHorizontalAngle() * 180 / PI) + "*", (int) line.calculateCentroid().getX(),
                    (int) line.calculateCentroid().getY(),
                    HersheyFont.ROMAN_DUPLEX,
                    20, RGBColour.BLUE);
        }
    }
}
