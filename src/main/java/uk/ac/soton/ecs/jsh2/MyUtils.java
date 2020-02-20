package uk.ac.soton.ecs.jsh2;

import Jama.Matrix;
import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.processing.resize.ResizeProcessor;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;
import org.openimaj.math.geometry.shape.Polygon;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class MyUtils {
    static int width = 720;
    static int height = 1080;
    public static final int THRESHOLD_STEP = 50;
    public static final float NEW_THRESHOLD = 0.1f;


    private static List<Line2d> selectRawLines(Point2d center, List<Line2d> lines) {
        Line2d top = null, right = null, bottom = null, left = null;
        top = findTop(width, height, center, lines);
        lines.remove(top);

        right = findRight(width, height, center, lines);
        lines.remove(right);

        bottom = findBottom(width, height, center, lines);
        lines.remove(bottom);

        left = findLeft(width, height, center, lines);
        lines.remove(left);

        List<Line2d> l = new ArrayList<>();
        l.add(top);
        l.add(right);
        l.add(bottom);
        l.add(left);
        return l;
    }

    private static Tetragram getBounding(Point2d center, List<Line2d> lines) {

        Line2d top = null, right = null, bottom = null, left = null;
        top = findTop(width, height, center, lines);
        lines.remove(top);

        right = findRight(width, height, center, lines);
        lines.remove(right);

        bottom = findBottom(width, height, center, lines);
        lines.remove(bottom);

        left = findLeft(width, height, center, lines);
        lines.remove(left);


//        System.out.println("Top: "+top);
//        System.out.println("Bottom: "+bottom);
//        System.out.println("Right: "+right);
//        System.out.println("Left: "+left);


        Point2d topLeft = findLinesIntersection(top, left);
        Point2d topRight = findLinesIntersection(top, right);
        Point2d bottomRight = findLinesIntersection(bottom, right);
        Point2d bottomLeft = findLinesIntersection(bottom, left);

        return new Tetragram(topLeft, topRight, bottomRight, bottomLeft);
    }

    public static Point2d findLinesIntersection(Line2d line1, Line2d line2) {
        Point2d p = new Point2dImpl(0, 0);
        Line2d.IntersectionResult result = MyHelper.findIntersection(line1, line2);
        if (result.type == Line2d.IntersectionType.INTERSECTING) {
            p = result.intersectionPoint;
            if (p.getX() < 0) p.setX(0);
            else if (p.getX() >= width) p.setX(width - 1);
            if (p.getY() < 0) p.setY(0);
            else if (p.getY() >= height) p.setY(height - 1);
        } else {
            System.out.println("Not intersecting: " + line1 + " with " + line2);
        }
        return p;
    }

    private static Line2d findTop(int width, int height, Point2d center, List<Line2d> lines) {
        double max = -1;
        int threshold = 0;
        Line2d top = null;
        while (top == null && threshold < height / 3) {
            max = -1;
            for (Line2d l : lines) {
                if (l.begin.getY() < center.getY() + threshold
                        && l.end.getY() < center.getY() + threshold
//                        && Math.abs(l.begin.getY() - l.end.getY()) < height/1.5f) {
                ) {
                    if (l.calculateLength() > max) {
                        max = l.calculateLength();
                        top = l;
                    }
                }
            }
            threshold += THRESHOLD_STEP;
        }
        if (top == null)
            top = new Line2d(0, 0, width - 1, 0);
        return top;
    }

    private static Line2d findRight(int width, int height, Point2d center, List<Line2d> lines) {
        double max;
        int threshold = 0;
        Line2d right = null;
        while (right == null && threshold < width / 3) {
            max = -1;
            for (Line2d l : lines) {
                if (l.begin.getX() > center.getX() - threshold
                        && l.end.getX() > center.getX() - threshold
//                        && Math.abs(l.begin.getX() - l.end.getX()) < width/1.5f) {
                ) {
                    if (l.calculateLength() > max) {
                        max = l.calculateLength();
                        right = l;
                    }
                }
            }
            threshold += THRESHOLD_STEP;
        }
        if (right == null)
            right = new Line2d(width - 1, 0, width - 1, height - 1);
        return right;
    }

    private static Line2d findBottom(int width, int height, Point2d center, List<Line2d> lines) {
        double max;
        int threshold = 0;
        Line2d bottom = null;
        while (bottom == null && threshold < height / 3) {
            max = -1;
            for (Line2d l : lines) {
                if (l.begin.getY() > center.getY() - threshold
                        && l.end.getY() > center.getY() - threshold
//                        && Math.abs(l.begin.getY()-l.end.getY())<height/1.5f) {
                ) {
                    if (l.calculateLength() > max) {
                        max = l.calculateLength();
                        bottom = l;
                    }
                }
            }
            threshold += THRESHOLD_STEP;
        }
        if (bottom == null)
            bottom = new Line2d(0, height - 1, width - 1, height - 1);
        return bottom;
    }

    private static Line2d findLeft(int width, int height, Point2d center, List<Line2d> lines) {
        double max;
        int threshold = 0;
        Line2d left = null;
        while (left == null && threshold < width / 3) {
            max = -1;
            for (Line2d l : lines) {
                if (l.begin.getX() < center.getX() + threshold
                        && l.end.getX() < center.getX() + threshold
//                        && Math.abs(l.begin.getX() - l.end.getX()) < width/1.5f) {
                ) {
                    if (l.calculateLength() > max) {
                        max = l.calculateLength();
                        left = l;
                    }
                }
            }
            threshold += THRESHOLD_STEP;
        }
        if (left == null)
            left = new Line2d(0, 0, 0, height - 1);
        return left;
    }


//    private static Tetragram detectWithCustomProcessing(File fin, File folder) throws IOException {
//
//        MBFImage frame = ImageUtilities.readMBF(fin);
//        scaleFactor = STANDARD_WIDTH /(float)frame.getWidth();
//        scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;
//
//        frame = frame.process(new ResizeProcessor(scaleFactor));
//        width = frame.getWidth();
//        height = frame.getHeight();
//
//        System.out.println("processing: " + fin.getName() +"...");
//
//        FImage grey = applyCustomPreproccessing(frame);
//        ImageUtilities.write(grey, new File(folder.getAbsolutePath()+"/edged/"+fin.getName()));
//
//        List<Point2d> contour = getContour(grey);
//
//        System.out.println("contour: "+contour.size());
//        frame.drawPoints(contour, RGBColour.GREEN, 4);
//        ImageUtilities.write(frame, new File(folder.getAbsolutePath()+"/contour/"+fin.getName()));
//
//        Point2d center = new Polygon(contour).calculateCentroid();
//
//        frame.drawPoint(center, RGBColour.RED, 20);
//        List<Line2d> lines = getLinesUsingHoughTransformP(contour);
//
//
//        lines = removeSimilarAndNoiseLines(lines);
//
//        drawLines(frame, center, lines, RGBColour.BLUE);
//
//
//        Tetragram bounding = getBounding(center, lines);//.scale(1.0f/scaleFactor);
////        System.out.println("bounding scaled: "+bounding.toString());
////        drawBounds(frame, bounding);
//
////        bounding = bounding.scale(1.0f/scaleFactor);
//        System.out.println("bounding: "+bounding.toString());
////        drawBounds(frame, bounding);
//
//        ImageUtilities.write(frame, new File(folder.getAbsolutePath()+"/out/"+fin.getName()));
//        return bounding;
//    }

    private static void drawBounds(MBFImage frame, Tetragram bounding, Float[] lineColor) {
        frame.drawLine(bounding.getTopLeft(), bounding.getTopRight(), 3, lineColor);
        frame.drawLine(bounding.getTopLeft(), bounding.getBottomLeft(), 3, lineColor);
        frame.drawLine(bounding.getTopRight(), bounding.getBottomRight(), 3, lineColor);
        frame.drawLine(bounding.getBottomLeft(), bounding.getBottomRight(), 3, lineColor);
        frame.drawPoints(bounding.toList(), RGBColour.RED, 8);
    }

    private static Tetragram findDestinationRectangle(Tetragram bound) {
        Point2d tl = bound.getTopLeft();
        Point2d tr = bound.getTopRight();
        Point2d br = bound.getBottomRight();
        Point2d bl = bound.getBottomLeft();

        int newWidth = (int) Math.max(App.distance(br, bl), App.distance(tr, tl));
        int newHeight = (int) Math.max(App.distance(tr, br), App.distance(tl, bl));

        float[][] rect = new float[][]{
                {0, 0},
                {newWidth - 1, 0},
                {newWidth - 1, newHeight - 1},
                {0, newHeight - 1}
        };

        return new Tetragram(rect);
    }


    private static MBFImage transformImage(MBFImage frame, Matrix transformMatrix, int newWidth, int newHeight) {
        MBFImage img = new MBFImage(newWidth, newHeight);
//        img.processInplace(new ;
        return null;
    }

    private static Matrix getTransformMatrix(Tetragram originTetra, Tetragram destTetra) {
        float[][] src = originTetra.toArray();
        float[][] dst = destTetra.toArray();

        double A[][] = new double[8][8];
        // A x h = b -> find h = (A)^-1 x b
        for (int i = 0; i < 4; i++) {
            A[2 * i][0] = src[i][0];
            A[2 * i][1] = src[i][1];
            A[2 * i][2] = 1;
            A[2 * i][3] = 0;
            A[2 * i][4] = 0;
            A[2 * i][5] = 0;
            A[2 * i][6] = -src[i][0] * dst[i][0];
            A[2 * i][7] = -dst[i][0] * src[i][1];

            A[2 * i + 1][0] = 0;
            A[2 * i + 1][1] = 0;
            A[2 * i + 1][2] = 0;
            A[2 * i + 1][3] = src[i][0];
            A[2 * i + 1][4] = src[i][1];
            A[2 * i + 1][5] = 1;
            A[2 * i + 1][6] = -src[i][0] * dst[i][1];
            A[2 * i + 1][7] = -src[i][1] * dst[i][1];
        }

        // build A (8x8)
        Matrix matrixA = new Matrix(A);

        // build b (8x1)
        double[] b = new double[8];
        for (int i = 0; i < 4; i++) {
            b[2 * i] = dst[i][0];
            b[2 * i + 1] = dst[i][1];
        }

        Matrix matrixB = new Matrix(b, 8);

        Matrix trsf = matrixA.inverse().times(matrixB);

        // reshape to 3x3
        double[] flat = trsf.getColumnPackedCopy();
        Matrix trf3x3 = new Matrix(3, 3);
        for (int i = 0; i < 8; i++)
            trf3x3.set(i / 3, i % 3, flat[i]);
        trf3x3.set(2, 2, 1);

        return trf3x3;
    }

}