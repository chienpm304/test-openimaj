package uk.ac.soton.ecs.jsh2;

import Jama.Matrix;
import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.analysis.algorithm.HoughLines;
import org.openimaj.image.analysis.algorithm.histogram.HistogramAnalyser;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.colour.Transforms;
import org.openimaj.image.contour.Contour;
import org.openimaj.image.contour.SuzukiContourProcessor;
import org.openimaj.image.processing.edges.CannyEdgeDetector;
import org.openimaj.image.processing.resize.ResizeProcessor;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;
import org.openimaj.math.geometry.shape.Polygon;
import org.openimaj.math.geometry.shape.Rectangle;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * OpenIMAJ Hello world!
 *
 */
public class App {
    public static final float THRESHOLD_BIN_INV = 0.07133f;
    public static final float STANDARD_WIDTH = 600;
    public static final int S_CHANNEL_ID = 1;
    public static final String WINDOW_DIR = "D:/detect/input/all";
    private static final String LINUX_DIR = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/all";
    private static float scaleFactor = 1.0f;

    public static void main( String[] args ) throws IOException {
        File folder = new File(LINUX_DIR);
        if(folder.exists() && folder.isDirectory())
            for (final File file : folder.listFiles()) {
                if(file.isFile())
                    detectBox(file, folder);
            }
//        File file = new File("D:/detect/input/4/rsz_t2.jpg");
//        System.out.println("Processing "+file.getName());
//        MBFImage frame = ImageUtilities.readMBF(file);
//
//        Tetragram originTetra = detectBox(frame, null);
//        System.out.println("Bound:");
//        System.out.println(originTetra);
//
//        Tetragram destTetra = findDestinationRectangle(originTetra);
//        System.out.println("destimation rect:");
//        System.out.println(destTetra);
//
//        int newWidth = (int) destTetra.getBottomRight().getX() + 1;
//        int newHeight = (int) destTetra.getBottomRight().getY() + 1;
//        System.out.println("New size" + newWidth+"x"+newHeight);
//
//        Matrix transformMatrix = getTransformMatrix(originTetra, destTetra);
//
//        System.out.println(transformMatrix.getArray());
//
//        MBFImage dstImg = transformImage(frame, transformMatrix, newWidth, newHeight);
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
        for(int i =0; i < 4; i++){
            A[2*i][0] = src[i][0];
            A[2*i][1] = src[i][1];
            A[2*i][2] = 1;
            A[2*i][3] = 0;
            A[2*i][4] = 0;
            A[2*i][5] = 0;
            A[2*i][6] = -src[i][0]*dst[i][0];
            A[2*i][7] = -dst[i][0]*src[i][1];

            A[2*i + 1][0] = 0;
            A[2*i + 1][1] = 0;
            A[2*i + 1][2] = 0;
            A[2*i + 1][3] = src[i][0];
            A[2*i + 1][4] = src[i][1];
            A[2*i + 1][5] = 1;
            A[2*i + 1][6] = -src[i][0]*dst[i][1];
            A[2*i + 1][7] = -src[i][1]*dst[i][1];
        }

        // build A (8x8)
        Matrix matrixA = new Matrix(A);

        // build b (8x1)
        double []b = new double[8];
        for(int i = 0; i<4; i++){
            b[2*i] = dst[i][0];
            b[2*i+1] = dst[i][1];
        }

        Matrix matrixB = new Matrix(b, 8);

        Matrix trsf = matrixA.inverse().times(matrixB);

        // reshape to 3x3
        double[] flat = trsf.getColumnPackedCopy();
        Matrix trf3x3 = new Matrix(3,3);
        for(int i = 0; i<8; i++)
            trf3x3.set(i/3, i%3, flat[i]);
        trf3x3.set(2,2, 1);

        return trf3x3;
    }

    private static Tetragram findDestinationRectangle(Tetragram bound) {
        Point2d tl = bound.getTopLeft();
        Point2d tr = bound.getTopRight();
        Point2d br = bound.getBottomRight();
        Point2d bl = bound.getBottomLeft();

        int newWidth = (int) Math.max(distance(br, bl), distance(tr, tl));
        int newHeight = (int) Math.max(distance(tr, br), distance(tl, bl));

        float[][] rect = new float[][]{
                {0, 0},
                {newWidth-1, 0},
                {newWidth-1, newHeight-1},
                {0, newHeight-1}
        };

        return new Tetragram(rect);
    }

    private static Tetragram detectBox(File fin, File folder) throws IOException {

        MBFImage frame = ImageUtilities.readMBF(fin);
        scaleFactor = STANDARD_WIDTH /(float)frame.getWidth();
        scaleFactor = scaleFactor>1?1.0f:scaleFactor;
        System.out.println("processing: " + fin.getName() +"...");
//        System.out.println("scale: "+scaleFactor);
//        frame.processInplace(new ResizeProcessor(scaleFactor));

        FImage grey = applyCustomPreproccessing(frame);
//        ImageUtilities.write(grey, new File(folder.getAbsolutePath()+"/edged/"+fin.getName()));

        List<Point2d> contour = getContour(grey);
        System.out.println("contour: "+contour.size());
        frame.drawPoints(contour, RGBColour.GREEN, 4);

        HoughLinesP ht = new HoughLinesP(contour, grey.width, grey.height, 1,
                Math.PI/360, 50, 230, 250, 50);

        List<Line2d> lines = ht.getLines();

        System.out.println(lines.size() + " lines");

        for(Line2d line: lines){
            System.out.println(line.toString());
            frame.drawLine(line, 2, RGBColour.RED);
        }

        ImageUtilities.write(frame, new File(folder.getAbsolutePath()+"/out/"+fin.getName()));
        return null;
    }

    private static List<Point2d> getContour(FImage grey) {
        Contour contour = SuzukiContourProcessor.findContours(grey);
        Rectangle max = new Rectangle(0,0, 1, 1);
        Polygon fit = contour.clone();

        FImage tmp = new FImage(grey.width, grey.height);
        for(int i = 0; i < contour.children.size(); i++){
            Contour cnt = contour.children.get(i);
            Rectangle box = cnt.calculateRegularBoundingBox();
            double area = box.calculateArea();
            if(area > max.calculateArea() && area>1 && area < contour.calculateRegularBoundingBox().calculateArea()) {
                max = box.clone();
                fit = cnt.clone();
            }
        }
        return fit.points;
    }

    private static Tetragram findBoundingBoxByHough(MBFImage frame, FImage grey) {
        TestHough houghLines = new TestHough();
        houghLines.analyseImage(grey);

        Iterator iterator = houghLines.iterator();
        while(iterator.hasNext()){
            Line2d line = houghLines.next();
            System.out.println(line.toString());
            frame.drawLine(houghLines.next(), 2, RGBColour.RED);
        }

//        List<Line2d> lines = houghLines.getBestLines( 4);
//        for(Line2d l: lines){
//            frame.drawLine(l, 2, RGBColour.RED);
//        }
        return null;
    }

    private static Tetragram findBoundingBoxByContour(MBFImage frame, FImage grey) {
        Contour contour = SuzukiContourProcessor.findContours(grey);
        Rectangle max = new Rectangle(0,0, 1, 1);
        Polygon fit = contour.clone();


        FImage tmp = new FImage(grey.width, grey.height);
        for(int i = 0; i < contour.children.size(); i++){

            Contour cnt = contour.children.get(i);

            Rectangle box = cnt.calculateRegularBoundingBox();
            double area = box.calculateArea();
            if(area > max.calculateArea() && area>1 && area < contour.calculateRegularBoundingBox().calculateArea()) {
                max = box.clone();
                fit = cnt.clone();
//                tmp.drawPoints(cnt.points, 1f, 1);
            }
        }
//        ImageUtilities.write(tmp, new File(folder.getAbsolutePath()+"/contour/"+fin.getName()));

        System.out.println(max.calculateArea()+"| "+max.toString());

        Tetragram bound = findBounding(fit.points, scaleFactor);

        frame.drawPoint(fit.asPolygon().calculateCentroid(), RGBColour.RED, 20);


//        frame.drawShape(max.scale(1.0f/scaleFactor);, 8, RGBColour.BLUE);

//        frame.drawLine(bound.getTopLeft(), bound.getTopRight(), 10, RGBColour.BLUE);
//        frame.drawLine(bound.getTopRight(), bound.getBottomRight(), 10, RGBColour.BLUE);
//        frame.drawLine(bound.getBottomRight(), bound.getBottomLeft(), 10, RGBColour.BLUE);
//        frame.drawLine(bound.getBottomLeft(), bound.getTopLeft(), 10, RGBColour.BLUE);

        frame.drawPoints(bound.toList(), RGBColour.BLUE, 10);

        frame.drawPoints(fit.points, RGBColour.GREEN, 10);

        System.out.println(bound);
        return bound;
    }

    private static FImage applyCannyDetector(MBFImage frame) {
        CannyEdgeDetector canny = new CannyEdgeDetector(0.05f, 0.15f, 3.0f);
        FImage grey = frame.flatten();

//        scaleFactor = STANDARD_WIDTH /(float)frame.getWidth();
//        scaleFactor = scaleFactor>1?1.0f:scaleFactor;
//        System.out.println("scale: "+scaleFactor);
//        grey.processInplace(new ResizeProcessor(scaleFactor));

        canny.processImage(grey);
        return grey;
    }

    private static FImage applyCustomPreproccessing(MBFImage frame) {
        //convert to hsv
        MBFImage hsv = Transforms.RGB_TO_HSV(frame);

        // get S channel and RESIZE
//        scaleFactor = STANDARD_WIDTH /(float)frame.getWidth();

//        scaleFactor = scaleFactor>1?1.0f:scaleFactor;
//        System.out.println("scale: "+scaleFactor);
        FImage s = hsv.getBand(S_CHANNEL_ID);//.processInplace(new ResizeProcessor(scaleFactor));

//        s.analyseWith(new HistogramAnalyser(64));

        return s.threshold(THRESHOLD_BIN_INV).inverse();
    }

    private static List<Point2d> getOriginScale(List<Point2d> points, float scaleFactor) {
        for (Point2d p: points) {
            p.setX(p.getX() / scaleFactor);
            p.setY(p.getY()/scaleFactor);
        }
        return points;
    }

    private static Tetragram findBounding(List<Point2d> points, float scale_factor) {
        float x_min = Float.MAX_VALUE;
        float x_max = -1;
        float y_min = Float.MAX_VALUE;
        float y_max = -1;
        float x,y;
        for(Point2d p: points){
            x = p.getX();
            y = p.getY();
            if(x > x_max) x_max = x;
            if(x < x_min) x_min = x;
            if(y < y_min) y_min = y;
            if(y > y_max) y_max = y;

        }

        Point2dImpl
                tl = new Point2dImpl(x_min, y_max),
                tr =new Point2dImpl(x_max, y_max),
                bl = new Point2dImpl(x_min,y_min),
                br = new Point2dImpl(x_max, y_min);

        Point2dImpl r_tl = br.clone(), r_tr = bl.clone(), r_bl = tr.clone(), r_br = tl.clone();


        float width = x_max - x_min, height = y_max - y_min;

        for(Point2d p: points){
            x = p.getX();
            y = p.getY();
            //top left
            if(x < width/2 && y>height/2){
                if(distance(p, tl) < distance(r_tl, tl)){
                    r_tl.x = p.getX();
                    r_tl.y = p.getY();
                }
            }
            else if(x > width/2 && y > height/2) {// top right
                if(distance(p, tr) < distance(r_tr, tr)){
                    r_tr.x = p.getX();
                    r_tr.y = p.getY();
                }
            }
            else if(x < width/2 && y < height/2) {// bottom left
                if(distance(p, bl) < distance(r_bl, bl)){
                    r_bl.x = p.getX();
                    r_bl.y = p.getY();
                }
            }
            else {//bottom right
                if(distance(p, br) < distance(r_br, br)){
                    r_br.x = p.getX();
                    r_br.y = p.getY();
                }
            }

        }

        // small add]just
        if(distance(r_bl, bl) > width/3) r_bl = bl;
        if(distance(r_tl, tl) > width/3) r_tl = tl;
        if(distance(r_br, br) > width/3) r_br = br;
        if(distance(r_tr, tr) > width/3) r_tr = tr;


        // return original scale
        r_tl.x = (int)(r_tl.x / scale_factor);
        r_tl.y = (int)(r_tl.y / scale_factor);

        r_tr.x = (int)(r_tr.x / scale_factor);
        r_tr.y = (int)(r_tr.y / scale_factor);

        r_bl.x = (int)(r_bl.x / scale_factor);
        r_bl.y = (int)(r_bl.y / scale_factor);

        r_br.x = (int)(r_br.x / scale_factor);
        r_br.y = (int)(r_br.y / scale_factor);


        Tetragram bound = new Tetragram(r_tl, r_tr, r_br, r_bl);
        return bound;
    }

    static Random rnd = new Random();

    private static Float[] getRandomColor() {
        return RGBColour.RGB(rnd.nextInt(256), rnd.nextInt(256), rnd.nextInt(256));
    }


    static double distance(Point2d p1, Point2d p2){

        float x1 = p1.getX(), x2 = p2.getX(), y1 = p1.getY(), y2 = p2.getY();

        return Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }
}
