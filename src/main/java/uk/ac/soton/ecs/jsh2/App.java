package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.analysis.algorithm.histogram.HistogramAnalyser;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.colour.Transforms;
import org.openimaj.image.contour.Contour;
import org.openimaj.image.contour.SuzukiContourProcessor;
import org.openimaj.image.processing.convolution.FGaussianConvolve;
import org.openimaj.image.processing.edges.CannyEdgeDetector;
import org.openimaj.image.processing.resize.ResizeProcessor;
import org.openimaj.image.typography.hershey.HersheyFont;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;
import org.openimaj.math.geometry.shape.Polygon;
import org.openimaj.math.geometry.shape.Rectangle;

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
    public static final float STANDARD_WIDTH = 720;
    public static final int S_CHANNEL_ID = 1;

    public static final String WINDOW_DIR = "D:/detect/input/AZdoc/in";
    public static final String WINDOW_OUT_DIR = "D:/detect/input/AZdoc/java";

    private static final String LINUX_DIR = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/AZdoc/in";
    private static final String LINUX_DIR_out = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/AZdoc/java";
    public static final int THRESHOLD_STEP = 50;
    public static final float NEW_THRESHOLD = 0.1f;



    private static float scaleFactor = 1.0f;
    private static int height;
    private static int width;

    public static final float GAUSSIAN_BLUR_SIGMA = 5f;

    public static final float CANNY_LOW_THRESH = 0.02f;
    public static final float CANNY_HIGH_THRESH = 0.03f;
    public static final float CANNY_SIGMA = 5f;

    public static final int HOUGHLINE_THRESHOLD = 80;
    public static final int HOUGH_LINE_LENGTH = 200;

    private static final int LINE_GAP_REMOVAL = 25;
    private static final int BOUNDING_GAP_REMOVAL = 0;

    public static final int HISTOGRAM_NBINS = 64;


    public static void testIntersect(){
        Line2d.IntersectionResult res =
                MyHelper.findIntersection(
                        new Line2d(6, 5, 8, 5),
                        new Line2d(1, 5, 3, 5));
        System.out.println(res.type);
        System.out.println(res.intersectionPoint);
    }

    public static void main( String[] args ) throws IOException {
        testDetectBox();
//        testIntersect();
    }

    private static void testDetectBox() throws IOException {
        File folder = new File(LINUX_DIR);
        File out = new File(LINUX_DIR_out);
        if(folder.exists() && folder.isDirectory())
            for (final File file : folder.listFiles()) {
                if(file.isFile())
                    detectWithCanny(file, out);
            }
    }

    private static Tetragram detectWithCanny(File fin, File folder) throws IOException{
        MBFImage frame = ImageUtilities.readMBF(fin);
        scaleFactor = STANDARD_WIDTH /(float)frame.getWidth();
        scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;

        frame = frame.process(new ResizeProcessor(scaleFactor));
        width = frame.getWidth();
        height = frame.getHeight();

        System.out.println("processing: " + fin.getName() +"...");

        FImage edges = applyCannyDetector(frame);
        ImageUtilities.write(edges, new File(folder.getAbsolutePath()+"/edged/"+fin.getName()));

        List<Line2d> lines = getLinesUsingHoughTransformP(edges);
        drawLines(frame, new Point2dImpl(width/2, height/2), lines, RGBColour.BLUE);

        lines = removeSimilarAndNoiseLines(lines);

        Point2dImpl center =  new Point2dImpl(width/2, height/2);

//        drawLines(frame,center, lines, RGBColour.BLUE);

        drawLines(frame, center, selectRawLines(center, lines), RGBColour.GREEN);

//        Tetragram bounding = getBounding(new Point2dImpl(width/2, height/2), lines);
//        System.out.println("bounding scaled: "+bounding.toString());
//        drawBounds(frame, bounding, RGBColour.GREEN);

//        bounding = bounding.scale(1.0f/scaleFactor);
//        System.out.println("bounding: "+bounding.toString());
//        drawBounds(frame, bounding);

        ImageUtilities.write(frame, new File(folder.getAbsolutePath()+"/out/"+fin.getName()));
        return null;
    }

    private static List<Line2d> getLinesUsingHoughTransformP(FImage image) {
        List<Point2d> points = new ArrayList<>();
        for(int i = 0; i < image.width; i++)
            for(int j = 0; j < image.height; j++){
                if(image.getPixel(i, j) > 0.5f)
                    points.add(new Point2dImpl(i, j));
            }
        return getLinesUsingHoughTransformP(points);
    }

    private static Tetragram detectBox(File fin, File folder) throws IOException {

        MBFImage frame = ImageUtilities.readMBF(fin);
        scaleFactor = STANDARD_WIDTH /(float)frame.getWidth();
        scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;

        frame = frame.process(new ResizeProcessor(scaleFactor));
        width = frame.getWidth();
        height = frame.getHeight();

        System.out.println("processing: " + fin.getName() +"...");

        FImage grey = applyCustomPreproccessing(frame);
        ImageUtilities.write(grey, new File(folder.getAbsolutePath()+"/edged/"+fin.getName()));

        List<Point2d> contour = getContour(grey);

        System.out.println("contour: "+contour.size());
        frame.drawPoints(contour, RGBColour.GREEN, 4);
        ImageUtilities.write(frame, new File(folder.getAbsolutePath()+"/contour/"+fin.getName()));

        Point2d center = new Polygon(contour).calculateCentroid();

        frame.drawPoint(center, RGBColour.RED, 20);
        List<Line2d> lines = getLinesUsingHoughTransformP(contour);


        lines = removeSimilarAndNoiseLines(lines);

        drawLines(frame, center, lines, RGBColour.BLUE);


        Tetragram bounding = getBounding(center, lines);//.scale(1.0f/scaleFactor);
//        System.out.println("bounding scaled: "+bounding.toString());
//        drawBounds(frame, bounding);

//        bounding = bounding.scale(1.0f/scaleFactor);
        System.out.println("bounding: "+bounding.toString());
//        drawBounds(frame, bounding);

        ImageUtilities.write(frame, new File(folder.getAbsolutePath()+"/out/"+fin.getName()));
        return bounding;
    }

    private static void drawBounds(MBFImage frame, Tetragram bounding, Float[] lineColor) {
        frame.drawLine(bounding.getTopLeft(), bounding.getTopRight(), 3, lineColor);
        frame.drawLine(bounding.getTopLeft(), bounding.getBottomLeft(), 3, lineColor);
        frame.drawLine(bounding.getTopRight(), bounding.getBottomRight(), 3, lineColor);
        frame.drawLine(bounding.getBottomLeft(), bounding.getBottomRight(), 3, lineColor);
        frame.drawPoints(bounding.toList(), RGBColour.RED, 8);
    }


    private static List<Line2d> removeSimilarAndNoiseLines(List<Line2d> lines) {
        for(int i = 0; i < lines.size(); i++){
            Line2d l = lines.get(i);

            if(l.begin.getX() < BOUNDING_GAP_REMOVAL && l.end.getX() < BOUNDING_GAP_REMOVAL //LEFT

                || l.begin.getY() < BOUNDING_GAP_REMOVAL && l.end.getY() < BOUNDING_GAP_REMOVAL //TOP

                || l.begin.getX() > width - BOUNDING_GAP_REMOVAL && l.end.getX() > width - BOUNDING_GAP_REMOVAL // RIGHT

                || l.begin.getY() > height - BOUNDING_GAP_REMOVAL && l.end.getY() > height - BOUNDING_GAP_REMOVAL) {

                   lines.remove(i);
                   i--;
            }
        }

        Line2d l1, l2;
        for(int i = 0; i < lines.size()-1; i++){
            l1 = lines.get(i);
            for(int j = i+1; j < lines.size(); j++){
                l2 = lines.get(j);
                if(l1 != l2
                        && l2.distanceToLine(l1.begin) < LINE_GAP_REMOVAL
                        && l2.distanceToLine(l1.end) < LINE_GAP_REMOVAL
                        && l1.distanceToLine(l2.begin) < LINE_GAP_REMOVAL
                        && l1.distanceToLine(l2.end) < LINE_GAP_REMOVAL
                ){
                    if(l1.calculateLength() < l2.calculateLength()) {
                        Collections.swap(lines, i, j);
                        lines.remove(j);
                        j = i;
                        l1 = lines.get(i);
                    }else{
                        lines.remove(j);
                        j--;
                    }
                }
            }
        }
        return lines;
    }

    private static List<Line2d> selectRawLines(Point2d center, List<Line2d> lines){
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

    private static Point2d findLinesIntersection(Line2d line1, Line2d line2) {
        Point2d p = new Point2dImpl(0,0);
        Line2d.IntersectionResult result = MyHelper.findIntersection(line1, line2);
        if(result.type == Line2d.IntersectionType.INTERSECTING){
            p = result.intersectionPoint;
            if(p.getX() < 0) p.setX(0);
            else if(p.getX() >= width) p.setX(width-1);
            if(p.getY() < 0) p.setY(0);
            else if(p.getY() >= height) p.setY(height-1);
        }else{
            System.out.println("Not intersecting: "+line1+" with "+line2);
        }
        return p;
    }

    private static Line2d findTop(int width, int height, Point2d center, List<Line2d> lines) {
        double max = -1;
        int threshold = 0;
        Line2d top = null;
        while (top == null && threshold < height/3){
            max = -1;
            for (Line2d l : lines) {
                if (l.begin.getY() < center.getY() + threshold
                        && l.end.getY() < center.getY() + threshold
//                        && Math.abs(l.begin.getY() - l.end.getY()) < height/1.5f) {
                ){
                    if (l.calculateLength() > max) {
                        max = l.calculateLength();
                        top = l;
                    }
                }
            }
            threshold += THRESHOLD_STEP;
        }
        if(top == null)
            top = new Line2d(0,0, width-1, 0);
        return top;
    }

    private static Line2d findRight(int width, int height, Point2d center, List<Line2d> lines) {
        double max;
        int threshold = 0;
        Line2d right = null;
        while (right == null && threshold < width/3){
            max = -1;
            for (Line2d l : lines) {
                if (l.begin.getX() > center.getX() - threshold
                        && l.end.getX() > center.getX()-threshold
//                        && Math.abs(l.begin.getX() - l.end.getX()) < width/1.5f) {
                ){
                    if (l.calculateLength() > max) {
                        max = l.calculateLength();
                        right = l;
                    }
                }
            }
            threshold += THRESHOLD_STEP;
        }
        if(right == null)
            right = new Line2d(width-1,0, width-1, height-1);
        return right;
    }

    private static Line2d findBottom(int width, int height, Point2d center, List<Line2d> lines) {
        double max;
        int threshold = 0;
        Line2d bottom = null;
        while (bottom == null && threshold < height/3){
            max = -1;
            for (Line2d l : lines) {
                if (l.begin.getY() > center.getY() - threshold
                        && l.end.getY() > center.getY() - threshold
//                        && Math.abs(l.begin.getY()-l.end.getY())<height/1.5f) {
                ){
                    if (l.calculateLength() > max) {
                        max = l.calculateLength();
                        bottom = l;
                    }
                }
            }
            threshold += THRESHOLD_STEP;
        }
        if(bottom == null)
            bottom = new Line2d(0,height-1, width-1, height-1);
        return bottom;
    }

    private static Line2d findLeft(int width, int height, Point2d center, List<Line2d> lines) {
        double max;
        int threshold = 0;
        Line2d left = null;
        while (left == null && threshold < width/3){
            max = -1;
            for (Line2d l : lines) {
                if (l.begin.getX() < center.getX() + threshold
                        && l.end.getX() < center.getX() + threshold
//                        && Math.abs(l.begin.getX() - l.end.getX()) < width/1.5f) {
                ){
                    if (l.calculateLength() > max) {
                        max = l.calculateLength();
                        left = l;
                    }
                }
            }
            threshold += THRESHOLD_STEP;
        }
        if(left == null)
            left = new Line2d(0, 0, 0, height-1);
        return left;
    }


    private static void drawLines(MBFImage frame, Point2d center, List<Line2d> lines, Float[] lineColor) {
        System.out.println(lines.size() + " lines");
        frame.drawText(lines.size()+"", (int)center.getX(),
                (int)center.getY(),
                HersheyFont.ROMAN_DUPLEX,
                40,
                RGBColour.RED);

        for(Line2d line: lines){
//            line.scale(1.0f/scaleFactor);
//            System.out.println(line.toString());
            frame.drawLine(line, 2, lineColor);
            frame.drawPoint(line.begin, RGBColour.RED, 5);
            frame.drawPoint(line.end, RGBColour.YELLOW, 5);
            frame.drawText(lines.indexOf(line)+1+"", (int)line.calculateCentroid().getX(),
                    (int)line.calculateCentroid().getY()+80,
                    HersheyFont.ROMAN_DUPLEX,
                    20,
                    lineColor);
        }
    }

    private static List<Line2d> getLinesUsingHoughTransformP(List<Point2d> contour) {
        HoughLinesP ht = new HoughLinesP(contour, width, height,  HOUGHLINE_THRESHOLD, HOUGH_LINE_LENGTH, 200);
        return ht.getLines();
    }

    private static List<Point2d> getContour(FImage grey) {
        Contour contour = SuzukiContourProcessor.findContours(grey);
        Rectangle max = new Rectangle(0,0, 1, 1);
        Polygon fit = contour.clone();


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
        CannyEdgeDetector canny = new CannyEdgeDetector(CANNY_LOW_THRESH, CANNY_HIGH_THRESH, CANNY_SIGMA);
        FImage grey = frame.flatten();


        grey.processInplace(new FGaussianConvolve(GAUSSIAN_BLUR_SIGMA));

//        HistogramAnalyser analyser = new HistogramAnalyser(HISTOGRAM_NBINS);
//        analyser.analyseImage(grey);

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
//
//        scaleFactor = scaleFactor>1?1.0f:scaleFactor;
//        System.out.println("scale: "+scaleFactor);
        FImage s = hsv.getBand(S_CHANNEL_ID);//.processInplace(new ResizeProcessor(scaleFactor));
//        DisplayUtilities.display(s);
//        s.analyseWith(new HistogramAnalyser(64));

        return s.threshold(NEW_THRESHOLD).inverse();
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
