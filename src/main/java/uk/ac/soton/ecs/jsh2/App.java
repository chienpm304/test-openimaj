package uk.ac.soton.ecs.jsh2;

import com.google.gson.internal.$Gson$Preconditions;
import org.apache.lucene.search.grouping.CollectedSearchGroup;
import org.openimaj.image.DisplayUtilities;
import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.colour.Transforms;
import org.openimaj.image.contour.Contour;
import org.openimaj.image.contour.SuzukiContourProcessor;
import org.openimaj.image.processing.resize.ResizeProcessor;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;
import org.openimaj.math.geometry.shape.Polygon;
import org.openimaj.math.geometry.shape.Rectangle;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * OpenIMAJ Hello world!
 *
 */
public class App {
    String root = "D:/detect";

    public static void main( String[] args ) throws IOException {
        File folder = new File("D:/detect/input/5");
        if(folder.exists() && folder.isDirectory())
            for (final File file : folder.listFiles()) {
                if(file.isFile())
                    detectBox(file, folder);
            }

//        detectBox(new File("D:/detect/input/4/3.jpg"), null);

    }

    private static void detectBox(File fin, File folder) throws IOException {
        System.out.println("Processing "+fin.getName());
        MBFImage frame = ImageUtilities.readMBF(fin);

        //convert to hsv
        MBFImage hsv = Transforms.RGB_TO_HSV(frame);




        // get S channel and RESIZE
        float scaleFactor = 800f/(float)frame.getWidth();
        FImage s = hsv.getBand(1).processInplace(new ResizeProcessor(scaleFactor));

        FImage grey = s.threshold(0.07133f).inverse();
//        display(grey);

        Contour contour = SuzukiContourProcessor.findContours(grey);
//        System.out.println(contour);


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

        System.out.println(max.calculateArea()+"| "+max.toString());

        List<Point2d> bound = findBounding(fit.points, scaleFactor);

        frame.drawPoints(fit.points, RGBColour.GREEN, 4);

//        frame.drawShape(max, 8, RGBColour.BLUE);

        frame.drawLine(bound.get(0), bound.get(1), 10, RGBColour.BLUE);
        frame.drawLine(bound.get(1), bound.get(2), 10, RGBColour.BLUE);
        frame.drawLine(bound.get(2), bound.get(3), 10, RGBColour.BLUE);
        frame.drawLine(bound.get(3), bound.get(0), 10, RGBColour.BLUE);

        frame.drawPoints(bound, RGBColour.RED, 10);


        System.out.println(bound);

//        DisplayUtilities.display(frame);
        ImageUtilities.write(frame, new File(folder.getAbsolutePath()+"/out/"+fin.getName()));
    }

    private static List<Point2d> findBounding(List<Point2d> points, float scale_factor) {
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
        r_tl.x /= scale_factor;
        r_tl.y /= scale_factor;
        r_tr.x /= scale_factor;
        r_tr.y /= scale_factor;
        r_bl.x /= scale_factor;
        r_bl.y /= scale_factor;
        r_br.x /= scale_factor;
        r_br.y /= scale_factor;


        List<Point2d> res = new ArrayList<>();
        res.add(r_tl);
        res.add(r_tr);
        res.add(r_br);
        res.add(r_bl);
        return res;
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
