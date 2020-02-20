package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.FImage;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.colour.Transforms;
import org.openimaj.image.contour.Contour;
import org.openimaj.image.contour.SuzukiContourProcessor;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;
import org.openimaj.math.geometry.shape.Polygon;
import org.openimaj.math.geometry.shape.Rectangle;

import java.util.List;

import static uk.ac.soton.ecs.jsh2.App.*;
import static uk.ac.soton.ecs.jsh2.MyUtils.NEW_THRESHOLD;

public class OldAlgorithm {


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

}
