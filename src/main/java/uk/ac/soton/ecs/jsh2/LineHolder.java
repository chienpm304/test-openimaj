package uk.ac.soton.ecs.jsh2;

import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;
import org.openimaj.math.geometry.shape.Polygon;
import uk.ac.soton.ecs.jsh2.old.MyHelper;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.min;
import static uk.ac.soton.ecs.jsh2.old.MyHelper.distance;

public class LineHolder implements Comparable<LineHolder>{
//    public List<Line2d> lines;
    public Line2d ver1 =null, ver2 =null, hoz1 =null, hoz2 =null;
    public double rank = 0;

    public static int LEFT_INDEX = 0;
    public static int TOP_INDEX = 1;
    public static int RIGHT_INDEX = 2;
    public static int BOTTOM_INDEX = 3;
    public Tetragram tetragram = null;
    public double area;
    private Line2d top, left, right, bottom;
    public double gap;

    public LineHolder(Line2d ver1, Line2d hoz1, Line2d ver2, Line2d hoz2){
        this.ver1 = ver1;
        this.ver2 = ver2;
        this.hoz1 = hoz1;
        this.hoz2 = hoz2;
    }

    public LineHolder(){
//        lines = new ArrayList<>(4);
    }

//    public void setTop(Line2d line){
//        top = line;
//        lines.set(TOP_INDEX, top);
//    }
//
//    public


    public void compute(int width, int height){
        this.tetragram = getBounding2(width, height);
        this.area = new Polygon(tetragram.toList()).calculateArea();
        this.gap = calculateMinLineGap(top, left) + calculateMinLineGap(top, right) + calculateMinLineGap(bottom, left) + calculateMinLineGap(bottom,right);
    }

    double calculateMinLineGap(Line2d l1, Line2d l2){
        double d1 = distance(l1.begin, l2.begin);
        double d2 = distance(l1.begin, l2.end);
        double d3 = distance(l1.end, l2.begin);
        double d4 = distance(l1.end, l2.end);
        d1 = Math.min(d1, d2);
        d3 = Math.min(d3, d4);
        return min(d1, d3);
    }

    @Override
    public int compareTo(LineHolder other) {
        return Double.compare(this.rank, other.rank);
    }

    @Override
    public String toString() {
        return "LineHolder{" +
                "top=" + ver1 +
                ", bottom=" + ver2 +
                ", left=" + hoz1 +
                ", right=" + hoz2 +
                ", rank=" + rank +
                '}';
    }

    //    public Tetragram getBounding(int width, int height) {
//        if(lines.size() < 4) return null;
//        Line2d base1 =  lines.get(0);
//        Line2d fit1 = lines.get(1);
//        Line2d base2 = lines.get(2);
//        Line2d fit2 = lines.get(3);
//
//        Line2d top, left, right, bottom;
//        if(base1.calculateHorizontalAngle() < PI/4){
//            //base1 and base2 are top and bottom
//            if(base1.calculateCentroid().getY() < base2.calculateCentroid().getY()){
//                top = base1;
//                bottom = base2;
//            }else{
//                top = base2;
//                bottom = base1;
//            }
//
//            if(fit1.calculateCentroid().getX() < fit2.calculateCentroid().getX()){
//                left = fit1;
//                right = fit2;
//            }else{
//                left = fit2;
//                right = fit1;
//            }
//
//        }else{
//            //fit1 and fit2 are top and bottom
//            if(fit1.calculateCentroid().getY() < fit2.calculateCentroid().getY()){
//                top = fit1;
//                bottom = fit2;
//            }else{
//                top = fit2;
//                bottom = fit1;
//            }
//
//            if(base1.calculateCentroid().getX() < base2.calculateCentroid().getX()){
//                left = base1;
//                right = base2;
//            }else{
//                left = base2;
//                right = base1;
//            }
//        }
//
//        Point2d topLeft = findLinesIntersection(top, left, width, height);
//        Point2d topRight = findLinesIntersection(top, right, width, height);
//        Point2d bottomRight = findLinesIntersection(bottom, right, width, height);
//        Point2d bottomLeft = findLinesIntersection(bottom, left, width, height);
//
//        return new Tetragram(topLeft, topRight, bottomRight, bottomLeft);
//    }

    public Tetragram getBounding2(int width, int height) {

//        top, left, right, bottom;

        if(hoz1.calculateCentroid().getY() < hoz2.calculateCentroid().getY()){
            top = hoz1;
            bottom = hoz2;
        }else{
            top = hoz2;
            bottom = hoz1;
        }

        if(ver1.calculateCentroid().getX() < ver2.calculateCentroid().getX()){
            left = ver1;
            right = ver2;
        }else{
            left = ver2;
            right = ver1;
        }
        Point2d topLeft = findLinesIntersection(top, left, width, height);
        Point2d topRight = findLinesIntersection(top, right, width, height);
        Point2d bottomRight = findLinesIntersection(bottom, right, width, height);
        Point2d bottomLeft = findLinesIntersection(bottom, left, width, height);

        return new Tetragram(topLeft, topRight, bottomRight, bottomLeft);
    }

    public static Point2d findLinesIntersection(Line2d line1, Line2d line2, int width, int height) {
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

    public List<Line2d> toList() {
        ArrayList<Line2d> lines = new ArrayList<>();
        if(ver1!=null) lines.add(ver1);
        if(ver2!=null) lines.add(ver2);
        if(hoz1!=null) lines.add(hoz1);
        if(hoz2!=null) lines.add(hoz2);
        return lines;
    }
}
