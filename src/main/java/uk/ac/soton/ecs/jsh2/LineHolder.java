package uk.ac.soton.ecs.jsh2;

import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;

import javax.sound.sampled.Line;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.zip.DeflaterOutputStream;

import static java.lang.StrictMath.PI;
import static uk.ac.soton.ecs.jsh2.MyUtils.findLinesIntersection;

public class LineHolder implements Comparable<LineHolder>{
    public List<Line2d> lines;
    public double rank = 0;
    public LineHolder(Line2d line1, Line2d line2, Line2d line3, Line2d line4){
        lines = new ArrayList<>();
        lines.add(line1);
        lines.add(line2);
        lines.add(line3);
        lines.add(line4);
    }

    public LineHolder(){
        lines = new ArrayList<>();
    }

    @Override
    public int compareTo(LineHolder other) {
        return Double.compare(this.rank, other.rank);
    }

    @Override
    public String toString() {
        return "LineHolder{" +
                "lines=" + lines +
                ", rank=" + rank +
                '}';
    }

    public Tetragram getBounding(int width, int height) {
        if(lines.size() < 4) return null;
        Line2d base1 = lines.get(0);
        Line2d fit1 = lines.get(1);
        Line2d base2 = lines.get(2);
        Line2d fit2 = lines.get(3);
        Line2d top, left, right, bottom;
        if(base1.calculateHorizontalAngle() < PI/4){
            //base1 and base2 are top and bottom
            if(base1.calculateCentroid().getY() < base2.calculateCentroid().getY()){
                top = base1;
                bottom = base2;
            }else{
                top = base2;
                bottom = base1;
            }

            if(fit1.calculateCentroid().getX() < fit2.calculateCentroid().getX()){
                left = fit1;
                right = fit2;
            }else{
                left = fit2;
                right = fit1;
            }

        }else{
            //fit1 and fit2 are top and bottom
            if(fit1.calculateCentroid().getY() < fit2.calculateCentroid().getY()){
                top = fit1;
                bottom = fit2;
            }else{
                top = fit2;
                bottom = fit1;
            }

            if(base1.calculateCentroid().getX() < base2.calculateCentroid().getX()){
                left = base1;
                right = base2;
            }else{
                left = base2;
                right = base1;
            }
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

}
