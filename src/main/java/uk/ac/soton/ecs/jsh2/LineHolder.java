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

public class LineHolder implements Comparable<LineHolder> {
    //    public List<Line2d> lines;
    public Line2d left = null, right = null, top = null, bottom = null;
    public double rank = 0;

    public static int LEFT_INDEX = 0;
    public static int TOP_INDEX = 1;
    public static int RIGHT_INDEX = 2;
    public static int BOTTOM_INDEX = 3;
    public Tetragram tetragram = null;
    public double area;
    public double gap;

    public LineHolder(Line2d left, Line2d top, Line2d right, Line2d bottom) {
        this.left = left;
        this.right = right;
        this.top = top;
        this.bottom = bottom;
    }

    public LineHolder() {
//        lines = new ArrayList<>(4);
    }

//    public void setTop(Line2d line){
//        top = line;
//        lines.set(TOP_INDEX, top);
//    }
//
//    public


    public void compute(int width, int height) {
        this.tetragram = getBounding2(width, height);
        this.area = new Polygon(tetragram.toList()).calculateArea();

//        this.gap =
//                calculateMinLineGap(top, left)
//                + calculateMinLineGap(top, right)
//                + calculateMinLineGap(bottom, left)
//                + calculateMinLineGap(bottom, right);
        this.gap = distance(tetragram.getTopLeft(), top.begin)
                + distance(tetragram.getTopLeft(), left.begin)
                + distance(tetragram.getTopRight(), top.end)
                + distance(tetragram.getTopRight(), right.begin)
                + distance(tetragram.getBottomRight(), right.end)
                + distance(tetragram.getBottomRight(), bottom.end)
                + distance(tetragram.getBottomLeft(), bottom.begin)
                + distance(tetragram.getBottomLeft(), left.end);
    }

    private double calculateMinLineGap(Line2d l1, Line2d l2) {
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
        StringBuilder st = new StringBuilder("LineHolder{ [");
        if(left!=null)
            st.append(Math.round(DetectorUtils.getUnsignedAngleInDegree(left)));
        st.append(", ");
        if(right!=null)
            st.append(Math.round(DetectorUtils.getUnsignedAngleInDegree(right)));
        st.append("] [");
        if(top!=null)
            st.append(Math.round(DetectorUtils.getUnsignedAngleInDegree(top)));
        st.append(", ");
        if(bottom!=null)
            st.append(Math.round(DetectorUtils.getUnsignedAngleInDegree(bottom)));
        st.append("]}");
        return st.toString();
    }

    public Tetragram getBounding2(int width, int height) {

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

    public List<Line2d> toRawList() {
        ArrayList<Line2d> lines = new ArrayList<>();
        if (left != null) lines.add(left);
        if (right != null) lines.add(right);
        if (top != null) lines.add(top);
        if (bottom != null) lines.add(bottom);
        return lines;
    }

    public boolean isSatisfyingShape() {
        if(!checkGap()){
            System.out.println("rejected by check gap");
            return false;
        }
//        if(!checkRawLines())
//            return false;

        //check size
        List<Line2d> lines = tetragram.toLineList(); // top, right, bot, left
        double top = lines.get(0).calculateLength();
        double right = lines.get(1).calculateLength();
        double bottom = lines.get(2).calculateLength();
        double left = lines.get(3).calculateLength();

        if (top / bottom > 3f || bottom / top > 3f) {
            System.out.println("rejected by check top/bottom ratio");
            return false;
        }
        if (left / right > 3f || right / left > 3f) {
            System.out.println("rejected by check left/right ratio");
            return false;
        }

        double a1, a2, gapAngle;

        //check angles
        for (int i = 0; i < 4; i++) {
            a1 = DetectorUtils.getSignedAngleInDegree(lines.get(i % 4));
            a2 = DetectorUtils.getSignedAngleInDegree(lines.get((i + 1) % 4));
            gapAngle = DetectorUtils.calcAngleDiffInDegree(a1, a2);// Math.abs(a1 - a2);
            if(gapAngle < 30) {
                System.out.println("rejected by check angle");
                return false;
            }
        }

        //check convex
        if(!new Polygon(tetragram.toList()).isConvex()) {
            System.out.println("rejected by check convex");
            return false;
        }

        return true;
    }

    private boolean checkGap() {
        double perimeter =
                left.calculateLength()
                + right.calculateLength()
                + top.calculateLength()
                + bottom.calculateLength();
        return gap < perimeter / 3f;
    }

    private boolean checkRawLines() {

        double verticalDistance = left.distanceToLine(right.calculateCentroid());
        double horizontalDistance = top.distanceToLine(bottom.calculateCentroid());

        // check horizontal lines are located at the same side of vertical lines
        if (left.distanceToLine(top.begin) > verticalDistance &&
                left.distanceToLine(top.end) > verticalDistance)
            return false;

        if (left.distanceToLine(bottom.begin) > verticalDistance &&
                left.distanceToLine(bottom.end) > verticalDistance)
            return false;

        if (right.distanceToLine(top.begin) > verticalDistance &&
                right.distanceToLine(top.end) > verticalDistance)
            return false;

        if (right.distanceToLine(bottom.begin) > verticalDistance &&
                right.distanceToLine(bottom.end) > verticalDistance)
            return false;


        // check top and bottom length
//        if (top.calculateLength() < verticalDistance / 2f &&
//                bottom.calculateLength() < verticalDistance / 2f)
//            return false;

        Point2d topPoint = left.begin.getY() < right.begin.getY() ? left.begin : right.begin;
        Point2d bottomPoint = left.end.getY() > right.end.getY() ? left.end : right.end;

        // check top and bottom lines are located at 2 separated position
        float select_line_const = 30;
        if ((top.distanceToLine(topPoint) < select_line_const &&
                bottom.distanceToLine(bottomPoint) < select_line_const)
                ||
                (top.distanceToLine(bottomPoint) < select_line_const &&
                        bottom.distanceToLine(topPoint) < select_line_const)
        ) {
            return true;
        }

        return false;
    }
}
