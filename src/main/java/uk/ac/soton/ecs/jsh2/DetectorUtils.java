package uk.ac.soton.ecs.jsh2;

import org.apache.commons.lang.mutable.MutableFloat;
import org.openimaj.image.MBFImage;
import org.openimaj.image.processing.resize.ResizeProcessor;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;

import java.util.List;

import static java.lang.Math.PI;
import static java.lang.Math.min;

/**
 * Package com.chienpm.library.detector
 * Created by @chienpm on 20/02/2020
 */

public class DetectorUtils {

    public static Line2d.IntersectionResult findIntersection(Line2d l1, Line2d l2) {
        float dx1 = l1.begin.getX() - l1.end.getX();
        float dy1 = l1.begin.getY() - l1.end.getY();
        float ma = l1.begin.getX() * l1.end.getY();
        float mb = l1.end.getX() * l1.begin.getY();

        float dx2 = l2.begin.getX() - l2.end.getX();
        float dy2 = l2.begin.getY() - l2.end.getY();
        float mc = l2.begin.getX() * l2.end.getY();
        float md = l2.end.getX() * l2.begin.getY();

        float du = dx2 * dy1 - dx1 * dy2;
        if (du == 0.0f || (dx1 == 0.0f && dx2 == 0.0f)) {
            return new Line2d.IntersectionResult(Line2d.IntersectionType.NOT_INTERESECTING);
        } else {
            float x = (dx1 * (mc - md) - dx2 * (ma - mb)) / du;
            float y;
            if (dx1 != 0.0f) {
                return new Line2d.IntersectionResult(new Point2dImpl(Math.round(x), Math.round((ma - mb + x * dy1) / dx1)));
            } else if (dx2 != 0) {
                return new Line2d.IntersectionResult(new Point2dImpl(Math.round(x), Math.round((mc - md + x * dy2) / dx2)));
            } else {
                return new Line2d.IntersectionResult(Line2d.IntersectionType.NOT_INTERESECTING);
            }
        }
    }

    static void recoveryOriginalScale(Tetragram tetragram, MutableFloat scaleFactor) {
        float val = (float) scaleFactor.getValue();
        if (val != 1.0f) {
            for (Point2d p : tetragram.toList()) {
                p.setX((int) (p.getX() / val));
                p.setY((int) (p.getY() / val));
            }
        }
    }

    static double getHorizontalAngleInDegree(Line2d l1) {
        return calcAngleDiffInDegree(getAngleInDegree(l1), 0);
    }

    static float calculateLineGap(Line2d l1, Line2d l2) {

        if(l1.isInLine(l2.begin, Constants.MERGE_MAX_LINE_DISTANCE)
                || l1.isInLine(l2.end, Constants.MERGE_MAX_LINE_DISTANCE)
                || l2.isInLine(l1.begin, Constants.MERGE_MAX_LINE_DISTANCE)
                || l2.isInLine(l1.end, Constants.MERGE_MAX_LINE_DISTANCE))
            return 0;

        float d1 = (float) Line2d.distance(l1.begin, l2.begin);
        float d2 = (float) Line2d.distance(l1.begin, l2.end);
        float d3 = (float) Line2d.distance(l1.end, l2.begin);
        float d4 = (float) Line2d.distance(l1.end, l2.end);
        d1 = Math.min(d1, d2);
        d3 = Math.min(d3, d4);
        return min(d1, d3);
    }

    static double getAngleInDegree(Line2d line) {
        return line.calculateHorizontalAngle() * 180 / PI;
    }

    static double calcHorizontalAngleDiff(Line2d l1, Line2d l2){
        double b1 = getHorizontalAngleInDegree(l1);
        double b2 = getHorizontalAngleInDegree(l2);
        return Math.abs(b1 - b2);
    }

    static boolean checkIfMayOnTheSameLine(List<Line2d> lines, Line2d line, int maxLineDistance) {
        for(Line2d l: lines){
            if(!isOnTheSameLine(l, line, maxLineDistance))
                return false;
        }
        return true;
    }

    private static boolean isOnTheSameLine(Line2d l1, Line2d l2, int threshold) {
        return
                l1.isOnLine(l2.begin, threshold)
                        && l1.isOnLine(l2.end, threshold)
                        && l2.isOnLine(l1.begin, threshold)
                        && l2.isOnLine(l1.end, threshold)
                        && calcAngleDiffInDegree(l1.calculateHorizontalAngle(), l2.calculateHorizontalAngle()) <= Constants.MIN_ANGLE;
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

    static void sortByYAxis(Line2d line) {
        if (line.begin.getY() > line.end.getY()) {
            Point2d tmp = line.begin;
            line.setBeginPoint(line.end);
            line.setEndPoint(tmp);
        }
    }

    static void sortByXAxis(Line2d line) {
        if (line.begin.getX() > line.end.getX()) {
            Point2d tmp = line.begin;
            line.setBeginPoint(line.end);
            line.setEndPoint(tmp);
        }
    }

}