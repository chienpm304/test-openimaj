package uk.ac.soton.ecs.jsh2;

import org.openimaj.math.geometry.point.Point2d;

import java.util.List;

public class MyHoughLinesP {

    List<Point2d> points;
    double rho;
    double theta;
    int threshold;
    double minLineLength;
    double maxGap;

    public MyHoughLinesP(List<Point2d> points,
                         double rho,
                         double theta,
                         int threshold,
                         double minLineLength,
                         double maxGap) {
        if(!(rho > 0 && theta >0)){
            throw new RuntimeException("rho and theta must > 0");
        }


        this.points = points;
        this.theta = theta;
        this.theta = threshold;
        this.rho = rho;
        this.minLineLength = minLineLength;
        this.maxGap = maxGap;
    }

    public List<Point2d> getLines(){


        return null;
    }




}
