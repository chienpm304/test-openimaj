package uk.ac.soton.ecs.jsh2;

import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class HoughLinesP {
    private List<Point2d> nonZeroPoints;
    private int width;
    private int height;
    private double rho = 1;
    private double theta = Math.PI / 180d;
    private int threshold = 50;
    private int lineLength = 150;
    private int lineGap = 175;
    private int linesMax = 20;
    private int numangle;
    private float[] cosCache;
    private float[] sinCache;
    private boolean[] mask;
    private int[][] accum;
    private int numrho;

    public HoughLinesP(List<Point2d> points,
                       int width,
                       int height,
                       int threshold,
                       int lineLength,
                       int linesMax) {
        this.nonZeroPoints = points;
        this.width = width;
        this.height = height;
        this.threshold = threshold;
        this.lineLength = lineLength;
        this.linesMax = linesMax;

        initialize();
    }

    private void initialize() {
        this.numangle = (int)Math.round(Math.PI / theta);
        this.numrho = (int) Math.round(((width + height) * 2 + 1) / rho);

        this.accum = new int[numangle][];

        this.mask = new boolean[width * height];
        Arrays.fill(this.mask, false);

        //nonZeroPoints
        this.cosCache = new float[numangle];
        this.sinCache = new float[numangle];

        for(int i = 0; i < numangle; i++){
            cosCache[i] = (float) (Math.cos((double)i * theta) / rho);
            sinCache[i] = (float) (Math.sin((double)i * theta) / rho);
        }

        // stage 1. collect non-zero image points
        for(Point2d p: nonZeroPoints){
            int x = (int) p.getX();
            int y = (int) p.getY();
            mask[x + y*width] = true;
        }

        // Shuffle the array randomly
//        Collections.shuffle(nonZeroPoints);


    }

    List<Line2d> getLines(){
        List<Line2d> lines = new ArrayList<>();

        // stage 2. process all the points in random order
        for(Point2d p: nonZeroPoints){
            int col = (int) p.getX();
            int row = (int) p.getY();

            if(!mask[row*width+col]) {
                continue;
            }

            int maxVal = threshold - 1;
            int maxThetaIndex = 0;

            //update accumulator, find the most probable line
            for(int thetaIndex = 0; thetaIndex < numangle; thetaIndex++){
                int rho = (int)Math.round(col * cosCache[thetaIndex] + row * sinCache[thetaIndex]);
                rho += (numrho - 1) / 2;

                if(rho >= numrho)
                    throw new RuntimeException("rho out of index");

                if (accum[thetaIndex]==null) {
                    accum[thetaIndex] = new int[numrho];
                    Arrays.fill(accum[thetaIndex], 0);
                }
                if (accum[thetaIndex][rho] == 0) {
                    accum[thetaIndex][rho] = 1;
                } else {
                    accum[thetaIndex][rho]++;
                }
                int val = accum[thetaIndex][rho];

                if (maxVal < val) {
                    maxVal = val;
                    maxThetaIndex = thetaIndex;
                }
            }

            if(maxVal < threshold)
                continue;

            // from the current point walk in each direction
            // along the found line and extract the line segment
            int[][] lineEnds = new int[2][2];
            int shift = 16;
            double a = -sinCache[maxThetaIndex];
            double b = cosCache[maxThetaIndex];
            int x0 = col;
            int y0 = row;
            int dx0;
            int dy0;
            boolean isWalkingX;
            if (Math.abs(a) > Math.abs(b)) {
                isWalkingX = true;
                dx0 = a > 0 ? 1 : -1;
                dy0 = (int) Math.round(b * (1 << shift) / Math.abs(a));
                y0 = (y0 << shift) + (1 << (shift - 1));
            } else {
                isWalkingX = false;
                dy0 = b > 0 ? 1 : -1;
                dx0 = (int) Math.round(a * (1 << shift) / Math.abs(b));
                x0 = (x0 << shift) + (1 << (shift - 1));
            }

            for(int k = 0; k < 2; k++){
                int gap = 0;
                int x = x0;
                int y = y0;
                int dx = dx0;
                int dy = dy0;

                //Walk in the opposite direction for the second point
                if(k>0){
                    dx = -dx;
                    dy = -dy;
                }

                // walk along the line using fixed-point arithmetics
                for (; ; x += dx, y += dy) {
                    int i1, j1;

                    if (isWalkingX) {
                        j1 = x;
                        i1 = y >> shift;
                    } else {
                        j1 = x >> shift;
                        i1 = y;
                    }

                    // stop at the image border or in case of too big gap
                    if (j1 < 0 || j1 >= width || i1 < 0 || i1 >= height) {
                        break;
                    }

                    // for each non-zero point:
                    //    update line end,
                    //    clear the mask element
                    //    reset the gap
                    if (mask[i1 * width + j1]) {
                        gap = 0;
                        lineEnds[k][0] = j1;
                        lineEnds[k][1] = i1; // x, y of kth point
                    } else if (++gap > lineGap) {
                        break;
                    }
                }
            }

            final boolean
                    goodLine = Math.abs(lineEnds[1][0] - lineEnds[0][0]) >= lineLength ||
                    Math.abs(lineEnds[1][1] - lineEnds[0][1]) >= lineLength;

//            goodLine = new Line2d(lineEnds[0][0], lineEnds[0][1], lineEnds[1][0], lineEnds[1][1]).calculateLength() >= lineLength;

            for (int k = 0; k < 2; k++) {
                int x = x0;
                int y = y0;
                int dx = dx0;
                int dy = dy0;

                if (k > 0) {
                    dx = -dx;
                    dy = -dy;
                }

                // walk along the line using fixed-point arithmetics,
                // stop at the image border or in case of too big gap
                for (; ; x += dx, y += dy) {
                    int i1, j1;

                    if (isWalkingX) {
                        j1 = x;
                        i1 = y >> shift;
                    } else {
                        j1 = x >> shift;
                        i1 = y;
                    }

                    // for each non-zero point:
                    //    update line end,
                    //    clear the mask element
                    //    reset the gap
                    if (mask[i1 * width + j1]) {
                        if (goodLine) {
                            // Since we decided on this line as authentic, remove this pixel's
                            // weights for all possible angles from the accumulator array
                            for (int thetaIndex = 0; thetaIndex < numangle; thetaIndex++) {
                                int rho = (int) Math.round(
                                        j1 * cosCache[thetaIndex] + i1 * sinCache[thetaIndex]
                                );
                                rho += (numrho - 1) / 2;
                                if (accum[thetaIndex]!=null && accum[thetaIndex][rho] != 0) {
                                    accum[thetaIndex][rho]--;
                                }
                            }
                        }

                        mask[i1 * width + j1] = false;
                    }

                    if (i1 == lineEnds[k][1] && j1 == lineEnds[k][0]) break;
                }
            }

            if(goodLine){
                lines.add(new Line2d(lineEnds[0][0], lineEnds[0][1],
                        lineEnds[1][0], lineEnds[1][1]));
                if(lines.size() >= linesMax) {
                    return lines;
                }
            }
        }

        return lines;
    }
}
