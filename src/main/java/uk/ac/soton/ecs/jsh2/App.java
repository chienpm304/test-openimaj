package uk.ac.soton.ecs.jsh2;

import org.apache.logging.log4j.core.net.TcpSocketManager;
import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.processing.algorithm.GammaCorrection;
import org.openimaj.image.processing.convolution.FGaussianConvolve;
import org.openimaj.image.processing.edges.CannyEdgeDetector;
import org.openimaj.image.processing.resize.ResizeProcessor;
import org.openimaj.image.typography.hershey.HersheyFont;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;

import static java.lang.Math.*;

/**
 * OpenIMAJ Hello world!
 */
public class App {


    public static float scaleFactor = 1.0f;
    static int height;
    static int width;


    public static void main(String[] args) throws IOException {
        testDetectBox();
    }


    private static void testDetectBox() throws IOException {
        File fin = new File(Constants.LINUX_DIR_IN);
        File fout = new File(Constants.LINUX_DIR_OUT);
        if (fin.exists() && fin.isDirectory())
            for (final File file : fin.listFiles()) {
                if (file.isFile())
                    detectWithCanny(file, fout);
            }
    }

    static FImage applyCannyDetector(FImage grey) {
        CannyEdgeDetector canny = new CannyEdgeDetector(Constants.CANNY_LOW_THRESH, Constants.CANNY_HIGH_THRESH, Constants.CANNY_SIGMA);
        grey.processInplace(new FGaussianConvolve(Constants.GAUSSIAN_BLUR_SIGMA));
        canny.processImage(grey);
        return grey;
    }

    private static Tetragram detectWithCanny(File fin, File fout) throws IOException {

        MBFImage frame = ImageUtilities.readMBF(fin);
        if(frame.getWidth() > frame.getHeight())
            scaleFactor = Constants.STANDARD_WIDTH / (float) frame.getHeight();
        else
            scaleFactor = Constants.STANDARD_WIDTH / (float) frame.getWidth();

        scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;

        if(scaleFactor!=1.0f) {
            frame.processInplace(new ResizeProcessor(scaleFactor));
        }

        width = frame.getWidth();
        height = frame.getHeight();

        System.out.println("processing: " + fin.getName() + "...");

        GammaCorrection gc = new GammaCorrection(Constants.GAMMA);
        for (int i = 0; i < frame.numBands(); i++) {
            frame.getBand(i).processInplace(gc);
        }

        Point2dImpl center = new Point2dImpl(width / 2, height / 2);

        FImage edges = applyCannyDetector(frame.flattenMax());

        ImageUtilities.write(edges, new File(fout.getAbsolutePath() + "/edged/" + fin.getName()));

        List<Line2d> lines = getLinesUsingHoughTransformP(edges);

        MBFImage tmp = frame.clone();
        drawLines(tmp,center, lines, RGBColour.GREEN, false);
        ImageUtilities.write(tmp, new File(fout.getAbsolutePath() + "/raw_detected/" + fin.getName()));


        System.out.println("Before merge: " + lines.size());

        removeNoiseLines(lines, Constants.MERGE_MAX_LINE_DISTANCE, Constants.MERGE_MAX_LINE_GAP);
//        removeNoiseLines(lines, Constants.MERGE_MAX_LINE_DISTANCE/2, Constants.MERGE_MAX_LINE_GAP/3);

        System.out.println("After merge: " + lines.size());

        drawLines(frame, center, lines, RGBColour.GREEN, true);

//        List<LineHolder> results = findBounds(lines);
//
//        if (!results.isEmpty())
//            drawBound(frame, center, results.get(0).lines, RGBColour.GREEN, RGBColour.YELLOW);


        ImageUtilities.write(frame, new File(fout.getAbsolutePath() + "/out/" + fin.getName()));
        return null;
    }

    static List<LineHolder> findBounds(List<Line2d> lines) {

        ArrayList<LineHolder> res = new ArrayList<>();

        Line2d base1, base2, fit1, fit2;
        double b1, b2, f1, f2;
        int n = lines.size();
        boolean addDummyLine = false;
        int angle_step = 5;
        boolean detected = false;

        while (res.isEmpty() && angle_step <= 45) {
            for (int i = 0; i < n - 1; i++) {
                base1 = lines.get(i);
                b1 = getAngleInDegree(base1);
                for (int j = i + 1; j < n; j++) {
                    base2 = lines.get(j);
                    b2 = getAngleInDegree(base2);
                    detected = false;
                    if (calcAngleDiffInDegree(b1, b2) <= angle_step) {

                        for (int k = 0; k < lines.size() - 1; k++) {
                            if (k == i || k == j) continue;

                            fit1 = lines.get(k);

                            for (int l = k + 1; l < lines.size(); l++) {
                                if (l == i || l == j) continue;
                                fit2 = lines.get(l);

                                LineHolder lh = considerToAddBound(base1, base2, fit1, fit2, angle_step);
                                if (lh != null) {
                                    res.add(lh);
                                    detected = true;
                                }
                            }
                        }

                        if (!detected) {
                            Line2d tmp1 = base1.clone();
                            Line2d tmp2 = base2.clone();
                            // added 2 dummy lines
                            if (b1 < 45) {
                                sortByXAxis(tmp1);
                                sortByXAxis(tmp2);
                            } else {
                                sortByYAxis(tmp1);
                                sortByYAxis(tmp2);
                            }

                            fit1 = new Line2d(base1.begin, base2.begin);
                            fit2 = new Line2d(base1.end, base2.end);

                            LineHolder lh = considerToAddBound(base1, base2, fit1, fit2, angle_step);
                            if (lh != null) {
                                res.add(lh);
                                detected = true;
                            }
                        }
                    }
                }
            }
            angle_step += 5;
        }
        Collections.sort(res);
        System.out.println("Detected " + res.size() + " bounds");
//        for (LineHolder lh : res) {
//            System.out.println(lh.toString());
//        }


        return res;
    }

    private static LineHolder considerToAddBound(Line2d base1, Line2d base2, Line2d fit1, Line2d fit2, double angle_step) {
        double f1 = getAngleInDegree(fit1);
        double f2 = getAngleInDegree(fit2);

        double b1 = getAngleInDegree(base1);
        double b2 = getAngleInDegree(base2);

        double bases_distance = (base1.distanceToLine(base2.begin) + base1.distanceToLine(base2.end)) / 2d;

        if (calcAngleDiffInDegree(b1, f1) >= 90 - angle_step * 2 && calcAngleDiffInDegree(b2, f1) >= 90 - angle_step * 2) {
            if (calcAngleDiffInDegree(b1, f2) >= 90 - angle_step * 2 && calcAngleDiffInDegree(b2, f2) >= 90 - angle_step * 2) {
                LineHolder lh = new LineHolder(base1, fit1, base2, fit2);
                double base1_fit1 = Math.min(fit1.distanceToLine(base1.begin), fit1.distanceToLine(base1.end));
                double base1_fit2 = Math.min(fit2.distanceToLine(base1.begin), fit2.distanceToLine(base1.end));
                double base2_fit1 = Math.min(fit1.distanceToLine(base2.begin), fit1.distanceToLine(base2.end));
                double base2_fit2 = Math.min(fit2.distanceToLine(base2.begin), fit2.distanceToLine(base2.end));

                double fit1_base1 = Math.min(base1.distanceToLine(fit1.begin), base1.distanceToLine(fit1.end));
                double fit1_base2 = Math.min(base2.distanceToLine(fit1.begin), base2.distanceToLine(fit1.end));
                double fit2_base1 = Math.min(base1.distanceToLine(fit2.begin), base1.distanceToLine(fit2.end));
                double fit2_base2 = Math.min(base2.distanceToLine(fit2.begin), base2.distanceToLine(fit2.end));

                double fits_distance = (fit1.distanceToLine(fit2.begin) + fit1.distanceToLine(fit2.end)) / 2d;


                lh.rank =
                        base1_fit1 + base1_fit2 + base2_fit1 + base2_fit2
                                - fits_distance - bases_distance +
                                fit1_base1 + fit1_base2 + fit2_base1 + fit2_base2;

//                                    lh.rank = Math.min(base1_fit1, fit1_base1) +
//                                            Math.min(base1_fit2, fit1_base2) +
//                                            Math.min(base2_fit1, fit2_base1) +
//                                            Math.min(base2_fit2, fit2_base2);

                return lh;
            }
        }
        return null;
    }

    private static List<Line2d> getLinesUsingHoughTransformP(FImage image) {
        List<Point2d> points = new ArrayList<>();
        for (int i = 0; i < image.width; i++)
            for (int j = 0; j < image.height; j++) {
                if (image.getPixel(i, j) > 0.5f)
                    points.add(new Point2dImpl(i, j));
            }
        HoughLinesP ht = new HoughLinesP(
                points,
                width,
                height,
                Constants.HOUGH_LINE_RHO,
                Constants.HOUGH_LINE_THETA,
                Constants.HOUGH_LINE_THRESHOLD,
                Constants.HOUGH_LINE_MAX_LINE_GAP,
                Constants.HOUGH_MIN_LINE_LENGTH,
                Constants.HOUGH_MAX_NUM_LINE);
        return ht.getLines();
    }

    static void removeNoiseLines(List<Line2d> lines, int maxLineDistance, int maxLineGap) {

        removeLinesNearbyBounding(lines, Constants.BOUNDING_GAP_REMOVAL);

//        Collections.sort(lines, Comparator.comparingDouble(Line2d::calculateLength).reversed());

        removeLineByMergingX(lines, maxLineDistance, maxLineGap);
        removeLineByMergingY(lines, maxLineDistance, maxLineGap);
//        removeNoiseX(lines, maxLineDistance, maxLineGap);
//        removeNoiseY(lines, maxLineDistance, maxLineGap);

        removeLineWithShorterThanThreshold(lines, Constants.REMOVE_AFTER_MERGE_THRESHOLD);
//        removeLineNotMakeSquareAngle(lines);
    }

    private static void removeLineNotMakeSquareAngle(List<Line2d> lines) {
        int count;
        for(int i = 0; i < lines.size(); i++){
            count = 0;
            for(int j = 0; j < lines.size(); j++){
                if(i==j) continue;;
                double f1 = getHorizontalAngleInDegree(lines.get(i));
                double f2 = getHorizontalAngleInDegree(lines.get(j));
                if(Math.abs(f1-f2) >= 60){
                    count++;
                }
            }
            if(count < 2){
                lines.remove(i);
                i--;
            }

        }

    }

    private static void removeLineByMergingX(List<Line2d> lines, int maxLineDistance, int maxLineGap) {
        for(Line2d l: lines) sortByXAxis(l);
        lines.sort(Comparator.comparingDouble(Line2d::calculateLength).reversed());

        Line2d l1, l2;
        int beforeMerge, afterMerge;
        List<Line2d> tmpLines = new ArrayList<>();
        System.out.print("Merging: ");
        for(int i = 0; i < lines.size(); i++){
            tmpLines.clear();
            l1 = lines.get(i);
            if(getHorizontalAngleInDegree(l1) >= 45)
                continue;
            tmpLines.add(l1);
            for(int j = i + 1; j < lines.size(); j++){
                l2 = lines.get(j);
                if(checkIfMayOnTheSameLine(tmpLines, l2, maxLineDistance)){
                    tmpLines.add(l2);
                }
            }
            if(tmpLines.size() > 1){
                lines.removeAll(tmpLines);
                System.out.print(" "+tmpLines.size());
                beforeMerge = tmpLines.size();
                mergeLines2(tmpLines, maxLineGap);
                afterMerge = tmpLines.size();
                lines.addAll(tmpLines);
                if(afterMerge < beforeMerge)
                    i--;
            }
        }
        System.out.println("");
    }
    private static void removeLineByMergingY(List<Line2d> lines, int maxLineDistance, int maxLineGap) {
        for(Line2d l: lines) sortByYAxis(l);
        lines.sort(Comparator.comparingDouble(Line2d::calculateLength).reversed());

        Line2d l1, l2;
        int beforeMerge, afterMerge;
        List<Line2d> tmpLines = new ArrayList<>();
        System.out.print("Merging: ");
        for(int i = 0; i < lines.size(); i++){
            tmpLines.clear();
            l1 = lines.get(i);
            if(getHorizontalAngleInDegree(l1) < 45)
                continue;
            tmpLines.add(l1);
            for(int j = i + 1; j < lines.size(); j++){
                l2 = lines.get(j);
                if(checkIfMayOnTheSameLine(tmpLines, l2, maxLineDistance)){
                    tmpLines.add(l2);
                }
            }
            if(tmpLines.size() > 1){
                lines.removeAll(tmpLines);
                System.out.print(" "+tmpLines.size());
                beforeMerge = tmpLines.size();
                mergeLines2(tmpLines, maxLineGap);
                afterMerge = tmpLines.size();
                lines.addAll(tmpLines);
                if(afterMerge < beforeMerge)
                    i--;
            }
        }
        System.out.println("");
    }

    private static void mergeLines2(List<Line2d> lines, int maxLineGap) {
        boolean mergeX = getHorizontalAngleInDegree(lines.get(0)) < 45;
        if(mergeX){
            for(Line2d l: lines) sortByXAxis(l);
            lines.sort(Comparator.comparingDouble(Line2d::calculateLength).reversed());
        }else{
            for(Line2d l: lines) sortByYAxis(l);
            lines.sort(Comparator.comparingDouble(Line2d::calculateLength).reversed());
        }

//        Line2d candidateLine = findCandidateLineByMinimumLineDistance(lines).clone();

        Line2d l1, l2;
        int idx;
        for(int i = 0; i < lines.size(); i++){
            l1 = lines.get(i);
            while((idx = findClosestLineByLineGap(lines, l1, maxLineGap)) != -1){
                if(i>=lines.size() || idx >= lines.size())
                    break;;

                l2 = lines.get(idx);

//                if(sumLinesDistance(l1, candidateLine) < sumLinesDistance(l2, candidateLine)){
                if(l1.calculateLength() < l2.calculateLength()){
                    Collections.swap(lines, i, idx);
                    l1 = lines.get(i); // keep this
                    l2 = lines.get(idx);
                }


                if(mergeX) mergeX(l1, l2);
                else mergeY(l1, l2);

                lines.remove(idx);

            }
        }


    }

    private static double sumLinesDistance(Line2d l1, Line2d l2){
        return l1.distanceToLine(l2.begin) + l1.distanceToLine(l2.end);
    }

    private static Line2d findCandidateLineByMinimumLineDistance(List<Line2d> lines) {
        int idx = -1;
        double dist;
        double minDistance = Double.MAX_VALUE;

        Line2d l1, l2;
        for(int i = 0; i < lines.size(); i++){
            dist = 0;
            l1 = lines.get(i);
            for(int j = 0; j < lines.size(); j++){
                if (i==j) continue;
                l2 = lines.get(j);
                dist += l1.distanceToLine(l2.begin) + l1.distanceToLine(l2.end);
            }
            if(dist < minDistance){
                minDistance = dist;
                idx = i;
            }
        }
        return lines.get(idx);
    }

    private static int findClosestLineByLineGap(List<Line2d> lines, Line2d line, int maxLineGap) {
        double minGap = maxLineGap + 1;
        int idx = -1;
        Line2d l1, l2;
        for(int i = 0; i < lines.size(); i++){
            l1 = lines.get(i);
            if(l1!=line){
                double gap = calculateLineGap(l1, line);
                if(gap < minGap && gap < maxLineGap){
                    minGap = gap;
                    idx = i;
                }
            }
        }
        return idx;
    }

    private static boolean checkIfMayOnTheSameLine(List<Line2d> lines, Line2d line, int maxLineDistance) {
        for(Line2d l: lines){
            if(!isOnTheSameLine(l, line, maxLineDistance))
                return false;
        }
        return true;
    }

    private static void removeNoiseX(List<Line2d> lines, int maxLineDistance, int maxLineGap) {
        for(Line2d l: lines) sortByXAxis(l);
        lines.sort(Comparator.comparingDouble(Line2d::minX));

        Line2d l1, l2;
        List<Line2d> tmpLines = new ArrayList<>();

        for (int i = 0; i < lines.size() - 1; i++) {
            tmpLines.clear();
            l1 = lines.get(i);
            if(getHorizontalAngleInDegree(l1) >= 45)
                continue;
            tmpLines.add(l1);

            for (int j = i + 1; j < lines.size(); j++) {
                l2 = lines.get(j);

                if (l1 != l2 && isTheSameLineGroup(tmpLines, l2, maxLineDistance, maxLineGap)) {
                    l2 = lines.remove(j);
                    tmpLines.add(l2);
                    j--;
                }
            }
            if(tmpLines.size() > 1){
                System.out.println("merging "+tmpLines.size());
                Line2d newLine = mergeLines(tmpLines);
                lines.set(i, newLine);
            }
        }
    }

    private static void removeNoiseY(List<Line2d> lines, int maxLineDistance, int maxLineGap) {
        for(Line2d l: lines) sortByYAxis(l);
        lines.sort(Comparator.comparingDouble(Line2d::minY));

        Line2d l1, l2;
        List<Line2d> tmpLines = new ArrayList<>();

        for (int i = 0; i < lines.size() - 1; i++) {
            tmpLines.clear();
            l1 = lines.get(i);
            if(getHorizontalAngleInDegree(l1) < 45)
                continue;
            tmpLines.add(l1);

            for (int j = i + 1; j < lines.size(); j++) {
                l2 = lines.get(j);

                if (l1 != l2 && isTheSameLineGroup(tmpLines, l2, maxLineDistance, maxLineGap)) {
                    l2 = lines.remove(j);
                    tmpLines.add(l2);
                    j--;
                }
            }
            if(tmpLines.size() > 1){
                System.out.println("merging "+tmpLines.size());
                Line2d newLine = mergeLines(tmpLines);
                lines.set(i, newLine);
            }
        }
    }

    private static double getHorizontalAngleInDegree(Line2d l1) {
        return calcAngleDiffInDegree(getAngleInDegree(l1), 0);
    }

    private static Line2d mergeLines(List<Line2d> lines) {
        Line2d keep, remove;
        //find the candidate line in lines
        int idx = findClosestLineByAverageAngle(lines);

        keep = lines.remove(idx);

        for(int i = 0; i < lines.size(); i++){
            remove = lines.get(i);
            mergeLine(keep, remove);
        }
        return keep;
    }

    private static int findClosestLineByAverageAngle(List<Line2d> lines) {
        double avgAngle = 0;
        for(Line2d l: lines)
            avgAngle+= getHorizontalAngleInDegree(l);

        avgAngle = avgAngle/lines.size();

        int idx = 0;

        double minAngleGap = Math.abs(avgAngle - getHorizontalAngleInDegree(lines.get(idx)));

        for(int i = 1; i < lines.size(); i++){
            double t = Math.abs(avgAngle - getHorizontalAngleInDegree(lines.get(i)));
            if(t < minAngleGap){
                minAngleGap = t;
                idx = i;
            }
        }
        return idx;
    }


    private static boolean isTheSameLineGroup(List<Line2d> lines, Line2d line, int lineDistanceThreshold, int lineGapThreshold) {
        boolean gap = false;
        for(Line2d l: lines){
            if(!isOnTheSameLine(l, line, lineDistanceThreshold))
                return false;
            if(calculateLineGap(l, line) < lineGapThreshold)
                gap = true;
        }
        return gap;
    }

    private static float calculateLineGap(Line2d l1, Line2d l2) {
        float d1 = distance(l1.begin, l2.begin);
        float d2 = distance(l1.begin, l2.end);
        float d3 = distance(l1.end, l2.begin);
        float d4 = distance(l1.end, l2.end);
        d1 = Math.min(d1, d2);
        d3 = Math.min(d3, d4);
        return min(d1, d3);
    }

    private static boolean isOnTheSameLine(Line2d l1, Line2d l2, int threshold) {
        return
                   l1.isOnLine(l2.begin, threshold)
                && l1.isOnLine(l2.end, threshold)
                && l2.isOnLine(l1.begin, threshold)
                && l2.isOnLine(l1.end, threshold)
                && calcAngleDiffInDegree(l1.calculateHorizontalAngle(), l2.calculateHorizontalAngle()) <= Constants.MIN_ANGLE;
    }


    static void removeSimilarAndNoiseLines(List<Line2d> lines) {

//        removeLineWithShorterThanTheshold(lines, 40);

        removeLinesNearbyBounding(lines, Constants.BOUNDING_GAP_REMOVAL);

//        merge

        Line2d l1, l2;
        Line2d kpLine, rmLine;
        for (int i = 0; i < lines.size() - 1; i++) {
            l1 = lines.get(i);
            for (int j = i + 1; j < lines.size(); j++) {
                l2 = lines.get(j);
                if (l1 != l2
                        && l2.distanceToLine(l1.begin) < Constants.MERGE_MAX_LINE_DISTANCE
                        && l2.distanceToLine(l1.end) < Constants.MERGE_MAX_LINE_DISTANCE
                        && l1.distanceToLine(l2.begin) < Constants.MERGE_MAX_LINE_DISTANCE
                        && l1.distanceToLine(l2.end) < Constants.MERGE_MAX_LINE_DISTANCE
                ) {
                    if (l1.calculateLength() < l2.calculateLength()) { // keep line i
                        Collections.swap(lines, i, j);
                        rmLine = lines.remove(j);
                        kpLine = lines.get(i);
                        j = i;

                    } else {
                        rmLine = lines.remove(j);
                        kpLine = lines.get(i);
                        j--;
                    }

                    //merge rmLine to kpLine
                    mergeLine(kpLine, rmLine);
                    lines.set(i, kpLine);
                    l1 = lines.get(i);
                }
            }
        }

        removeLineWithShorterThanThreshold(lines, 80);

    }

    private static void removeLinesNearbyBounding(List<Line2d> lines, int threshold) {
        for (int i = 0; i < lines.size(); i++) {
            Line2d l = lines.get(i);

            if (l.begin.getX() < threshold && l.end.getX() < threshold //LEFT

                    || l.begin.getY() < threshold && l.end.getY() < threshold //TOP

                    || l.begin.getX() > width - threshold && l.end.getX() > width - threshold // RIGHT

                    || l.begin.getY() > height - threshold && l.end.getY() > height - threshold) {

                lines.remove(i);
                i--;
            }
        }
    }

    private static void removeLineWithShorterThanThreshold(List<Line2d> lines, int threshold) {
        //remove
        for (int i = 0; i < lines.size(); i++) {
            Line2d l = lines.get(i);

            if (l.calculateLength() < threshold) {

                lines.remove(i);
                i--;
            }
        }
    }

    static void mergeLine(Line2d keep, Line2d remove) {
        //determinate the merge axis
        double hAngle = getHorizontalAngleInDegree(keep);

        Point2d expandedPoint;
        if (hAngle < 45) {//x-axis
            mergeX(keep, remove);

        } else { // y-axis
            mergeY(keep, remove);
        }
    }

    private static void mergeY(Line2d keep, Line2d remove) {
        // begin.x < end.x
        sortByYAxis(keep);
        sortByYAxis(remove);
        float a = keep.end.getY() - keep.begin.getY();
        float b = keep.end.getX() - keep.begin.getX();
        float c = keep.end.getX() * keep.begin.getY() - keep.begin.getX() * keep.end.getY();

        float newX, newY;

        if (remove.begin.getY() < keep.begin.getY()) {
            // expand top
            newY = remove.begin.getY();
            newX = (b * newY - c) / (a);

            if (newX < 0) {
                newX = 0;
                newY = (a * newX + c) / (b);
            } else if (newX > width) {
                newX = width - 1;
                newY = (b * newY - c) / (a);
            }

            keep.begin.setX(newX);
            keep.begin.setY(newY);
        } else if (remove.end.getY() > keep.end.getY()) {
            // expand bottom
            newY = remove.end.getY();
            newX = (b * newY - c) / (a);

            if (newX < 0) {
                newX = 0;
                newY = (a * newX + c) / (b);
            } else if (newX > width) {
                newX = width - 1;
                newY = (b * newY - c) / (a);
            }

            keep.end.setX(newX);
            keep.end.setY(newY);
        }
    }

    private static void mergeX(Line2d keep, Line2d remove) {
        sortByXAxis(keep);
        sortByXAxis(remove);
        float a = keep.end.getY() - keep.begin.getY();
        float b = keep.end.getX() - keep.begin.getX();
        float c = keep.end.getX() * keep.begin.getY() - keep.begin.getX() * keep.end.getY();

        float newX, newY;

        if (remove.begin.getX() < keep.begin.getX()) {
            // expand left
            newX = remove.begin.getX();
            newY = (a * newX + c) / (b);

            if (newY < 0) {
                newY = 0;
                newX = (b * newY - c) / (a);
            } else if (newY > height) {
                newY = height - 1;
                newX = (b * newY - c) / (a);
            }

            keep.begin.setX(newX);
            keep.begin.setY(newY);
        } else if (remove.end.getX() > keep.end.getX()) {
            // expand right
            newX = remove.end.getX();
            newY = (a * newX + c) / (b);

            if (newY < 0) {
                newY = 0;
                newX = (b * newY - c) / (a);
            } else if (newY > height) {
                newY = height - 1;
                newX = (b * newY - c) / (a);
            }

            keep.end.setX(newX);
            keep.end.setY(newY);
        }
    }


    static Random rnd = new Random();

    private static Float[] getRandomColor() {
        return RGBColour.RGB(rnd.nextInt(256), rnd.nextInt(256), rnd.nextInt(256));
    }

    static float distance(Point2d p1, Point2d p2) {

        float x1 = p1.getX(), x2 = p2.getX(), y1 = p1.getY(), y2 = p2.getY();

        return (float) Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }

    static void drawBound(MBFImage frame, Point2d center, List<Line2d> lines, Float[] baseColor, Float[] fitColor) {
        Float[] lineColor;
        for (Line2d line : lines) {
            if (lines.indexOf(line) % 2 == 0)
                lineColor = baseColor;
            else
                lineColor = fitColor;

            frame.drawLine(line, 2, lineColor);
            frame.drawPoint(line.begin, RGBColour.RED, 5);
            frame.drawPoint(line.end, RGBColour.YELLOW, 5);
        }
    }

    protected static void drawLines(MBFImage frame, Point2d center, List<Line2d> lines, Float[] lineColor, boolean drawAngle) {
        System.out.println(lines.size() + " lines");

        Float[] orgColor = lineColor;

        for (Line2d line : lines) {
            frame.drawLine(line, 2, lineColor);
            frame.drawPoint(line.begin, RGBColour.RED, 3);
            frame.drawPoint(line.end, RGBColour.YELLOW, 4);
            if(drawAngle)
                frame.drawText(
                    Math.round(getHorizontalAngleInDegree(line)) + "*",
                    (int) line.calculateCentroid().getX(),
                    (int) line.calculateCentroid().getY(),
                    HersheyFont.ROMAN_DUPLEX,
                    20, RGBColour.BLUE);
        }

        frame.drawText(lines.size() + "", (int) center.getX(),
                (int) center.getY(),
                HersheyFont.ROMAN_DUPLEX,
                40,
                RGBColour.RED);

    }
    private static double getAngleInDegree(Line2d line) {
        return line.calculateHorizontalAngle() * 180 / PI;
    }

    static double calcAngleDiffInDegree(double a1, double a2) {
        double gap = 0;
        if (a1 * a2 >= 0)
            return Math.abs(a1 - a2);

        if (a1 < 0) gap = a2 - a1;
        else gap = a1 - a2;

        if (gap >= 90)
            gap = 180 - gap;

        return gap;
    }

    private static void sortByYAxis(Line2d line) {
        if (line.begin.getY() > line.end.getY()) {
            Point2d tmp = line.begin;
            line.setBeginPoint(line.end);
            line.setEndPoint(tmp);
        }
    }

    private static void sortByXAxis(Line2d line) {
        if (line.begin.getX() > line.end.getX()) {
            Point2d tmp = line.begin;
            line.setBeginPoint(line.end);
            line.setEndPoint(tmp);
        }
    }


}
