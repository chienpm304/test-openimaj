package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.processing.algorithm.GammaCorrection;
import org.openimaj.image.processing.convolution.FGaussianConvolve;
import org.openimaj.image.processing.edges.CannyEdgeDetector;
import org.openimaj.image.processing.resize.ResizeProcessor;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import static uk.ac.soton.ecs.jsh2.App.*;

public class DetectByCanny {

    public static float _CANNY_LOW_THRESH = 0.05f;
    public static float _CANNY_HIGH_THRESH = 0.25f;
    public static float _CANNY_SIGMA = 5f;

    public static final double _HOUGH_LINE_RHO = 1;
    public static final double _HOUGH_LINE_THETA = Math.PI / 180d;;

    public static int _HOUGH_LINE_MAX_LINE_GAP = 30;
    public static int _HOUGH_LINE_THRESHOLD = 30;
    public static int _HOUGH_MIN_LINE_LENGTH = 20;
    public static int _HOUGH_MAX_NUM_LINE = 1000;

    public static final double _GAMMA = 2.2d;
    public static final float _GAUSSIAN_BLUR_SIGMA = 2.5f;


    private static float scaleFactor;
    public static int width;
    public static int height;

    public static Tetragram detectWithCanny(File fin, File fout) throws IOException {
        MBFImage frame = ImageUtilities.readMBF(fin);
        DetectByCanny.initSizeAndResizeImage(frame);
//        if(frame.getWidth() > frame.getHeight())
//            scaleFactor = Constants.STANDARD_WIDTH / (float) frame.getHeight();
//        else
//            scaleFactor = Constants.STANDARD_WIDTH / (float) frame.getWidth();
//
//        scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;
//
//        if(scaleFactor!=1.0f) {
//            frame.processInplace(new ResizeProcessor(scaleFactor));
//        }
//
//        width = frame.getWidth();
//        height = frame.getHeight();

        System.out.println("processing: " + fin.getName() + "...");

        DetectByCanny.enhanceInputImage(frame);
//        GammaCorrection gc = new GammaCorrection(Constants.GAMMA);
//        for (int i = 0; i < frame.numBands(); i++) {
//            frame.getBand(i).processInplace(gc);
//        }

        Point2dImpl center = new Point2dImpl(width / 2, height / 2);

        FImage edges = DetectByCanny.applyCannyDetector(frame.flattenMax());

        ImageUtilities.write(edges, new File(fout.getAbsolutePath() + "/edged/" + fin.getName()));

        List<Line2d> lines = getLinesUsingHoughTransformP(edges);

        MBFImage tmp = frame.clone();
        drawLines(tmp,center, lines, RGBColour.GREEN, false);
        ImageUtilities.write(tmp, new File(fout.getAbsolutePath() + "/raw_detected/" + fin.getName()));

        System.out.println("Before merge: " + lines.size());

        App.removeNoiseLines(lines, Constants.MERGE_MAX_LINE_DISTANCE, Constants.MERGE_MAX_LINE_GAP);

        System.out.println("After merge: " + lines.size());

        MBFImage merged = frame.clone();
        drawLines(merged, center, lines, RGBColour.GREEN, true);
        ImageUtilities.write(merged, new File(fout.getAbsolutePath() + "/merged/" + fin.getName()));

        Tetragram bound = findBounds2(lines);

        if (bound!=null) {
            drawBound(frame, center, bound.toLineList(), RGBColour.GREEN, RGBColour.YELLOW);
        }

        ImageUtilities.write(frame, new File(fout.getAbsolutePath() + "/out/" + fin.getName()));
        return null;
    }
    static void removeNoiseLines(List<Line2d> lines, int maxLineDistance, int maxLineGap) {

        removeLinesNearbyBounding(lines, Constants.BOUNDING_GAP_REMOVAL);

        DetectByCanny.removeLineByMergingX(lines, maxLineDistance, maxLineGap);
        DetectByCanny.removeLineByMergingY(lines, maxLineDistance, maxLineGap);

        removeLineWithShorterThanThreshold(lines, Constants.REMOVE_AFTER_MERGE_THRESHOLD);
    }


    public static void removeLineByMergingX(List<Line2d> lines, int maxLineDistance, int maxLineGap) {
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
    public static void removeLineByMergingY(List<Line2d> lines, int maxLineDistance, int maxLineGap) {
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


    public static boolean checkIfMayOnTheSameLine(List<Line2d> lines, Line2d line, int maxLineDistance) {
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
                        && Math.abs(getHorizontalAngleInDegree(l1) - getHorizontalAngleInDegree(l2)) <= Constants.MIN_ANGLE;
    }


    public static Line2d mergeLineNew(List<Line2d> lines) {
        if(lines.isEmpty()) return null;


        Line2d keep, remove;
        keep = lines.remove(0);
        boolean mergX = getHorizontalAngleInDegree(keep) < 45;
        System.out.print(mergX);

        while(!lines.isEmpty()){
            remove = lines.remove(0);
            if(keep.calculateLength() < remove.calculateLength()){
                swap2Lines(keep, remove);
            }
            if(mergX)
                mergeX(keep, remove);
            else
                mergeY(keep, remove);
        }
        return keep;
    }

    private static void swap2Lines(Line2d keep, Line2d remove) {
        Point2d bg = keep.begin;
        Point2d ed = keep.end;

        keep.begin = remove.begin;
        keep.end = remove.end;

        remove.begin = bg;
        remove.end = ed;
    }

    public static void enhanceInputImage(MBFImage frame) {
        // process gamma correction
        GammaCorrection gc = new GammaCorrection(_GAMMA);
        for (int i = 0; i < frame.numBands(); i++) {
            frame.getBand(i).processInplace(gc);

        }
        frame.processInplace(new FGaussianConvolve(_GAUSSIAN_BLUR_SIGMA));
    }

    public static void initSizeAndResizeImage(MBFImage frame) {
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
    }

    static FImage applyCannyDetector(FImage grey) {
        CannyEdgeDetector canny = new CannyEdgeDetector(_CANNY_LOW_THRESH, _CANNY_HIGH_THRESH, _CANNY_SIGMA);
        canny.processImage(grey);
        return grey;
    }

    public static List<Line2d> getLinesUsingHoughTransformP(FImage image) {
        List<Point2d> points = new ArrayList<>();
        for (int i = 0; i < image.width; i++)
            for (int j = 0; j < image.height; j++) {
                if (image.getPixel(i, j) > 0.5f)
                    points.add(new Point2dImpl(i, j));
            }
//        HoughLinesP ht = new HoughLinesP(
//                points,
//                width,
//                height,
//                _HOUGH_LINE_RHO,
//                _HOUGH_LINE_THETA,
//                _HOUGH_LINE_THRESHOLD,
//                _HOUGH_LINE_MAX_LINE_GAP,
//                _HOUGH_MIN_LINE_LENGTH,
//                _HOUGH_MAX_NUM_LINE);
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

}
