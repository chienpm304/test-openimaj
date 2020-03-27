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
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import static uk.ac.soton.ecs.jsh2.DetectorUtils.*;

/**
 * Package com.chienpm.library.detector
 * Created by @chienpm on 28/02/2020
 */

public class BoundDetector {

    private static int width;
    private static int height;
    private static float scaleFactor;
    static File mfout;
    static File mfin;

    public static Tetragram detectBound(File fin, File fout, boolean shuffle) throws IOException {
        System.out.println("processing " + fin.getName());
        mfout = fout;
        mfin = fin;
        MBFImage frame = ImageUtilities.readMBF(fin);

        initSizeAndResizeImage(frame);

        width = frame.getWidth();
        height = frame.getHeight();

        Point2dImpl center = new Point2dImpl(width / 2, height / 2);

        enhanceInputImage(frame);

        FImage edges = applyCannyDetector(frame.flatten());
        ImageUtilities.write(edges, new File(fout.getAbsolutePath() + "/edged/" + fin.getName()));


        List<Line2d> lines = getLinesUsingHoughTransformP(edges, shuffle);
        saveToFile(fin, fout, frame, center, lines, false, "/raw_detected/");
        System.out.println("Before merge: " + lines.size());

        sortLinesPointFieldByAngle(lines);

        removeNoiseLines(lines, Constants.MERGE_MAX_LINE_DISTANCE, Constants.MERGE_MAX_LINE_GAP, Constants.MIN_ANGLE);
        saveToFile(fin, fout, frame, center, lines, true, "/merged/");

//        Tetragram bound = findBounds2(lines);
        List<LineHolder> bounds = findBounds3(lines);
        if (!bounds.isEmpty()) {
            for (LineHolder lh : bounds) {
                Tetragram bound = lh.tetragram;
                if (bound != null) {
                    DrawUtils.drawBound(
                            frame,
                            center,
                            bound.toLineList(),
                            RGBColour.GREEN,
                            RGBColour.YELLOW,
                            RGBColour.BLUE,
                            RGBColour.RED,
                            true);
                }
            }
        }
        ImageUtilities.write(frame, new File(fout.getAbsolutePath() + "/out/" + fin.getName()));
        return null;
    }

    private static void saveToFile(File fin, File fout, MBFImage frame, Point2dImpl center, List<Line2d> lines, boolean isDrawAngle, String s) throws IOException {
        MBFImage merged2 = frame.clone();
        DrawUtils.drawLines(merged2, center, lines, RGBColour.GREEN, isDrawAngle);
        ImageUtilities.write(merged2, new File(fout.getAbsolutePath() + s + fin.getName()));
    }

    public static void initSizeAndResizeImage(MBFImage frame) throws IOException {
        if (frame.getWidth() > frame.getHeight())
            scaleFactor = Constants.STANDARD_WIDTH / (float) frame.getHeight();
        else
            scaleFactor = Constants.STANDARD_WIDTH / (float) frame.getWidth();

        scaleFactor = scaleFactor > 1 ? 1.0f : scaleFactor;

        if (scaleFactor != 1.0f) {
            frame.processInplace(new ResizeProcessor(scaleFactor));
            ImageUtilities.write(frame, new File(mfout.getAbsolutePath() + "/new/" + mfin.getName()));
        }

        width = frame.getWidth();
        height = frame.getHeight();
    }

    private static void removeNoiseLines(List<Line2d> lines, int maxLineDistance, int maxLineGap, int minAngle) {
        removeLineByMergingX(lines, maxLineDistance, maxLineGap, minAngle);
        removeLineByMergingY(lines, maxLineDistance, maxLineGap, minAngle);
        removeLineWithShorterThanThreshold(lines, Constants.REMOVE_AFTER_MERGE_THRESHOLD);
    }

    private static void sortLinesPointFieldByAngle(List<Line2d> lines) {
        for (Line2d l : lines) {
            if (getUnsignedAngleInDegree(l) < 45) {
                sortByXAxis(l);
            } else {
                sortByYAxis(l);
            }
        }
    }

    private static void enhanceInputImage(MBFImage frame) {
//        frame = Transforms.RGB_TO_HSV(frame);
//        float[][] H = frame.getBand(1).pixels;
//
//        for(int y = 0; y < frame.getBand(1).height; ++y) {
//            for(int x = 0; x < frame.getBand(1).width; ++x) {
//                H[y][x] = Math.min(0.5f + H[y][x], 1.0f);
//            }
//        }
        applyGammaCorrection(frame);
//        applyGaussianConvolution(frame);
    }

    private static void applyGaussianConvolution(MBFImage frame) {
        FGaussianConvolve gauss = new FGaussianConvolve(Constants.GAUSSIAN_BLUR_SIGMA);
        for (int i = 0; i < frame.numBands(); i++) {
//            gauss.processImage(frame.getBand(i));
            frame.getBand(i).processInplace(gauss);
        }
    }

    private static void applyGammaCorrection(MBFImage frame) {
        GammaCorrection gc = new GammaCorrection(Constants.GAMMA_CORRECTION_VALUE);
        for (int i = 0; i < frame.numBands(); i++) {
//            gc.processImage(frame.getBand(i));
            frame.getBand(i).processInplace(gc);
        }
    }

    private static FImage applyCannyDetector(FImage grey) {
        CannyEdgeDetector canny =
                new CannyEdgeDetector(
                        Constants.CANNY_LOW_THRESH,
                        Constants.CANNY_HIGH_THRESH,
                        Constants.CANNY_SIGMA);

        canny.processImage(grey);
        return grey;
    }

    private static List<Line2d> getLinesUsingHoughTransformP(FImage image, boolean shuffle) {
        List<Point2d> points = new ArrayList<>();
        for (int i = 0; i < image.width; i++)
            for (int j = 0; j < image.height; j++) {
                if (image.getPixel(i, j) > 0.5f)
                    points.add(new Point2dImpl(i, j));
            }

        if (shuffle)
            Collections.shuffle(points);

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

    private static List<LineHolder> findBounds3(List<Line2d> lines) {

        ArrayList<LineHolder> verticals = new ArrayList<>();
        ArrayList<LineHolder> horizontals = new ArrayList<>();
        // paring
        ArrayList<LineHolder> res = new ArrayList<>();

        LineHolder holder;

        Line2d base1, base2;
        double b1, b2;
        int angle_step = 20;

        // find vertical line (left, right)
        while (verticals.isEmpty() && angle_step <= Constants.ANGLE_DIFF_THRESHOLD) {

            for (int i = 0; i < lines.size() - 1; i++) {

                base1 = lines.get(i);
                b1 = getUnsignedAngleInDegree(base1);

                for (int j = i + 1; j < lines.size(); j++) {
                    base2 = lines.get(j);
                    b2 = getUnsignedAngleInDegree(base2);

                    if (Math.abs(b1 - b2) <= angle_step) {
                        holder = new LineHolder();
                        if (b1 > 45) { // left, right
                            if (base1.calculateCentroid().getX() < base2.calculateCentroid().getX()) {
                                holder.left = base1;
                                holder.right = base2;
                            } else {
                                holder.left = base2;
                                holder.right = base1;
                            }
                            verticals.add(holder);
                        } else { // top, bottom
                            if (base1.calculateCentroid().getY() < base2.calculateCentroid().getY()) {
                                holder.top = base1;
                                holder.bottom = base2;
                            } else {
                                holder.top = base2;
                                holder.bottom = base1;
                            }
                            horizontals.add(holder);
                        }
                    }
                }
            }
            angle_step += 5;
        }

        int select_line_const = 30;
        double a1, a2, gapAngle;

            for (LineHolder v : verticals) {
                for (LineHolder h : horizontals) {
                    a1 = DetectorUtils.getSignedAngleInDegree(v.left);
                    a2 = DetectorUtils.getSignedAngleInDegree(h.top);

                    gapAngle = DetectorUtils.calcAngleDiffInDegree(a1, a2);

                    if (gapAngle >= 90 - angle_step) {
                        LineHolder lh = new LineHolder(v.left, h.top, v.right, h.bottom);
                        System.out.println("Checking " + lh.toString());
                        lh.compute(width, height);
                        if (lh.isSatisfyingShape()) {
                            res.add(lh);
                            System.out.println("Accepted");
                        }
                    }
                }
            }
        System.out.println("Detected " + res.size() + " bounds");

        if (res.isEmpty())
            return res;

        Collections.sort(res, new Comparator<LineHolder>() {
            @Override
            public int compare(LineHolder o1, LineHolder o2) {
                return Double.compare(o2.area, o1.area); //reserved
            }
        });

        /**
         * Nếu có rất nhiều bounds thì chọn bound lớn nhất, it bound thì duyệt chọn bound nhỏ nhưng ngon
         * hoặc là chọn 2 thằng lớn nhất ra so sánh gap
         *
         */
        int idx = 0;
        double min = res.get(0).gap;
//        for (int i = 1; i < 5 && i < res.size(); i++) {
//            if (res.get(i).gap < min
//                    && res.get(i).area > res.get(idx).area * 0.8) {
//                min = res.get(i).gap;
//                idx = i;
//            }
//        }
        return res;//.get(idx).tetragram;
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

    private static void removeLineByMergingX(List<Line2d> lines, int maxLineDistance, int maxLineGap, int minAngle) {
        Collections.sort(lines, new Comparator<Line2d>() {
            @Override
            public int compare(Line2d o1, Line2d o2) {
                return Double.compare(o2.calculateLength(), o1.calculateLength()); //reserved
            }
        });

        Line2d l1, l2;
        int beforeMerge, afterMerge;
        List<Line2d> tmpLines = new ArrayList<>();
//        System.out.print("Merging: ");
        for (int i = 0; i < lines.size(); i++) {
            tmpLines.clear();
            l1 = lines.get(i);
            if (getUnsignedAngleInDegree(l1) >= 45)
                continue;
            tmpLines.add(l1);
            for (int j = i + 1; j < lines.size(); j++) {
                l2 = lines.get(j);
                if (checkIfMayOnTheSameLine(tmpLines, l2, maxLineDistance, minAngle)) {
                    tmpLines.add(l2);
                }
            }
            if (tmpLines.size() > 1) {
                lines.removeAll(tmpLines);
                System.out.print(" " + tmpLines.size());
                beforeMerge = tmpLines.size();
                mergeLines2(tmpLines, maxLineGap);
                afterMerge = tmpLines.size();
                lines.addAll(tmpLines);
                if (afterMerge < beforeMerge)
                    i--;
            }
        }
        System.out.println();
    }

    private static void removeLineByMergingY(List<Line2d> lines, int maxLineDistance, int maxLineGap, int minAngle) {
        lines.sort(new Comparator<Line2d>() {
            @Override
            public int compare(Line2d o1, Line2d o2) {
                return Double.compare(o2.calculateLength(), o1.calculateLength()); //reserved
            }
        });

        Line2d l1, l2;
        int beforeMerge, afterMerge;
        List<Line2d> tmpLines = new ArrayList<>();
//        System.out.print("Merging: ");
        for (int i = 0; i < lines.size(); i++) {
            tmpLines.clear();
            l1 = lines.get(i);
            if (getUnsignedAngleInDegree(l1) < 45)
                continue;
            tmpLines.add(l1);
            for (int j = i + 1; j < lines.size(); j++) {
                l2 = lines.get(j);
                if (checkIfMayOnTheSameLine(tmpLines, l2, maxLineDistance, minAngle)) {
                    tmpLines.add(l2);
                }
            }
            if (tmpLines.size() > 1) {
                lines.removeAll(tmpLines);
                System.out.print(" " + tmpLines.size());
                beforeMerge = tmpLines.size();
                mergeLines2(tmpLines, maxLineGap);
                afterMerge = tmpLines.size();
                lines.addAll(tmpLines);
                if (afterMerge < beforeMerge)
                    i--;
            }
        }
        System.out.println();
    }

    private static void mergeLines2(List<Line2d> lines, int maxLineGap) {
        boolean mergeX = getUnsignedAngleInDegree(lines.get(0)) < 45;

        Collections.sort(lines, new Comparator<Line2d>() {
            @Override
            public int compare(Line2d o1, Line2d o2) {
                return Double.compare(o2.calculateLength(), o1.calculateLength()); //reserved
            }
        });

        Line2d l1, l2;
        int idx;
        for (int i = 0; i < lines.size(); i++) {
            l1 = lines.get(i);
            while ((idx = findClosestLineByLineGap(lines, l1, maxLineGap)) != -1) {
                if (i >= lines.size() || idx >= lines.size())
                    break;

                l2 = lines.get(idx);

                if (l1.calculateLength() < l2.calculateLength()) {
                    Collections.swap(lines, i, idx);
                    l1 = lines.get(i); // keep this
                    l2 = lines.get(idx);
                }

                if (mergeX) mergeX(l1, l2);
                else mergeY(l1, l2);

                lines.remove(idx);

            }
        }
    }

    private static int findClosestLineByLineGap(List<Line2d> lines, Line2d line, int maxLineGap) {
        double minGap = maxLineGap + 1;
        int idx = -1;
        Line2d l1;
        for (int i = 0; i < lines.size(); i++) {
            l1 = lines.get(i);
            if (l1 != line) {
                double gap = calculateLineGap(l1, line);
                if (gap < minGap && gap < maxLineGap) {
                    minGap = gap;
                    idx = i;
                }
            }
        }
        return idx;
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

    private static void mergeY(Line2d keep, Line2d remove) {
        // begin.x < end.x
//        sortByYAxis(keep);
//        sortByYAxis(remove);
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
//        sortByXAxis(keep);
//        sortByXAxis(remove);
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


}
