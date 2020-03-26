package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.typography.hershey.HersheyFont;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2d;

import java.util.List;
import java.util.Random;

public class DrawUtils {
    static void drawLineHolder(MBFImage frame, List<Line2d> lines, int code, Float[] baseColor, Float[] fitColor) {
        Float[] lineColor;
        for (Line2d line : lines) {
            if (lines.indexOf(line) % 2 == 0)
                lineColor = baseColor;
            else
                lineColor = fitColor;

            frame.drawLine(line, 2, lineColor);
            frame.drawPoint(line.begin, RGBColour.RED, 5);
            frame.drawPoint(line.end, RGBColour.YELLOW, 5);

            frame.drawText(
                    code + "*",
                    (int) line.calculateCentroid().getX(),
                    (int) line.calculateCentroid().getY(),
                    HersheyFont.ROMAN_DUPLEX,
                    20, RGBColour.BLUE);
        }
    }

    static void drawBound(MBFImage frame, Point2d center, List<Line2d> lines, Float[] topColor, Float[] rightColor,
                          Float[] bottomColor, Float[] leftColor, boolean drawAngle) {
        Float[] lineColor;
        for (Line2d line : lines) {
            switch (lines.indexOf(line)){
                case 0: lineColor = topColor; break;
                case 1: lineColor = rightColor; break;
                case 2: lineColor = bottomColor; break;
                default: lineColor = leftColor;
            }

            frame.drawLine(line, 2, lineColor);
            frame.drawPoint(line.begin, RGBColour.RED, 5);
            frame.drawPoint(line.end, RGBColour.YELLOW, 5);
            if (drawAngle)
                frame.drawText(
                        Math.round(DetectorUtils.getUnsignedAngleInDegree(line)) + "*",
                        (int) line.calculateCentroid().getX(),
                        (int) line.calculateCentroid().getY(),
                        HersheyFont.ROMAN_DUPLEX,
                        20, RGBColour.BLUE);
        }
    }

    protected static void drawLines(MBFImage frame, Point2d center, List<Line2d> lines, Float[] lineColor, boolean drawAngle) {
        System.out.println(lines.size() + " lines");

        Float[] orgColor = lineColor;

        for (Line2d line : lines) {
            frame.drawLine(line, 2, lineColor);
            frame.drawPoint(line.begin, RGBColour.RED, 5);
            frame.drawPoint(line.end, RGBColour.YELLOW, 5);
            if (drawAngle)
                frame.drawText(
                        Math.round(DetectorUtils.getUnsignedAngleInDegree(line)) + "*",
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
    static Random rnd = new Random();

    public static Float[] getRandomColor() {
        return RGBColour.RGB(rnd.nextInt(256), rnd.nextInt(256), rnd.nextInt(256));
    }
}