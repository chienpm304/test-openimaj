package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.FImage;
import org.openimaj.image.pixel.FValuePixel;
import org.openimaj.math.geometry.line.Line2d;
import org.openimaj.math.geometry.point.Point2dImpl;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class TestHough {

    private FImage accum;
    private int numberOfSegments;
    private FImage iteratorAccum;
    private FValuePixel iteratorCurrentPix;
    private float onValue;

    public TestHough() {
        this(360, 0.0F);
    }

    public TestHough(float onValue) {
        this(360, onValue);
    }

    public TestHough(int nSegments, float onValue) {
        this.accum = null;
        this.numberOfSegments = 360;
        this.iteratorAccum = null;
        this.iteratorCurrentPix = null;
        this.setNumberOfSegments(nSegments);
        this.onValue = onValue;
    }

    public void analyseImage(FImage image) {
        int amax = (int)Math.round(Math.sqrt((double)(image.getHeight() * image.getHeight() + image.getWidth() * image.getWidth())));
        if (this.accum != null && this.accum.height == amax && this.accum.width == this.getNumberOfSegments()) {
            this.accum.zero();
        } else {
            this.accum = new FImage(this.getNumberOfSegments(), amax);
        }

        for(int y = 0; y < image.getHeight(); ++y) {
            for(int x = 0; x < image.getWidth(); ++x) {
                if (image.getPixel(x, y) == this.onValue) {
                    for(int m = 0; m < this.getNumberOfSegments(); ++m) {
                        double mm = (double)m / (double)this.getNumberOfSegments() * 6.283185307179586D;
                        int a = (int)Math.round((double)x * Math.cos(mm) + (double)y * Math.sin(mm));
                        if (a < amax && a >= 0) {
                            this.accum.pixels[a][m]++;
                        }
                    }
                }
            }
        }

    }

    public Line2d getBestLine(FImage accumulatorSpace, int offset) {
        FValuePixel p = accumulatorSpace.maxPixel();
        int theta = p.x + offset;
        int dist = p.y;
        return this.getLineFromParams(theta, dist, -2000, 2000);
    }


    public Line2d getLineFromParams(int theta, int dist, int x1, int x2) {
        if (theta == 0) {
            return new Line2d(new Point2dImpl((float)dist, -2000.0F), new Point2dImpl((float)dist, 2000.0F));
        } else {
            double t = (double)theta * (360.0D / (double)this.getNumberOfSegments()) * 3.141592653589793D / 180.0D;
            return new Line2d(new Point2dImpl((float)x1, (float)((double)x1 * (-Math.cos(t) / Math.sin(t)) + (double)dist / Math.sin(t))), new Point2dImpl((float)x2, (float)((double)x2 * (-Math.cos(t) / Math.sin(t)) + (double)dist / Math.sin(t))));
        }
    }


    public boolean hasNext() {
        return this.iteratorAccum.maxPixel().value > 0.0F;
    }

    public Line2d next() {
        this.iteratorCurrentPix = this.iteratorAccum.maxPixel();
        Line2d l = this.getBestLine(this.iteratorAccum, 0);
        this.iteratorAccum.setPixel(this.iteratorCurrentPix.x, this.iteratorCurrentPix.y, 0.0F);
        return l;
    }


    public void setNumberOfSegments(int numberOfSegments) {
        this.numberOfSegments = numberOfSegments;
    }

    public int getNumberOfSegments() {
        return this.numberOfSegments;
    }
}
