package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.FImage;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.ColourSpace;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class MyHelper {
    public String root = "D:/detect";
    public static Random rnd = new Random();

    public MyHelper() {
    }

    public static Float[] getRandomColor() {
        return RGBColour.RGB(rnd.nextInt(256), rnd.nextInt(256), rnd.nextInt(256));
    }



    public static FImage thresholdInv(FImage img, final Float thresh) {
        final float fthresh = thresh;
        for (int r = 0; r < img.height; r++) {
            for (int c = 0; c < img.width; c++) {
                if (img.pixels[r][c] <= fthresh)
                    img.pixels[r][c] = 1;
                else
                    img.pixels[r][c] = 0;
            }
        }
        return img;
    }


    public static List<Point2d> findBounding(List<Point2d> points) {
        float x_min = Float.MAX_VALUE;
        float x_max = -1;
        float y_min = Float.MAX_VALUE;
        float y_max = -1;
        float x,y;
        for(Point2d p: points){
            x = p.getX();
            y = p.getY();
            if(x > x_max) x_max = x;
            if(x < x_min) x_min = x;
            if(y < y_min) y_min = y;
            if(y > y_max) y_max = y;

        }

        Point2dImpl
                tl = new Point2dImpl(x_min, y_max),
                tr =new Point2dImpl(x_max, y_max),
                bl = new Point2dImpl(x_min,y_min),
                br = new Point2dImpl(x_max, y_min);

        Point2dImpl r_tl = br.clone(), r_tr = bl.clone(), r_bl = tr.clone(), r_br = tl.clone();


        float width = x_max - x_min, height = y_max - y_min;

        for(Point2d p: points){
            x = p.getX();
            y = p.getY();
            //top left
            if(x < width/2 && y>height/2){
                if(MyHelper.distance(p, tl) < MyHelper.distance(r_tl, tl)){
                    r_tl.x = p.getX();
                    r_tl.y = p.getY();
                }
            }
            else if(x > width/2 && y > height/2) {// top right
                if(MyHelper.distance(p, tr) < MyHelper.distance(r_tr, tr)){
                    r_tr.x = p.getX();
                    r_tr.y = p.getY();
                }
            }
            else if(x < width/2 && y < height/2) {// bottom left
                if(MyHelper.distance(p, bl) < MyHelper.distance(r_bl, bl)){
                    r_bl.x = p.getX();
                    r_bl.y = p.getY();
                }
            }
            else {//bottom right
                if(MyHelper.distance(p, br) < MyHelper.distance(r_br, br)){
                    r_br.x = p.getX();
                    r_br.y = p.getY();
                }
            }

        }

        // small add]just
        if(MyHelper.distance(r_bl, bl) > width/3) r_bl = bl;
        if(MyHelper.distance(r_tl, tl) > width/3) r_tl = tl;
        if(MyHelper.distance(r_br, br) > width/3) r_br = br;
        if(MyHelper.distance(r_tr, tr) > width/3) r_tr = tr;

        List<Point2d> res = new ArrayList<>();
        res.add(r_tl);
        res.add(r_tr);
        res.add(r_br);
        res.add(r_bl);
        return res;
    }

    public static double distance(Point2d p1, Point2d p2) {

        float x1 = p1.getX(), x2 = p2.getX(), y1 = p1.getY(), y2 = p2.getY();

        return Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }


    public static FImage RGB_TO_S_CHANNEL(final MBFImage in) {
        if (in.colourSpace != ColourSpace.RGB && in.colourSpace != ColourSpace.RGBA)
            throw new IllegalArgumentException("RGB or RGBA colourspace is required");

        final int width = in.getWidth();
        final int height = in.getHeight();

        final FImage out = new FImage(width, height);

        final float[][] R = in.getBand(0).pixels;
        final float[][] G = in.getBand(1).pixels;
        final float[][] B = in.getBand(2).pixels;

        final float[][] S = out.pixels;

        final float[] pIn = new float[3];
        final float[] pOut = new float[3];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {

                pIn[0] = R[y][x];
                pIn[1] = G[y][x];
                pIn[2] = B[y][x];

                RGB_TO_HSV(pIn, pOut);

                S[y][x] = pOut[1];
            }
        }
        return out;
    }

    public static float[] RGB_TO_HSV(final float[] rgb, final float[] hsv) {
        float H, S, V;

        final float R = rgb[0];
        final float G = rgb[1];
        final float B = rgb[2];

        // Blue Is the dominant color
        if ((B > G) && (B > R)) {
            // Value is set as the dominant color
            V = B;
            if (V != 0) {
                float min;
                if (R > G)
                    min = G;
                else
                    min = R;

                // Delta is the difference between the most dominant
                // color
                // and the least dominant color. This will be used to
                // compute saturation.
                final float delta = V - min;
                if (delta != 0) {
                    S = (delta / V);
                    H = 4 + (R - G) / delta;
                } else {
                    S = 0;
                    H = 4 + (R - G);
                }

                // Hue is just the difference between the two least
                // dominant
                // colors offset by the dominant color. That is, here 4
                // puts
                // hue in the blue range. Then red and green just tug it
                // one
                // way or the other. Notice if red and green are equal,
                // hue
                // will stick squarely on blue
                H *= 60;
                if (H < 0)
                    H += 360;

                H /= 360;
            } else {
                S = 0;
                H = 0;
            }
        }
        // Green is the dominant color
        else if (G > R) {
            V = G;
            if (V != 0) {
                float min;
                if (R > B)
                    min = B;
                else
                    min = R;

                final float delta = V - min;

                if (delta != 0) {
                    S = (delta / V);
                    H = 2 + (B - R) / delta;
                } else {
                    S = 0;
                    H = 2 + (B - R);
                }
                H *= 60;
                if (H < 0)
                    H += 360;

                H /= 360;
            } else {
                S = 0;
                H = 0;
            }
        }
        // Red is the dominant color
        else {
            V = R;
            if (V != 0) {
                float min;
                if (G > B)
                    min = B;
                else
                    min = G;

                final float delta = V - min;
                if (delta != 0) {
                    S = (delta / V);
                    H = (G - B) / delta;
                } else {
                    S = 0;
                    H = (G - B);
                }
                H *= 60;

                if (H < 0)
                    H += 360;
                H /= 360;
            } else {
                S = 0;
                H = 0;
            }
        }

        hsv[0] = H;
        hsv[1] = S;
        hsv[2] = V;
        return hsv;
    }
}