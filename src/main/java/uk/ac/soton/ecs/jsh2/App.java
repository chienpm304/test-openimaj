package uk.ac.soton.ecs.jsh2;

import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.contour.Contour;
import org.openimaj.image.contour.SuzukiContourProcessor;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.point.Point2dImpl;
import org.openimaj.math.geometry.shape.Polygon;
import org.openimaj.math.geometry.shape.Rectangle;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * OpenIMAJ Hello world!
 *
 */
public class App {

    public static void main( String[] args ) throws IOException {
        File folder = new File("/home/cpu11427/chienpm/WhitePaper/test-threshold/input/dataset_03/original");
        if(folder.exists() && folder.isDirectory())
            for (final File file : folder.listFiles()) {
                if(file.isFile())
                    detectBox(file, folder);
            }


        File input_folder = new File("/home/cpu11427/chienpm/WhitePaper/test-threshold/input/dataset_03/original/rsz_t2.jpg");
        File out_folder = new File("/home/cpu11427/chienpm/WhitePaper/test-threshold/input/dataset_03/original");
        detectBox(input_folder, out_folder);

    }

    private static void detectBox(File fin, File folder) throws IOException {
        System.out.println("Processing "+fin.getName());

        MBFImage frame = ImageUtilities.readMBF(fin);



        //convert to hsv
        FImage s = MyHelper.RGB_TO_S_CHANNEL(frame);

        FImage grey = MyHelper.thresholdInv(s, 0.07133f);

        Contour contour = SuzukiContourProcessor.findContours(grey);

        Rectangle max = new Rectangle(0,0, 1, 1);
        Polygon fit = contour.clone();

        for(int i = 0; i < contour.children.size(); i++){

            Contour cnt = contour.children.get(i);

            Rectangle box = cnt.calculateRegularBoundingBox();
            double area = box.calculateArea();
            if(area > max.calculateArea() && area>1 && area < contour.calculateRegularBoundingBox().calculateArea()) {
                max = box.clone();
                fit = cnt.clone();

            }
        }

        System.out.println(max.calculateArea()+"| "+max.toString());

        List<Point2d> bound = MyHelper.findBounding(fit.points);

        frame.drawPoints(fit.points, RGBColour.GREEN, 4);

//        frame.drawShape(max, 8, RGBColour.BLUE);

        frame.drawLine(bound.get(0), bound.get(1), 10, RGBColour.BLUE);
        frame.drawLine(bound.get(1), bound.get(2), 10, RGBColour.BLUE);
        frame.drawLine(bound.get(2), bound.get(3), 10, RGBColour.BLUE);
        frame.drawLine(bound.get(3), bound.get(0), 10, RGBColour.BLUE);

        frame.drawPoints(bound, RGBColour.RED, 10);


        System.out.println(bound);

//        DisplayUtilities.display(frame);
        ImageUtilities.write(frame, new File(folder.getAbsolutePath()+"/out/"+fin.getName()));
    }




}
