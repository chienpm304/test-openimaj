package uk.ac.soton.ecs.jsh2;

public class Constants {
    public static final float STANDARD_WIDTH = 720;
    public static final int S_CHANNEL_ID = 1;
    public static final String WINDOW_DIR = "D:/detect/input/AZdoc/new";
    public static final String WINDOW_OUT_DIR = "D:/detect/input/AZdoc";
    public static final String LINUX_DIR_IN = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/AZdoc/new";
    public static String LINUX_DIR_OUT = "/home/cpu11427/chienpm/WhitePaper/test-threshold/input/AZdoc/";


    public static final float GAUSSIAN_BLUR_SIGMA = 2f;// min distance of 2 lines calculated by line.distanceToPoint

    public static final int MIN_ANGLE = 8;

    public static final int MERGE_MAX_LINE_DISTANCE = 30;// min gap of 2 line's points = min(l1.begin -> l2.begin, l1.begin->l2.end, l1.end->l2.begin, l1.end -> l2.end)
    public static final int MERGE_MAX_LINE_GAP = 100; //50is ok
    public static final int MERGE_GAP_PER_LINE = 10;
    public static final int BOUNDING_GAP_REMOVAL = 3;

    //0.05 - 0.1 is ok (0.01 - 0.05)
    public static float CANNY_LOW_THRESH = 0.03f;
    public static float CANNY_HIGH_THRESH = 0.1f;
    public static float CANNY_SIGMA = 3f;

    public static final double HOUGH_LINE_RHO = 1;
    public static final double HOUGH_LINE_THETA = Math.PI / 180d;;

    public static int HOUGH_LINE_MAX_LINE_GAP = 10;
    public static int HOUGH_LINE_THRESHOLD = 20;
    public static int HOUGH_MIN_LINE_LENGTH = 20;
    public static int HOUGH_MAX_NUM_LINE = 1000;

    public static final double GAMMA = 2.2d;

    public static final int REMOVE_AFTER_MERGE_THRESHOLD = 100;

}