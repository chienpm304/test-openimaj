package uk.ac.soton.ecs.jsh2;

/**
 * Package com.chienpm.library.detector
 * Created by @chienpm on 20/02/2020
 */

public class Constants {

    public static final int ANGLE_DIFF_THRESHOLD = 60;
    private static final float MAGIC_RATIO = (360f / 720f);

    static final float STANDARD_WIDTH = 720 * MAGIC_RATIO;

    static final int MIN_ANGLE = 8;

    static final int MERGE_MAX_LINE_DISTANCE = (int) (30 * MAGIC_RATIO);// min gap of 2 line's points = min(l1.begin -> l2.begin, l1.begin->l2.end, l1.end->l2.begin, l1.end -> l2.end)
    static final int MERGE_MAX_LINE_GAP = (int) (100 * MAGIC_RATIO); //50is ok

    static final int BOUNDING_GAP_REMOVAL = 3;

    static final double HOUGH_LINE_RHO = 1;
    static final double HOUGH_LINE_THETA = Math.PI / 180d;

    static final int HOUGH_LINE_MAX_LINE_GAP = 10;
    static final int HOUGH_LINE_THRESHOLD = (int) (20 * MAGIC_RATIO);
    static final int HOUGH_MIN_LINE_LENGTH = (int) (20 * MAGIC_RATIO);
    static final int HOUGH_MAX_NUM_LINE = 1000;

    static final float CANNY_LOW_THRESH = 0.05f;
    static final float CANNY_HIGH_THRESH = 0.25f;
    static final float CANNY_SIGMA = 5f;

    static final double GAMMA_CORRECTION_VALUE = 2.2d;
    static final float GAUSSIAN_BLUR_SIGMA = 2.5f;

    static final int REMOVE_AFTER_MERGE_THRESHOLD = (int) (100 * MAGIC_RATIO);

    /**
     *  Hints
     * Exposure Adjust the amount of light in the image.
     * Shadows Only adjust the brightness of the shadows, or darkest areas of the image.
     * Highlight Only adjust the brightness of the highlights, or brightest areas of the image.
     * Contrast Adjust the degree of difference between shadows and highlights.
     * Structure Increase the amount of details in the image. Structure uses a unique algorithm to bring out the texture of objects throughout the photo, without affecting the edges of the objects.
     * Temperature Adjust the balance between cool blue tones and warm yellow tones in the image.
     * Tint Adjust the balance between green and magenta tones in the image.
     */

}