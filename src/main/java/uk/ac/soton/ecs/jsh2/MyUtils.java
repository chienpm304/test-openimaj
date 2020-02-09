package uk.ac.soton.ecs.jsh2;

import Jama.Matrix;
import org.openimaj.image.MBFImage;
import org.openimaj.math.geometry.point.Point2d;

public class MyUtils {
    private static Tetragram findDestinationRectangle(Tetragram bound) {
        Point2d tl = bound.getTopLeft();
        Point2d tr = bound.getTopRight();
        Point2d br = bound.getBottomRight();
        Point2d bl = bound.getBottomLeft();

        int newWidth = (int) Math.max(App.distance(br, bl), App.distance(tr, tl));
        int newHeight = (int) Math.max(App.distance(tr, br), App.distance(tl, bl));

        float[][] rect = new float[][]{
                {0, 0},
                {newWidth - 1, 0},
                {newWidth - 1, newHeight - 1},
                {0, newHeight - 1}
        };

        return new Tetragram(rect);
    }


    private static MBFImage transformImage(MBFImage frame, Matrix transformMatrix, int newWidth, int newHeight) {
        MBFImage img = new MBFImage(newWidth, newHeight);
//        img.processInplace(new ;
        return null;
    }

    private static Matrix getTransformMatrix(Tetragram originTetra, Tetragram destTetra) {
        float[][] src = originTetra.toArray();
        float[][] dst = destTetra.toArray();

        double A[][] = new double[8][8];
        // A x h = b -> find h = (A)^-1 x b
        for(int i =0; i < 4; i++){
            A[2*i][0] = src[i][0];
            A[2*i][1] = src[i][1];
            A[2*i][2] = 1;
            A[2*i][3] = 0;
            A[2*i][4] = 0;
            A[2*i][5] = 0;
            A[2*i][6] = -src[i][0]*dst[i][0];
            A[2*i][7] = -dst[i][0]*src[i][1];

            A[2*i + 1][0] = 0;
            A[2*i + 1][1] = 0;
            A[2*i + 1][2] = 0;
            A[2*i + 1][3] = src[i][0];
            A[2*i + 1][4] = src[i][1];
            A[2*i + 1][5] = 1;
            A[2*i + 1][6] = -src[i][0]*dst[i][1];
            A[2*i + 1][7] = -src[i][1]*dst[i][1];
        }

        // build A (8x8)
        Matrix matrixA = new Matrix(A);

        // build b (8x1)
        double []b = new double[8];
        for(int i = 0; i<4; i++){
            b[2*i] = dst[i][0];
            b[2*i+1] = dst[i][1];
        }

        Matrix matrixB = new Matrix(b, 8);

        Matrix trsf = matrixA.inverse().times(matrixB);

        // reshape to 3x3
        double[] flat = trsf.getColumnPackedCopy();
        Matrix trf3x3 = new Matrix(3,3);
        for(int i = 0; i<8; i++)
            trf3x3.set(i/3, i%3, flat[i]);
        trf3x3.set(2,2, 1);

        return trf3x3;
    }

}