package uk.ac.soton.ecs.jsh2.old;

class LsdUtils {
    /**
     * Size of the table to store already computed inverse values.
     */
    static int TABSIZE = 100000;

    static double[] inv = new double[TABSIZE]; /*

    /** ln(10) */
    static double M_LN10 = 2.30258509299404568402;


    /**
     * PI
     */
    static double M_PI = 3.14159265358979323846;

    /**
     * Label for pixels with undefined gradient.
     */
    static double NOTDEF = -1024.0;

    /**
     * 3/2 pi
     */
    static double M_3_2_PI = 4.71238898038;

    /**
     * 2 pi
     */
    static double M_2__PI = 6.28318530718;

    static void error(String msg) {
        System.err.println("LSD Error: " + msg);
        throw new RuntimeException("major error");
    }

    static double angle_diff(double a, double b) {
        a -= b;
        while (a <= -M_PI)
            a += M_2__PI;
        while (a > M_PI)
            a -= M_2__PI;
        if (a < 0.0)
            a = -a;
        return a;
    }

    static boolean isAligned(int x, int y, ImageDouble angles, double theta,
                             double prec) {
        double a;

        /* check parameters */
        if (angles == null || angles.data == null)
            error("isaligned: invalid image 'angles'.");
        if (x < 0 || y < 0 || x >= (int) angles.width
                || y >= (int) angles.height)
            error("isaligned: (x,y) out of the image.");
        if (prec < 0.0)
            error("isaligned: 'prec' must be positive.");

        /* angle at pixel (x,y) */
        a = angles.data[x + y * angles.width];

        /*
         * pixels whose level-line angle is not defined are considered as
         * NON-aligned
         */
        if (a == NOTDEF)
            return false;

        /* it is assumed that 'theta' and 'a' are in the range [-pi,pi] */
        theta -= a;
        if (theta < 0.0)
            theta = -theta;
        if (theta > M_3_2_PI) {
            theta -= M_2__PI;
            if (theta < 0.0)
                theta = -theta;
        }

        return theta <= prec;
    }

    static double angle_diff_signed(double a, double b) {
        a -= b;
        while (a <= -M_PI)
            a += M_2__PI;
        while (a > M_PI)
            a -= M_2__PI;
        return a;
    }

    static double dist(double x1, double y1, double x2, double y2) {
        return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    }

    static double nfa(int n, int k, double p, double logNT) {

        double tolerance = 0.1; /* an error of 10% in the result is accepted */
        double log1term, term, bin_term, mult_term, bin_tail, err, p_term;
        int i;

        /* check parameters */
        if (n < 0 || k < 0 || k > n || p <= 0.0 || p >= 1.0)
            error("nfa: wrong n, k or p values.");

        /* trivial cases */
        if (n == 0 || k == 0)
            return -logNT;
        if (n == k)
            return -logNT - (double) n * Math.log10(p);

        /* probability term */
        p_term = p / (1.0 - p);

        log1term = log_gamma((double) n + 1.0) - log_gamma((double) k + 1.0)
                - log_gamma((double) (n - k) + 1.0) + (double) k * Math.log(p)
                + (double) (n - k) * Math.log(1.0 - p);
        term = Math.exp(log1term);

        /* in some cases no more computations are needed */
        if (double_equal(term, 0.0)) /* the first term is almost zero */ {
            if ((double) k > (double) n * p) /* at begin or end of the tail? */
                return -log1term / M_LN10 - logNT; /*
                 * end: use just the first
                 * term
                 */
            else
                return -logNT;
        }

        /* compute more terms if needed */
        bin_tail = term;
        for (i = k + 1; i <= n; i++) {
            bin_term = (double) (n - i + 1)
                    * (i < TABSIZE ? (inv[i] != 0.0 ? inv[i]
                    : (inv[i] = 1.0 / (double) i)) : 1.0 / (double) i);

            mult_term = bin_term * p_term;
            term *= mult_term;
            bin_tail += term;
            if (bin_term < 1.0) {
                err = term
                        * ((1.0 - Math.pow(mult_term, (double) (n - i + 1)))
                        / (1.0 - mult_term) - 1.0);

                if (err < tolerance * Math.abs(-Math.log10(bin_tail) - logNT)
                        * bin_tail)
                    break;
            }
        }
        return -Math.log10(bin_tail) - logNT;
    }

    private static double log_gamma_lanczos(double x) {
        double[] q = {75122.6331530, 80916.6278952, 36308.2951477,
                8687.24529705, 1168.92649479, 83.8676043424, 2.50662827511};
        double a = (x + 0.5) * Math.log(x + 5.5) - (x + 5.5);
        double b = 0.0;
        int n;

        for (n = 0; n < 7; n++) {
            a -= Math.log(x + (double) n);
            b += q[n] * Math.pow(x, (double) n);
        }
        return a + Math.log(b);
    }


    static double log_gamma_windschitl(double x) {
        return 0.918938533204673
                + (x - 0.5)
                * Math.log(x)
                - x
                + 0.5
                * x
                * Math.log(x * Math.sinh(1 / x) + 1
                / (810.0 * Math.pow(x, 6.0)));
    }


    static double log_gamma(double x) {
        return ((x) > 15.0 ? log_gamma_windschitl(x) : log_gamma_lanczos(x));
    }

    /**
     * Compare doubles by relative error.
     *
     * The resulting rounding error after floating point computations depend on
     * the specific operations done. The same number computed by different
     * algorithms could present different rounding errors. For a useful
     * comparison, an estimation of the relative rounding error should be
     * considered and compared to a factor times EPS. The factor should be
     * related to the cumulated rounding error in the chain of computation.
     * Here, as a simplification, a fixed factor is used.
     */
    /**
     * Doubles relative error factor
     */
    static double RELATIVE_ERROR_FACTOR = 100.0;

    static boolean double_equal(double a, double b) {
        double abs_diff, aa, bb, abs_max;

        /* trivial case */
        if (a == b)
            return true;

        abs_diff = Math.abs(a - b);
        aa = Math.abs(a);
        bb = Math.abs(b);
        abs_max = aa > bb ? aa : bb;

        /*
         * DBL_MIN is the smallest normalized number, thus, the smallest number
         * whose relative error is bounded by DBL_EPSILON. For smaller numbers,
         * the same quantization steps as for DBL_MIN are used. Then, for
         * smaller numbers, a meaningful "relative" error should be computed by
         * dividing the difference by DBL_MIN.
         */
        if (abs_max < -Double.MAX_VALUE)
            abs_max = -Double.MAX_VALUE;

        /* equal if relative error <= factor x eps */
        return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * calculateMachineEpsilonDouble());
    }

    public static double calculateMachineEpsilonDouble() {
        float machEps = 1.0f;

        do
            machEps /= 2.0f;
        while ((double) (1.0 + (machEps / 2.0)) != 1.0);

        return machEps;
    }
    /*----------------------------------------------------------------------------*/

    /**
     * Rectangle points iterator.
     * <p>
     * The integer coordinates of pixels inside a rectangle are
     * iteratively explored. This structure keep track of the fprocess and
     * functions ri_ini(), ri_inc(), ri_end(), and ri_del() are used in
     * the process. An example of how to use the iterator is as follows:
     * \code
     * <p>
     * struct rect * rec = XXX; // some rectangle
     * rect_iter * i;
     * for( i=ri_ini(rec); !ri_end(i); ri_inc(i) )
     * {
     * // your code, using 'i->x' and 'i->y' as coordinates
     * }
     * ri_del(i); // delete iterator
     * <p>
     * \endcode
     * The pixels are explored 'column' by 'column', where we call
     * 'column' a set of pixels with the same x value that are inside the
     * rectangle. The following is an schematic representation of a
     * rectangle, the 'column' being explored is marked by colons, and
     * the current pixel being explored is 'x,y'.
     * \verbatim
     * <p>
     * vx[1],vy[1]
     * *
     * *
     * *
     * ye
     * :  *
     * vx[0],vy[0]           :     *
     * :        *
     * x,y          *
     * :              *
     * :            vx[2],vy[2]
     * :                *
     * y                     ys              *
     * ^                        *           *
     * |                           *       *
     * |                              *   *
     * +---> x                      vx[3],vy[3]
     * <p>
     * \endverbatim
     * The first 'column' to be explored is the one with the smaller x
     * value. Each 'column' is explored starting from the pixel of the
     * 'column' (inside the rectangle) with the smallest y value.
     * <p>
     * The four corners of the rectangle are stored in order that rotates
     * around the corners at the arrays 'vx[]' and 'vy[]'. The first
     * point is always the one with smaller x value.
     * <p>
     * 'x' and 'y' are the coordinates of the pixel being explored. 'ys'
     * and 'ye' are the start and end values of the current column being
     * explored. So, 'ys' < 'ye'.
     */
    static class RectIterator {
        double[] vx; /* rectangle's corner X coordinates in circular order */
        double[] vy; /* rectangle's corner Y coordinates in circular order */
        double ys, ye; /* start and end Y values of current 'column' */
        int x, y; /* coordinates of currently explored pixel */

        RectIterator() {
            vx = new double[4];
            vy = new double[4];
        }
    }

    static class ImageDouble {
        float[] data;
        int width, height;

        ImageDouble(int width, int height) {
            this.width = width;
            this.height = height;

            data = new float[width * height];
        }

        ImageDouble(int width, int height, float[] data) {
            this.width = width;
            this.height = height;

            this.data = data;
        }

    }

    /*----------------------------------------------------------------------------*/

    /**
     * Chained list of coordinates.
     */
    static class CoorList {
        int x, y;
        CoorList next;
    }

    ;

    static class NTupleList {
        int size;
        int max_size;
        int dim;
        double[] values;


        NTupleList(int dim) {

            /* check parameters */
            if (dim == 0)
                error("new_NTupleList: 'dim' must be positive.");

            /* initialize list */
            size = 0;
            max_size = 1;
            this.dim = dim;

            values = new double[dim * max_size];

        }

    }


    static void enlarge_NTupleList(NTupleList n_tuple) {

        n_tuple.max_size *= 2;

        /* realloc memory */

        //System.out.println("THIS IS ACTUALLY WRONG!!!!!!!!!!");
        int oldlen = n_tuple.values.length;


        double[] arr = new double[n_tuple.dim * n_tuple.max_size];


        for (int i = 0; i < oldlen; i++) {
            arr[i] = n_tuple.values[i];
        }

        n_tuple.values = arr;

    }


    static void add_7tuple(NTupleList out, double v1, double v2, double v3,
                           double v4, double v5, double v6, double v7) {
        /* check parameters */
        if (out == null)
            error("add_7tuple: invalid n-tuple input.");
        if (out.dim != 7)
            error("add_7tuple: the n-tuple must be a 7-tuple.");

        /* if needed, alloc more tuples to 'out' */
        if (out.size == out.max_size)
            enlarge_NTupleList(out);
        if (out.values == null)
            error("add_7tuple: invalid n-tuple input.");

        /* add new 7-tuple */
        out.values[out.size * out.dim + 0] = v1;
        out.values[out.size * out.dim + 1] = v2;
        out.values[out.size * out.dim + 2] = v3;
        out.values[out.size * out.dim + 3] = v4;
        out.values[out.size * out.dim + 4] = v5;
        out.values[out.size * out.dim + 5] = v6;
        out.values[out.size * out.dim + 6] = v7;

        /* update number of tuples counter */
        out.size++;
    }

    static class Rect {
        double x1, y1, x2, y2; /* first and second point of the line segment */
        double width; /* rectangle width */
        double x, y; /* center of the rectangle */
        double theta; /* angle */
        double dx, dy; /* (dx,dy) is vector oriented as the line segment */
        double prec; /* tolerance angle */
        double p; /* probability of a point with angle within 'prec' */
    }

    static float[] gaussian_sampler(float[] in, int width, int height, double scale,
                                  double sigma_scale) {
        ImageDouble aux, out;
        NTupleList kernel;
        int N, M, h, n, x, y, i;
        int xc, yc, j, double_x_size, double_y_size;
        double sigma, xx, yy, sum, prec;

        /* check parameters */
        if (in == null || in == null || width == 0 || height == 0)
            error("gaussian_sampler: invalid image.");
        if (scale <= 0.0)
            error("gaussian_sampler: 'scale' must be positive.");
        if (sigma_scale <= 0.0)
            error("gaussian_sampler: 'sigma_scale' must be positive.");

        /* compute new image size and get memory for images */
        if (width * scale > (double) Integer.MAX_VALUE
                || height * scale > (double) Integer.MAX_VALUE)
            error("gaussian_sampler: the output image size exceeds the handled size.");
        N = (int) Math.ceil(width * scale);
        M = (int) Math.ceil(height * scale);

        aux = new ImageDouble(N, height);
        out = new ImageDouble(N, M);

        /* sigma, kernel size and memory for the kernel */
        sigma = scale < 1.0 ? sigma_scale / scale : sigma_scale;
        /*
         * The size of the kernel is selected to guarantee that the the first
         * discarded term is at least 10^prec times smaller than the central
         * value. For that, h should be larger than x, with e^(-x^2/2sigma^2) =
         * 1/10^prec. Then, x = sigma * sqrt( 2 * prec * ln(10) ).
         */
        prec = 3.0;
        h = (int) Math.ceil(sigma * Math.sqrt(2.0 * prec * Math.log(10.0)));
        n = 1 + 2 * h; /* kernel size */
        kernel = new NTupleList(n);

        /* auxiliary double image size variables */
        double_x_size = (int) (2 * width);
        double_y_size = (int) (2 * height);

        /* First subsampling: x axis */
        for (x = 0; x < aux.width; x++) {
            /*
             * x is the coordinate in the new image. xx is the corresponding
             * x-value in the original size image. xc is the integer value, the
             * pixel coordinate of xx.
             */
            xx = (double) x / scale;
            /*
             * coordinate (0.0,0.0) is in the center of pixel (0,0), so the
             * pixel with xc=0 get the values of xx from -0.5 to 0.5
             */
            xc = (int) Math.floor(xx + 0.5);
            gaussian_kernel(kernel, sigma, (double) h + xx - (double) xc);
            /*
             * the kernel must be computed for each x because the fine offset
             * xx-xc is different in each case
             */

            for (y = 0; y < aux.height; y++) {
                sum = 0.0;
                for (i = 0; i < kernel.dim; i++) {
                    j = xc - h + i;

                    /* symmetry boundary condition */
                    while (j < 0)
                        j += double_x_size;
                    while (j >= double_x_size)
                        j -= double_x_size;
                    if (j >= (int) width)
                        j = double_x_size - 1 - j;

                    sum += in[j + y * width] * kernel.values[i];
                }
                aux.data[x + y * aux.width] = (float)sum;
            }
        }

        /* Second subsampling: y axis */
        for (y = 0; y < out.height; y++) {
            /*
             * y is the coordinate in the new image. yy is the corresponding
             * x-value in the original size image. yc is the integer value, the
             * pixel coordinate of xx.
             */
            yy = (double) y / scale;
            /*
             * coordinate (0.0,0.0) is in the center of pixel (0,0), so the
             * pixel with yc=0 get the values of yy from -0.5 to 0.5
             */
            yc = (int) Math.floor(yy + 0.5);
            gaussian_kernel(kernel, sigma, (double) h + yy - (double) yc);
            /*
             * the kernel must be computed for each y because the fine offset
             * yy-yc is different in each case
             */

            for (x = 0; x < out.width; x++) {
                sum = 0.0;
                for (i = 0; i < kernel.dim; i++) {
                    j = yc - h + i;

                    /* symmetry boundary condition */
                    while (j < 0)
                        j += double_y_size;
                    while (j >= double_y_size)
                        j -= double_y_size;
                    if (j >= (int) height)
                        j = double_y_size - 1 - j;

                    sum += aux.data[x + j * aux.width] * kernel.values[i];
                }
                out.data[x + y * out.width] = (float)sum;
            }
        }

        return out.data;
    }

    static void gaussian_kernel(NTupleList kernel, double sigma, double mean) {
        double sum = 0.0;
        double val;
        int i;

        /* check parameters */
        if (kernel == null || kernel.values == null)
            error("gaussian_kernel: invalid n-tuple 'kernel'.");
        if (sigma <= 0.0)
            error("gaussian_kernel: 'sigma' must be positive.");

        /* compute Gaussian kernel */
        if (kernel.max_size < 1)
            enlarge_NTupleList(kernel);

        kernel.size = 1;
        for (i = 0; i < kernel.dim; i++) {
            val = ((double) i - mean) / sigma;
            kernel.values[i] = Math.exp(-0.5 * val * val);
            sum += kernel.values[i];
        }

        /* normalization */
        if (sum >= 0.0)
            for (i = 0; i < kernel.dim; i++)
                kernel.values[i] /= sum;

    }

}
