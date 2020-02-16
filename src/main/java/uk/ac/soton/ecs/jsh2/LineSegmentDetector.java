//package uk.ac.soton.ecs.jsh2;
//
//import org.openimaj.image.FImage;
//import org.openimaj.math.geometry.line.Line2d;
//
//import java.awt.Point;
//import java.util.ArrayList;
//import java.util.List;
//
//class rect_itr {
//    double[] vx; /* rectangle's corner X coordinates in circular order */
//    double[] vy; /* rectangle's corner Y coordinates in circular order */
//    double ys, ye; /* start and end Y values of current 'column' */
//    int x, y; /* coordinates of currently explored pixel */
//
//    rect_itr() {
//        vx = new double[4];
//        vy = new double[4];
//    }
//}
//
//class image_double {
//    double[] data;
//    int xsize, ysize;
//
//    image_double(int xsize, int ysize) {
//        this.xsize = xsize;
//        this.ysize = ysize;
//
//        data = new double[xsize * ysize];
//    }
//
//    image_double(int xsize, int ysize, double[] data) {
//        this.xsize = xsize;
//        this.ysize = ysize;
//
//        this.data = data;
//    }
//
//}
//
//
//class coorlist {
//    int x, y;
//    coorlist next;
//};
//
//public class LineSegmentDetector {
//
//    double[] inv = new double[TABSIZE];
//
//    double M_LN10 = 2.30258509299404568402;
//
//    double M_PI = 3.14159265358979323846;
//
//    int FALSE = 0;
//
//    int TRUE = 1;
//
//    /** Label for pixels with undefined gradient. */
//    double NOTDEF = -1024.0;
//
//    /** 3/2 pi */
//    double M_3_2_PI = 4.71238898038;
//
//    /** 2 pi */
//    double M_2__PI = 6.28318530718;
//
//    /** Label for pixels not used in yet. */
//    int NOTUSED = 0;
//
//    /** Label for pixels already used in detection. */
//    int USED = 1;
//
//
//    /*----------------------------------------------------------------------------*/;
//
//    void error(String msg) {
//        System.err.println("LSD Error: " + msg);
//        throw new RuntimeException("major error");
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Doubles relative error factor
//     */
//    double RELATIVE_ERROR_FACTOR = 100.0;
//
//    private double calculateMachineEpsilonDouble() {
//        float machEps = 1.0f;
//
//        do
//            machEps /= 2.0f;
//        while ((double) (1.0 + (machEps / 2.0)) != 1.0);
//
//        return machEps;
//    }
//
//    boolean double_equal(double a, double b) {
//        double abs_diff, aa, bb, abs_max;
//
//        /* trivial case */
//        if (a == b)
//            return true;
//
//        abs_diff = Math.abs(a - b);
//        aa = Math.abs(a);
//        bb = Math.abs(b);
//        abs_max = aa > bb ? aa : bb;
//
//        if (abs_max < -Double.MAX_VALUE)
//            abs_max = -Double.MAX_VALUE;
//
//        /* equal if relative error <= factor x eps */
//        return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * calculateMachineEpsilonDouble());
//    }
//
//
//    static double dist(double x1, double y1, double x2, double y2) {
//        return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
//    }
//
//
//    class ntuple_list {
//        int size;
//        int max_size;
//        int dim;
//        double[] values;
//
//        ntuple_list(int dim) {
//
//            /* check parameters */
//            if (dim == 0)
//                error("new_ntuple_list: 'dim' must be positive.");
//
//            /* initialize list */
//            size = 0;
//            max_size = 1;
//            this.dim = dim;
//
//            /* get memory for tuples */
//            // n_tuple.values = new ArrayList(); (double *) malloc( *
//            // sizeof(double)
//            // );
//            // if( n_tuple.values == NULL ) error("not enough memory.");
//
//            values = new double[dim * max_size];
//
//        }
//
//    }
//
//    void enlarge_ntuple_list(ntuple_list n_tuple) {
//        /* check parameters */
//        // if( n_tuple == null || n_tuple.values == null || n_tuple.max_size ==
//        // 0 )
//        // error("enlarge_ntuple_list: invalid n-tuple.");
//
//        /* duplicate number of tuples */
//        n_tuple.max_size *= 2;
//
//        /* realloc memory */
//
//        //System.out.println("THIS IS ACTUALLY WRONG!!!!!!!!!!");
//        int oldlen = n_tuple.values.length;
//
//
//        double [] arr = new double[n_tuple.dim * n_tuple.max_size];
//
//
//        for(int i=0; i < oldlen; i++){
//            arr[i] = n_tuple.values[i];
//        }
//
//        n_tuple.values = arr;
//
//    }
//
//
//    void add_7tuple(ntuple_list out, double v1, double v2, double v3,
//                    double v4, double v5, double v6, double v7) {
//        /* check parameters */
//        if (out == null)
//            error("add_7tuple: invalid n-tuple input.");
//        if (out.dim != 7)
//            error("add_7tuple: the n-tuple must be a 7-tuple.");
//
//        /* if needed, alloc more tuples to 'out' */
//        if (out.size == out.max_size)
//            enlarge_ntuple_list(out);
//        if (out.values == null)
//            error("add_7tuple: invalid n-tuple input.");
//
//        /* add new 7-tuple */
//        out.values[out.size * out.dim + 0] = v1;
//        out.values[out.size * out.dim + 1] = v2;
//        out.values[out.size * out.dim + 2] = v3;
//        out.values[out.size * out.dim + 3] = v4;
//        out.values[out.size * out.dim + 4] = v5;
//        out.values[out.size * out.dim + 5] = v6;
//        out.values[out.size * out.dim + 6] = v7;
//
//        /* update number of tuples counter */
//        out.size++;
//    }
//
//
//    class image_char {
//        int[] data;
//        int xsize, ysize;
//
//        image_char(int xsize, int ysize) {
//            this.xsize = xsize;
//            this.ysize = ysize;
//        }
//
//        image_char(int xsize, int ysize, int fill_value) {
//            data = new int[xsize * ysize]; /* create image */
//            int N = xsize * ysize;
//            int i;
//
//            /* initialize */
//            for (i = 0; i < N; i++)
//                data[i] = fill_value;
//            this.xsize = xsize;
//            this.ysize = ysize;
//        }
//    }
//
//    class image_int {
//        int[] data;
//        int xsize, ysize;
//
//        image_int(int xsize, int ysize) {
//            this.xsize = xsize;
//            this.ysize = ysize;
//        }
//
//        image_int(int xsize, int ysize, int fill_value) {
//            data = new int[xsize * ysize]; /* create image */
//            int N = xsize * ysize;
//            int i;
//
//            /* initialize */
//            for (i = 0; i < N; i++)
//                data[i] = fill_value;
//
//        }
//    }
//
//    class rect {
//        double x1, y1, x2, y2; /* first and second point of the line segment */
//        double width; /* rectangle width */
//        double x, y; /* center of the rectangle */
//        double theta; /* angle */
//        double dx, dy; /* (dx,dy) is vector oriented as the line segment */
//        double prec; /* tolerance angle */
//        double p; /* probability of a point with angle within 'prec' */
//    };
//
//    image_double new_image_double_ptr(int xsize, int ysize, double[] data) {
//
//        image_double image = new image_double(xsize, ysize);
//
//        /* check parameters */
//        if (xsize == 0 || ysize == 0)
//            error("new_image_double_ptr: invalid image size.");
//
//        /* set image */
//        image.data = data;
//
//        return image;
//    }
//
//
//    double log_gamma_lanczos(double x) {
//        double[] q = { 75122.6331530, 80916.6278952, 36308.2951477,
//                8687.24529705, 1168.92649479, 83.8676043424, 2.50662827511 };
//        double a = (x + 0.5) * Math.log(x + 5.5) - (x + 5.5);
//        double b = 0.0;
//        int n;
//
//        for (n = 0; n < 7; n++) {
//            a -= Math.log(x + (double) n);
//            b += q[n] * Math.pow(x, (double) n);
//        }
//        return a + Math.log(b);
//    }
//
//
//    double log_gamma_windschitl(double x) {
//        return 0.918938533204673
//                + (x - 0.5)
//                * Math.log(x)
//                - x
//                + 0.5
//                * x
//                * Math.log(x * Math.sinh(1 / x) + 1
//                / (810.0 * Math.pow(x, 6.0)));
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Computes the natural logarithm of the absolute value of the gamma
//     * function of x. When x>15 use log_gamma_windschitl(), otherwise use
//     * log_gamma_lanczos().
//     */
//    double log_gamma(double x) {
//        return ((x) > 15.0 ? log_gamma_windschitl(x) : log_gamma_lanczos(x));
//    }
//
//    static int TABSIZE = 100000;
//
//
//    double nfa(int n, int k, double p, double logNT) {
//
//        double tolerance = 0.1; /* an error of 10% in the result is accepted */
//        double log1term, term, bin_term, mult_term, bin_tail, err, p_term;
//        int i;
//
//        /* check parameters */
//        if (n < 0 || k < 0 || k > n || p <= 0.0 || p >= 1.0)
//            error("nfa: wrong n, k or p values.");
//
//        /* trivial cases */
//        if (n == 0 || k == 0)
//            return -logNT;
//        if (n == k)
//            return -logNT - (double) n * Math.log10(p);
//
//        /* probability term */
//        p_term = p / (1.0 - p);
//
//        log1term = log_gamma((double) n + 1.0) - log_gamma((double) k + 1.0)
//                - log_gamma((double) (n - k) + 1.0) + (double) k * Math.log(p)
//                + (double) (n - k) * Math.log(1.0 - p);
//        term = Math.exp(log1term);
//
//        /* in some cases no more computations are needed */
//        if (double_equal(term, 0.0)) /* the first term is almost zero */
//        {
//            if ((double) k > (double) n * p) /* at begin or end of the tail? */
//                return -log1term / M_LN10 - logNT; /*
//                 * end: use just the first
//                 * term
//                 */
//            else
//                return -logNT; /* begin: the tail is roughly 1 */
//        }
//
//        /* compute more terms if needed */
//        bin_tail = term;
//        for (i = k + 1; i <= n; i++) {
//            bin_term = (double) (n - i + 1)
//                    * (i < TABSIZE ? (inv[i] != 0.0 ? inv[i]
//                    : (inv[i] = 1.0 / (double) i)) : 1.0 / (double) i);
//
//            mult_term = bin_term * p_term;
//            term *= mult_term;
//            bin_tail += term;
//            if (bin_term < 1.0) {
//
//                err = term
//                        * ((1.0 - Math.pow(mult_term, (double) (n - i + 1)))
//                        / (1.0 - mult_term) - 1.0);
//
//                if (err < tolerance * Math.abs(-Math.log10(bin_tail) - logNT)
//                        * bin_tail)
//                    break;
//            }
//        }
//        return -Math.log10(bin_tail) - logNT;
//    }
//
//    void gaussian_kernel(ntuple_list kernel, double sigma, double mean) {
//        double sum = 0.0;
//        double val;
//        int i;
//
//        /* check parameters */
//        if (kernel == null || kernel.values == null)
//            error("gaussian_kernel: invalid n-tuple 'kernel'.");
//        if (sigma <= 0.0)
//            error("gaussian_kernel: 'sigma' must be positive.");
//
//        /* compute Gaussian kernel */
//        if (kernel.max_size < 1)
//            enlarge_ntuple_list(kernel);
//
//        kernel.size = 1;
//        for (i = 0; i < kernel.dim; i++) {
//            val = ((double) i - mean) / sigma;
//            kernel.values[i] = Math.exp(-0.5 * val * val);
//            sum += kernel.values[i];
//        }
//
//        /* normalization */
//        if (sum >= 0.0)
//            for (i = 0; i < kernel.dim; i++)
//                kernel.values[i] /= sum;
//    }
//
//    image_double gaussian_sampler(image_double in, double scale,
//                                  double sigma_scale) {
//        image_double aux, out;
//        ntuple_list kernel;
//        int N, M, h, n, x, y, i;
//        int xc, yc, j, double_x_size, double_y_size;
//        double sigma, xx, yy, sum, prec;
//
//        /* check parameters */
//        if (in == null || in.data == null || in.xsize == 0 || in.ysize == 0)
//            error("gaussian_sampler: invalid image.");
//        if (scale <= 0.0)
//            error("gaussian_sampler: 'scale' must be positive.");
//        if (sigma_scale <= 0.0)
//            error("gaussian_sampler: 'sigma_scale' must be positive.");
//
//        /* compute new image size and get memory for images */
//        if (in.xsize * scale > (double) Integer.MAX_VALUE
//                || in.ysize * scale > (double) Integer.MAX_VALUE)
//            error("gaussian_sampler: the output image size exceeds the handled size.");
//        N = (int) Math.ceil(in.xsize * scale);
//        M = (int) Math.ceil(in.ysize * scale);
//        aux = new image_double(N, in.ysize);
//        out = new image_double(N, M);
//
//        /* sigma, kernel size and memory for the kernel */
//        sigma = scale < 1.0 ? sigma_scale / scale : sigma_scale;
//
//        prec = 3.0;
//        h = (int) Math.ceil(sigma * Math.sqrt(2.0 * prec * Math.log(10.0)));
//        n = 1 + 2 * h; /* kernel size */
//        kernel = new ntuple_list(n);
//
//        /* auxiliary double image size variables */
//        double_x_size = (int) (2 * in.xsize);
//        double_y_size = (int) (2 * in.ysize);
//
//        /* First subsampling: x axis */
//        for (x = 0; x < aux.xsize; x++) {
//            /*
//             * x is the coordinate in the new image. xx is the corresponding
//             * x-value in the original size image. xc is the integer value, the
//             * pixel coordinate of xx.
//             */
//            xx = (double) x / scale;
//            /*
//             * coordinate (0.0,0.0) is in the center of pixel (0,0), so the
//             * pixel with xc=0 get the values of xx from -0.5 to 0.5
//             */
//            xc = (int) Math.floor(xx + 0.5);
//            gaussian_kernel(kernel, sigma, (double) h + xx - (double) xc);
//            /*
//             * the kernel must be computed for each x because the fine offset
//             * xx-xc is different in each case
//             */
//
//            for (y = 0; y < aux.ysize; y++) {
//                sum = 0.0;
//                for (i = 0; i < kernel.dim; i++) {
//                    j = xc - h + i;
//
//                    /* symmetry boundary condition */
//                    while (j < 0)
//                        j += double_x_size;
//                    while (j >= double_x_size)
//                        j -= double_x_size;
//                    if (j >= (int) in.xsize)
//                        j = double_x_size - 1 - j;
//
//                    sum += in.data[j + y * in.xsize] * kernel.values[i];
//                }
//                aux.data[x + y * aux.xsize] = sum;
//            }
//        }
//
//        /* Second subsampling: y axis */
//        for (y = 0; y < out.ysize; y++) {
//            /*
//             * y is the coordinate in the new image. yy is the corresponding
//             * x-value in the original size image. yc is the integer value, the
//             * pixel coordinate of xx.
//             */
//            yy = (double) y / scale;
//            /*
//             * coordinate (0.0,0.0) is in the center of pixel (0,0), so the
//             * pixel with yc=0 get the values of yy from -0.5 to 0.5
//             */
//            yc = (int) Math.floor(yy + 0.5);
//            gaussian_kernel(kernel, sigma, (double) h + yy - (double) yc);
//            /*
//             * the kernel must be computed for each y because the fine offset
//             * yy-yc is different in each case
//             */
//
//            for (x = 0; x < out.xsize; x++) {
//                sum = 0.0;
//                for (i = 0; i < kernel.dim; i++) {
//                    j = yc - h + i;
//
//                    /* symmetry boundary condition */
//                    while (j < 0)
//                        j += double_y_size;
//                    while (j >= double_y_size)
//                        j -= double_y_size;
//                    if (j >= (int) in.ysize)
//                        j = double_y_size - 1 - j;
//
//                    sum += aux.data[x + j * aux.xsize] * kernel.values[i];
//                }
//                out.data[x + y * out.xsize] = sum;
//            }
//        }
//
//        return out;
//    }
//
//    image_double modgrad;
//    coorlist[] mem_p;
//    coorlist list_p;
//
//    image_double ll_angle(image_double in, double threshold,
//            /* coorlist [] list_p, *//* void ** mem_p, */
//            /* image_double modgrad, */int n_bins) {
//        image_double g;
//        int n, p, x, y, adr, i;
//        double com1, com2, gx, gy, norm, norm2;
//        /*
//         * the rest of the variables are used for pseudo-ordering the gradient
//         * magnitude values
//         */
//        int list_count = 0;
//        coorlist[] list;
//
//        coorlist[] range_l_s; /* array of pointers to start of bin list */
//        coorlist[] range_l_e; /* array of pointers to end of bin list */
//        coorlist start;
//        coorlist end;
//        double max_grad = 0.0;
//
//        /* check parameters */
//        if (in == null || in.data == null || in.xsize == 0 || in.ysize == 0)
//            error("ll_angle: invalid image.");
//        if (threshold < 0.0)
//            error("ll_angle: 'threshold' must be positive.");
//        if (list_p == null) {
//            error("ll_angle: null pointer 'list_p'.");
//            // list_p = new coorlist();
//        }
//        // if (mem_p == null)
//        // error("ll_angle: null pointer 'mem_p'.");
//        // if (modgrad == null)
//        // error("ll_angle: null pointer 'modgrad'.");
//        if (n_bins == 0)
//            error("ll_angle: 'n_bins' must be positive.");
//
//        /* image size shortcuts */
//        n = in.ysize;
//        p = in.xsize;
//
//        list = new coorlist[n * p];
//        for (int z = 0; z < n * p; z++) {
//            list[z] = new coorlist();
//        }
//
//        mem_p = list;
//
//        /* allocate output image */
//        g = new image_double(in.xsize, in.ysize);
//
//        /* get memory for the image of gradient modulus */
//        modgrad = new image_double(in.xsize, in.ysize);
//
//        /* get memory for "ordered" list of pixels */
//        // list = new coorlist[n * p];// (struct coorlist *) calloc( (size_t)
//        // (n*p), sizeof(struct coorlist) );
//
//        range_l_s = new coorlist[n_bins];
//        range_l_e = new coorlist[n_bins];
//
//        if (list == null || range_l_s == null || range_l_e == null)
//            error("not enough memory.");
//
//        for (i = 0; i < n_bins; i++) {
//            range_l_s[i] = range_l_e[i] = null;
//        }
//
//        /* 'undefined' on the down and right boundaries */
//        for (x = 0; x < p; x++)
//            g.data[(n - 1) * p + x] = NOTDEF;
//
//        for (y = 0; y < n; y++)
//            g.data[p * y + p - 1] = NOTDEF;
//
//        /* compute gradient on the remaining pixels */
//        for (x = 0; x < p - 1; x++)
//            for (y = 0; y < n - 1; y++) {
//                adr = y * p + x;
//
//                /*
//                 * Norm 2 computation using 2x2 pixel window: A B C D and com1 =
//                 * D-A, com2 = B-C. Then gx = B+D - (A+C) horizontal difference
//                 * gy = C+D - (A+B) vertical difference com1 and com2 are just
//                 * to avoid 2 additions.
//                 */
//                com1 = in.data[adr + p + 1] - in.data[adr];
//                com2 = in.data[adr + 1] - in.data[adr + p];
//
//                gx = com1 + com2; /* gradient x component */
//                gy = com1 - com2; /* gradient y component */
//                norm2 = gx * gx + gy * gy;
//                norm = Math.sqrt(norm2 / 4.0); /* gradient norm */
//
//                modgrad.data[adr] = norm; /* store gradient norm */
//
//                // System.out.println("norm " + norm + " threshold " +
//                // threshold);
//
//                if (norm <= threshold) /* norm too small, gradient no defined */
//                    g.data[adr] = NOTDEF; /* gradient angle not defined */
//                else {
//                    // System.out.println("gradient angle --------------");
//                    /* gradient angle computation */
//                    g.data[adr] = Math.atan2(gx, -gy);
//
//                    /* look for the maximum of the gradient */
//                    if (norm > max_grad)
//                        max_grad = norm;
//                }
//            }
//
//        /* compute histogram of gradient values */
//        for (x = 0; x < p - 1; x++)
//            for (y = 0; y < n - 1; y++) {
//                norm = modgrad.data[y * p + x];
//
//                /* store the point in the right bin according to its norm */
//                i = (int) (norm * (double) n_bins / max_grad);
//                if (i >= n_bins)
//                    i = n_bins - 1;
//                if (range_l_e[i] == null) {
//                    // System.out.println("here1 "+list[list_count]);
//
//                    range_l_s[i] = range_l_e[i] = list[list_count++];
//                } else {
//                    // System.out.println("here2");
//                    range_l_e[i].next = list[list_count];
//                    range_l_e[i] = list[list_count++];
//                }
//                // System.out.println(i + " "+range_l_e[i]);
//                range_l_e[i].x = (int) x;
//                range_l_e[i].y = (int) y;
//                range_l_e[i].next = null;
//            }
//
//        /*
//         * Make the list of pixels (almost) ordered by norm value. It starts by
//         * the larger bin, so the list starts by the pixels with the highest
//         * gradient value. Pixels would be ordered by norm value, up to a
//         * precision given by max_grad/n_bins.
//         */
//        for (i = n_bins - 1; i > 0 && range_l_s[i] == null; i--) {
//
//        }
//
//        //	System.out.println("i val " + i);
//
//        start = range_l_s[i];
//        end = range_l_e[i];
//        if (start != null) {
//            //System.out.println("start not null");
//            while (i > 0) {
//                --i;
//                if (range_l_s[i] != null) {
//                    // System.out.println("range not null");
//
//                    end.next = range_l_s[i];
//                    end = range_l_e[i];
//                }
//            }
//        }
//        list_p = start;
//
//        return g;
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Absolute value angle difference.
//     */
//    double angle_diff(double a, double b) {
//        a -= b;
//        while (a <= -M_PI)
//            a += M_2__PI;
//        while (a > M_PI)
//            a -= M_2__PI;
//        if (a < 0.0)
//            a = -a;
//        return a;
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Signed angle difference.
//     */
//    double angle_diff_signed(double a, double b) {
//        a -= b;
//        while (a <= -M_PI)
//            a += M_2__PI;
//        while (a > M_PI)
//            a -= M_2__PI;
//        return a;
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Is point (x,y) aligned to angle theta, up to precision 'prec'?
//     */
//    boolean isaligned(int x, int y, image_double angles, double theta,
//                      double prec) {
//        double a;
//
//        /* check parameters */
//        if (angles == null || angles.data == null)
//            error("isaligned: invalid image 'angles'.");
//        if (x < 0 || y < 0 || x >= (int) angles.xsize
//                || y >= (int) angles.ysize)
//            error("isaligned: (x,y) out of the image.");
//        if (prec < 0.0)
//            error("isaligned: 'prec' must be positive.");
//
//        /* angle at pixel (x,y) */
//        a = angles.data[x + y * angles.xsize];
//
//        /*
//         * pixels whose level-line angle is not defined are considered as
//         * NON-aligned
//         */
//        if (a == NOTDEF)
//            return false; /*
//             * there is no need to call the function 'double_equal'
//             * here because there is no risk of problems related to
//             * the comparison doubles, we are only interested in the
//             * exact NOTDEF value
//             */
//
//        /* it is assumed that 'theta' and 'a' are in the range [-pi,pi] */
//        theta -= a;
//        if (theta < 0.0)
//            theta = -theta;
//        if (theta > M_3_2_PI) {
//            theta -= M_2__PI;
//            if (theta < 0.0)
//                theta = -theta;
//        }
//
//        return theta <= prec;
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Build a region of pixels that share the same angle, up to a tolerance
//     * 'prec', starting at point (x,y).
//     */
//
//    public int reg_size;
//    public double reg_angle;
//
//    void region_grow(int x, int y, image_double angles, Point[] reg,
//            /* int * reg_size, double * reg_angle, */image_char used, double prec) {
//        double sumdx, sumdy;
//        int xx, yy, i;
//
//        /* check parameters */
//        if (x < 0 || y < 0 || x >= (int) angles.xsize
//                || y >= (int) angles.ysize)
//            error("region_grow: (x,y) out of the image.");
//        if (angles == null || angles.data == null)
//            error("region_grow: invalid image 'angles'.");
//        if (reg == null)
//            error("region_grow: invalid 'reg'.");
//        // if( reg_size == null )
//        // error("region_grow: invalid pointer 'reg_size'.");
//        // if( reg_angle == NULL )
//        // error("region_grow: invalid pointer 'reg_angle'.");
//        // if( used == NULL || used->data == NULL )
//        // error("region_grow: invalid image 'used'.");
//
//        /* first point of the region */
//        reg_size = 1;
//        reg[0].x = x;
//        reg[0].y = y;
//        reg_angle = angles.data[x + y * angles.xsize]; /* region's angle */
//        sumdx = Math.cos(reg_angle);
//        sumdy = Math.sin(reg_angle);
//        used.data[x + y * used.xsize] = USED;
//
//        /* try neighbors as new region points */
//        for (i = 0; i < reg_size; i++)
//            for (xx = reg[i].x - 1; xx <= reg[i].x + 1; xx++)
//                for (yy = reg[i].y - 1; yy <= reg[i].y + 1; yy++) {
//
//                    //System.out.println("here2____"+(xx >= 0 && yy >= 0 && xx < (int) used.xsize
//                    //		&& yy < (int) used.ysize)/*
//                    //		&& used.data[xx + yy * used.xsize] != USED*/+"__________");
//
//                    //System.out.println("::::::::: "+xx+","+yy+"  "+used.xsize+","+used.ysize);
//
//                    if (xx >= 0 && yy >= 0 && xx < (int) used.xsize
//                            && yy < (int) used.ysize
//                            && used.data[xx + yy * used.xsize] != USED
//                            && isaligned(xx, yy, angles, reg_angle, prec)) {
//                        /* add point */
//                        used.data[xx + yy * used.xsize] = USED;
//                        reg[reg_size].x = xx;
//                        reg[reg_size].y = yy;
//                        ++(reg_size);
//
//                        /* update region's angle */
//                        sumdx += Math.cos(angles.data[xx + yy * angles.xsize]);
//                        sumdy += Math.sin(angles.data[xx + yy * angles.xsize]);
//                        reg_angle = Math.atan2(sumdy, sumdx);
//                    }
//                }
//        //System.out.println(">>>regsize " + reg_size);
//
//    }
//
//
//    double get_theta(Point[] reg, int reg_size, double x, double y,
//                     image_double modgrad, double reg_angle, double prec) {
//        double lambda, theta, weight;
//        double Ixx = 0.0;
//        double Iyy = 0.0;
//        double Ixy = 0.0;
//        int i;
//
//        /* check parameters */
//        if (reg == null)
//            error("get_theta: invalid region.");
//        if (reg_size <= 1)
//            error("get_theta: region size <= 1.");
//        if (modgrad == null || modgrad.data == null)
//            error("get_theta: invalid 'modgrad'.");
//        if (prec < 0.0)
//            error("get_theta: 'prec' must be positive.");
//
//        /* compute inertia matrix */
//        for (i = 0; i < reg_size; i++) {
//            weight = modgrad.data[reg[i].x + reg[i].y * modgrad.xsize];
//            Ixx += ((double) reg[i].y - y) * ((double) reg[i].y - y) * weight;
//            Iyy += ((double) reg[i].x - x) * ((double) reg[i].x - x) * weight;
//            Ixy -= ((double) reg[i].x - x) * ((double) reg[i].y - y) * weight;
//        }
//        if (double_equal(Ixx, 0.0) && double_equal(Iyy, 0.0)
//                && double_equal(Ixy, 0.0))
//            error("get_theta: null inertia matrix.");
//
//        /* compute smallest eigenvalue */
//        lambda = 0.5 * (Ixx + Iyy - Math.sqrt((Ixx - Iyy) * (Ixx - Iyy) + 4.0
//                * Ixy * Ixy));
//
//        /* compute angle */
//        theta = Math.abs(Ixx) > Math.abs(Iyy) ? Math.atan2(lambda - Ixx, Ixy)
//                : Math.atan2(Ixy, lambda - Iyy);
//
//        /*
//         * The previous procedure doesn't cares about orientation, so it could
//         * be wrong by 180 degrees. Here is corrected if necessary.
//         */
//        if (angle_diff(theta, reg_angle) > prec)
//            theta += M_PI;
//
//        return theta;
//    }
//
//    boolean reduce_region_radius(Point[] reg, /* int * reg_size, */
//                                 image_double modgrad, /*double reg_angle,*/ double prec, double p,
//                                 rect rec, image_char used, image_double angles, double density_th) {
//        double density, rad1, rad2, rad, xc, yc;
//        int i;
//
//        /* check parameters */
//        if (reg == null)
//            error("reduce_region_radius: invalid pointer 'reg'.");
//        // if( reg_size == null )
//        // error("reduce_region_radius: invalid pointer 'reg_size'.");
//        if (prec < 0.0)
//            error("reduce_region_radius: 'prec' must be positive.");
//        if (rec == null)
//            error("reduce_region_radius: invalid pointer 'rec'.");
//        if (used == null || used.data == null)
//            error("reduce_region_radius: invalid image 'used'.");
//        if (angles == null || angles.data == null)
//            error("reduce_region_radius: invalid image 'angles'.");
//
//        /* compute region points density */
//        density = (double) reg_size
//                / (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);
//
//        /* if the density criterion is satisfied there is nothing to do */
//        if (density >= density_th)
//            return true;
//
//        /* compute region's radius */
//        xc = (double) reg[0].x;
//        yc = (double) reg[0].y;
//        rad1 = dist(xc, yc, rec.x1, rec.y1);
//        rad2 = dist(xc, yc, rec.x2, rec.y2);
//        rad = rad1 > rad2 ? rad1 : rad2;
//
//        /* while the density criterion is not satisfied, remove farther pixels */
//        while (density < density_th) {
//            rad *= 0.75; /* reduce region's radius to 75% of its value */
//
//            /* remove points from the region and update 'used' map */
//            for (i = 0; i < reg_size; i++)
//                if (dist(xc, yc, (double) reg[i].x, (double) reg[i].y) > rad) {
//                    /* point not kept, mark it as NOTUSED */
//                    used.data[reg[i].x + reg[i].y * used.xsize] = NOTUSED;
//                    /* remove point from the region */
//                    reg[i].x = reg[reg_size - 1].x; /*
//                     * if i==*reg_size-1 copy
//                     * itself
//                     */
//                    reg[i].y = reg[reg_size - 1].y;
//                    --(reg_size);
//                    --i; /* to avoid skipping one point */
//                }
//
//            /*
//             * reject if the region is too small. 2 is the minimal region size
//             * for 'region2rect' to work.
//             */
//            if (reg_size < 2)
//                return false;
//
//            /* re-compute rectangle */
//            region2rect(reg, reg_size, modgrad, /*reg_angle,*/ prec, p, rec);
//
//            /* re-compute region points density */
//            density = (double) reg_size
//                    / (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);
//        }
//
//        /* if this point is reached, the density criterion is satisfied */
//        return true;
//    }
//
//
//    void region2rect(Point[] reg, int reg_size, image_double modgrad,
//            /*double reg_angle,*/ double prec, double p, rect rec) {
//        double x, y, dx, dy, l, w, theta, weight, sum, l_min, l_max, w_min, w_max;
//        int i;
//
//        /* check parameters */
//        if (reg == null)
//            error("region2rect: invalid region.");
//        if (reg_size <= 1)
//            error("region2rect: region size <= 1.");
//        if (modgrad == null || modgrad.data == null)
//            error("region2rect: invalid image 'modgrad'.");
//        if (rec == null)
//            error("region2rect: invalid 'rec'.");
//
//
//        x = y = sum = 0.0;
//        for (i = 0; i < reg_size; i++) {
//            weight = modgrad.data[reg[i].x + reg[i].y * modgrad.xsize];
//            x += (double) reg[i].x * weight;
//            y += (double) reg[i].y * weight;
//            sum += weight;
//        }
//        if (sum <= 0.0)
//            error("region2rect: weights sum equal to zero.");
//        x /= sum;
//        y /= sum;
//
//        /* theta */
//        theta = get_theta(reg, reg_size, x, y, modgrad, reg_angle, prec);
//
//
//        dx = Math.cos(theta);
//        dy = Math.sin(theta);
//        l_min = l_max = w_min = w_max = 0.0;
//        for (i = 0; i < reg_size; i++) {
//            l = ((double) reg[i].x - x) * dx + ((double) reg[i].y - y) * dy;
//            w = -((double) reg[i].x - x) * dy + ((double) reg[i].y - y) * dx;
//
//            if (l > l_max)
//                l_max = l;
//            if (l < l_min)
//                l_min = l;
//            if (w > w_max)
//                w_max = w;
//            if (w < w_min)
//                w_min = w;
//        }
//
//        /* store values */
//        rec.x1 = x + l_min * dx;
//        rec.y1 = y + l_min * dy;
//        rec.x2 = x + l_max * dx;
//        rec.y2 = y + l_max * dy;
//        rec.width = w_max - w_min;
//        rec.x = x;
//        rec.y = y;
//        rec.theta = theta;
//        rec.dx = dx;
//        rec.dy = dy;
//        rec.prec = prec;
//        rec.p = p;
//
//        /*
//         * we impose a minimal width of one pixel
//         *
//         * A sharp horizontal or vertical step would produce a perfectly
//         * horizontal or vertical region. The width computed would be zero. But
//         * that corresponds to a one pixels width transition in the image.
//         */
//        if (rec.width < 1.0)
//            rec.width = 1.0;
//    }
//
//
//    boolean refine(Point[] reg, /*int reg_size,*/ image_double modgrad,
//            /*double reg_angle,*/ double prec, double p, rect rec, image_char used,
//                   image_double angles, double density_th) {
//        double angle, ang_d, mean_angle, tau, density, xc, yc, ang_c, sum, s_sum;
//        int i, n;
//
//        //this.reg_size = reg_size;
//        //this.reg_angle = reg_angle;
//
//        /* check parameters */
//        if (reg == null)
//            error("refine: invalid pointer 'reg'.");
//        // if( reg_size == null ) error("refine: invalid pointer 'reg_size'.");
//        if (prec < 0.0)
//            error("refine: 'prec' must be positive.");
//        if (rec == null)
//            error("refine: invalid pointer 'rec'.");
//        if (used == null || used.data == null)
//            error("refine: invalid image 'used'.");
//        if (angles == null || angles.data == null)
//            error("refine: invalid image 'angles'.");
//
//        /* compute region points density */
//        density = (double) reg_size
//                / (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);
//
//        /* if the density criterion is satisfied there is nothing to do */
//        if (density >= density_th)
//            return true;
//
//        /*------ First try: reduce angle tolerance ------*/
//
//        /* compute the new mean angle and tolerance */
//        xc = (double) reg[0].x;
//        yc = (double) reg[0].y;
//        ang_c = angles.data[reg[0].x + reg[0].y * angles.xsize];
//        sum = s_sum = 0.0;
//        n = 0;
//        for (i = 0; i < this.reg_size; i++) {
//            used.data[reg[i].x + reg[i].y * used.xsize] = NOTUSED;
//            if (dist(xc, yc, (double) reg[i].x, (double) reg[i].y) < rec.width) {
//                angle = angles.data[reg[i].x + reg[i].y * angles.xsize];
//                ang_d = angle_diff_signed(angle, ang_c);
//                sum += ang_d;
//                s_sum += ang_d * ang_d;
//                ++n;
//            }
//        }
//        mean_angle = sum / (double) n;
//        tau = 2.0 * Math.sqrt((s_sum - 2.0 * mean_angle * sum) / (double) n
//                + mean_angle * mean_angle); /* 2 * standard deviation */
//
//        /*
//         * find a new region from the same starting point and new angle
//         * tolerance
//         */
//        region_grow(reg[0].x, reg[0].y, angles, reg /* ,reg_size,reg_angle, */,
//                used, tau);
//        /* if the region is too small, reject */
//        if (reg_size < 2)
//            return false;
//
//        /* re-compute rectangle */
//        region2rect(reg, reg_size, modgrad, /*reg_angle,*/ prec, p, rec);
//
//        /* re-compute region points density */
//        density = (double) reg_size
//                / (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);
//
//        /*------ Second try: reduce region radius ------*/
//        if (density < density_th)
//            return reduce_region_radius(reg, /* reg_size, */modgrad, /*reg_angle,*/
//                    prec, p, rec, used, angles, density_th);
//
//        /* if this point is reached, the density criterion is satisfied */
//        return true;
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Interpolate y value corresponding to 'x' value given, in the line 'x1,y1'
//     * to 'x2,y2'; if 'x1=x2' return the smaller of 'y1' and 'y2'.
//     *
//     * The following restrictions are required: - x1 <= x2 - x1 <= x - x <= x2
//     */
//    double inter_low(double x, double x1, double y1, double x2, double y2) {
//        /* check parameters */
//        if (x1 > x2 || x < x1 || x > x2)
//            error("inter_low: unsuitable input, 'x1>x2' or 'x<x1' or 'x>x2'.");
//
//        /* interpolation */
//        if (double_equal(x1, x2) && y1 < y2)
//            return y1;
//        if (double_equal(x1, x2) && y1 > y2)
//            return y2;
//        return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
//    }
//
//
//    double inter_hi(double x, double x1, double y1, double x2, double y2) {
//        /* check parameters */
//        if (x1 > x2 || x < x1 || x > x2)
//            error("inter_hi: unsuitable input, 'x1>x2' or 'x<x1' or 'x>x2'.");
//
//        /* interpolation */
//        if (double_equal(x1, x2) && y1 < y2)
//            return y2;
//        if (double_equal(x1, x2) && y1 > y2)
//            return y1;
//        return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
//    }
//
//
//    void ri_del(rect_itr iter) {
//        if (iter == null)
//            error("ri_del: NULL iterator.");
//        // free( (void *) iter );
//    }
//
//
//    boolean ri_end(rect_itr i) {
//        /* check input */
//        if (i == null)
//            error("ri_end: NULL iterator.");
//
//        /*
//         * if the current x value is larger than the largest x value in the
//         * rectangle (vx[2]), we know the full exploration of the rectangle is
//         * finished.
//         */
//        return (double) (i.x) > i.vx[2];
//    }
//
//
//    void ri_inc(rect_itr i) {
//        /* check input */
//        if (i == null)
//            error("ri_inc: NULL iterator.");
//
//        /*
//         * if not at end of exploration, increase y value for next pixel in the
//         * 'column'
//         */
//        if (!ri_end(i))
//            i.y++;
//
//        /*
//         * if the end of the current 'column' is reached, and it is not the end
//         * of exploration, advance to the next 'column'
//         */
//        while ((double) (i.y) > i.ye && !ri_end(i)) {
//            /* increase x, next 'column' */
//            i.x++;
//
//            /* if end of exploration, return */
//            if (ri_end(i))
//                return;
//
//            if ((double) i.x < i.vx[3])
//                i.ys = inter_low((double) i.x, i.vx[0], i.vy[0], i.vx[3],
//                        i.vy[3]);
//            else
//                i.ys = inter_low((double) i.x, i.vx[3], i.vy[3], i.vx[2],
//                        i.vy[2]);
//
//            if ((double) i.x < i.vx[1])
//                i.ye = inter_hi((double) i.x, i.vx[0], i.vy[0], i.vx[1],
//                        i.vy[1]);
//            else
//                i.ye = inter_hi((double) i.x, i.vx[1], i.vy[1], i.vx[2],
//                        i.vy[2]);
//
//            /* new y */
//            i.y = (int) Math.ceil(i.ys);
//        }
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Create and initialize a rectangle iterator.
//     *
//     * See details in \ref rect_iter
//     */
//    rect_itr ri_ini(rect r) {
//        double[] vx = new double[4];
//        double[] vy = new double[4];
//        int n, offset;
//        rect_itr i;
//
//        /* check parameters */
//        if (r == null)
//            error("ri_ini: invalid rectangle.");
//
//        /* get memory */
//        i = new rect_itr();
//        if (i == null)
//            error("ri_ini: Not enough memory.");
//
//        /*
//         * build list of rectangle corners ordered in a circular way around the
//         * rectangle
//         */
//        vx[0] = r.x1 - r.dy * r.width / 2.0;
//        vy[0] = r.y1 + r.dx * r.width / 2.0;
//        vx[1] = r.x2 - r.dy * r.width / 2.0;
//        vy[1] = r.y2 + r.dx * r.width / 2.0;
//        vx[2] = r.x2 + r.dy * r.width / 2.0;
//        vy[2] = r.y2 - r.dx * r.width / 2.0;
//        vx[3] = r.x1 + r.dy * r.width / 2.0;
//        vy[3] = r.y1 - r.dx * r.width / 2.0;
//
//        if (r.x1 < r.x2 && r.y1 <= r.y2)
//            offset = 0;
//        else if (r.x1 >= r.x2 && r.y1 < r.y2)
//            offset = 1;
//        else if (r.x1 > r.x2 && r.y1 >= r.y2)
//            offset = 2;
//        else
//            offset = 3;
//
//        /* apply rotation of index. */
//        for (n = 0; n < 4; n++) {
//            i.vx[n] = vx[(offset + n) % 4];
//            i.vy[n] = vy[(offset + n) % 4];
//        }
//
//        i.x = (int) Math.ceil(i.vx[0]) - 1;
//        i.y = (int) Math.ceil(i.vy[0]);
//        i.ys = i.ye = -Double.MAX_VALUE;
//
//        /* advance to the first pixel */
//        ri_inc(i);
//
//        return i;
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Compute a rectangle's NFA value.
//     */
//    double rect_nfa(rect rec, image_double angles, double logNT) {
//        rect_itr i;
//        int pts = 0;
//        int alg = 0;
//
//        /* check parameters */
//        if (rec == null)
//            error("rect_nfa: invalid rectangle.");
//        if (angles == null)
//            error("rect_nfa: invalid 'angles'.");
//
//        /* compute the total number of pixels and of aligned points in 'rec' */
//        for (i = ri_ini(rec); !ri_end(i); ri_inc(i))
//            /* rectangle iterator */
//            if (i.x >= 0 && i.y >= 0 && i.x < (int) angles.xsize
//                    && i.y < (int) angles.ysize) {
//                ++pts; /* total number of pixels counter */
//                if (isaligned(i.x, i.y, angles, rec.theta, rec.prec))
//                    ++alg; /* aligned points counter */
//            }
//        // ri_del(i); /* delete iterator */
//
//        return nfa(pts, alg, rec.p, logNT); /* compute NFA value */
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Copy one rectangle structure to another.
//     */
//    void rect_copy(rect in, rect out) {
//        /* check parameters */
//        if (in == null || out == null)
//            error("rect_copy: invalid 'in' or 'out'.");
//
//        /* copy values */
//        out.x1 = in.x1;
//        out.y1 = in.y1;
//        out.x2 = in.x2;
//        out.y2 = in.y2;
//        out.width = in.width;
//        out.x = in.x;
//        out.y = in.y;
//        out.theta = in.theta;
//        out.dx = in.dx;
//        out.dy = in.dy;
//        out.prec = in.prec;
//        out.p = in.p;
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * Try some rectangles variations to improve NFA value. Only if the
//     * rectangle is not meaningful (i.e., log_nfa <= log_eps).
//     */
//    double rect_improve(rect rec, image_double angles, double logNT,
//                        double log_eps) {
//        rect r = new rect();
//        double log_nfa, log_nfa_new;
//        double delta = 0.5;
//        double delta_2 = delta / 2.0;
//        int n;
//
//        log_nfa = rect_nfa(rec, angles, logNT);
//
//        if (log_nfa > log_eps)
//            return log_nfa;
//
//        /* try finer precisions */
//        rect_copy(rec, r);
//        for (n = 0; n < 5; n++) {
//            r.p /= 2.0;
//            r.prec = r.p * M_PI;
//            log_nfa_new = rect_nfa(r, angles, logNT);
//            if (log_nfa_new > log_nfa) {
//                log_nfa = log_nfa_new;
//                rect_copy(r, rec);
//            }
//        }
//
//        if (log_nfa > log_eps)
//            return log_nfa;
//
//        /* try to reduce width */
//        rect_copy(rec, r);
//        for (n = 0; n < 5; n++) {
//            if ((r.width - delta) >= 0.5) {
//                r.width -= delta;
//                log_nfa_new = rect_nfa(r, angles, logNT);
//                if (log_nfa_new > log_nfa) {
//                    rect_copy(r, rec);
//                    log_nfa = log_nfa_new;
//                }
//            }
//        }
//
//        if (log_nfa > log_eps)
//            return log_nfa;
//
//        /* try to reduce one side of the rectangle */
//        rect_copy(rec, r);
//        for (n = 0; n < 5; n++) {
//            if ((r.width - delta) >= 0.5) {
//                r.x1 += -r.dy * delta_2;
//                r.y1 += r.dx * delta_2;
//                r.x2 += -r.dy * delta_2;
//                r.y2 += r.dx * delta_2;
//                r.width -= delta;
//                log_nfa_new = rect_nfa(r, angles, logNT);
//                if (log_nfa_new > log_nfa) {
//                    rect_copy(r, rec);
//                    log_nfa = log_nfa_new;
//                }
//            }
//        }
//
//        if (log_nfa > log_eps)
//            return log_nfa;
//
//        /* try to reduce the other side of the rectangle */
//        rect_copy(rec, r);
//        for (n = 0; n < 5; n++) {
//            if ((r.width - delta) >= 0.5) {
//                r.x1 -= -r.dy * delta_2;
//                r.y1 -= r.dx * delta_2;
//                r.x2 -= -r.dy * delta_2;
//                r.y2 -= r.dx * delta_2;
//                r.width -= delta;
//                log_nfa_new = rect_nfa(r, angles, logNT);
//                if (log_nfa_new > log_nfa) {
//                    rect_copy(r, rec);
//                    log_nfa = log_nfa_new;
//                }
//            }
//        }
//
//        if (log_nfa > log_eps)
//            return log_nfa;
//
//        /* try even finer precisions */
//        rect_copy(rec, r);
//        for (n = 0; n < 5; n++) {
//            r.p /= 2.0;
//            r.prec = r.p * M_PI;
//            log_nfa_new = rect_nfa(r, angles, logNT);
//            if (log_nfa_new > log_nfa) {
//                log_nfa = log_nfa_new;
//                rect_copy(r, rec);
//            }
//        }
//
//        return log_nfa;
//    }
//
//    /**
//     * LSD full interface.
//     */
//    double[] runDetector(double[] img, int X, int Y,
//                                           double scale, double sigma_scale, double quant, double ang_th,
//                                           double log_eps, double density_th, int n_bins) {
//        image_double image;
//        ntuple_list out = new ntuple_list(7);
//
//        list_p = new coorlist();
//
//        double[] return_value;
//        image_double scaled_image, angles;
//        image_char used;
//        image_int region = null;
//
//        rect rec = new rect();
//
//        Point[] reg;
//        // int reg_size = 0,
//        int min_reg_size, i;
//        int xsize, ysize;
//        double rho, reg_angle = 0, prec, p, log_nfa, logNT;
//        int ls_count = 0; /* line segments are numbered 1,2,3,... */
//
//        /* check parameters */
//        if (img == null || X <= 0 || Y <= 0)
//            error("invalid image input.");
//        if (scale <= 0.0)
//            error("'scale' value must be positive.");
//        if (sigma_scale <= 0.0)
//            error("'sigma_scale' value must be positive.");
//        if (quant < 0.0)
//            error("'quant' value must be positive.");
//        if (ang_th <= 0.0 || ang_th >= 180.0)
//            error("'ang_th' value must be in the range (0,180).");
//        if (density_th < 0.0 || density_th > 1.0)
//            error("'density_th' value must be in the range [0,1].");
//        if (n_bins <= 0)
//            error("'n_bins' value must be positive.");
//
//        /* angle tolerance */
//        prec = M_PI * ang_th / 180.0;
//        p = ang_th / 180.0;
//        rho = quant / Math.sin(prec); /* gradient magnitude threshold */
//
//        modgrad = new image_double(X, Y);
//
//        /* load and scale image (if necessary) and compute angle at each pixel */
//        image = new image_double(X, Y, img);
//        if (scale != 1.0) {
//            scaled_image = gaussian_sampler(image, scale, sigma_scale);
//            angles = ll_angle(scaled_image,/* list_p, */rho, /*
//             * &list_p, &mem_p,*
//             * &modgrad,
//             */(int) n_bins);
//            // free_image_double(scaled_image);
//        } else
//            angles = ll_angle(image, /* list_p, */rho,/* list_p, &mem_p, &modgrad, */
//                    (int) n_bins);
//        xsize = angles.xsize;
//        ysize = angles.ysize;
//
//        /*
//         * Number of Tests - NT
//         *
//         * The theoretical number of tests is Np.(XY)^(5/2) where X and Y are
//         * number of columns and rows of the image. Np corresponds to the number
//         * of angle precisions considered. As the procedure 'rect_improve' tests
//         * 5 times to halve the angle precision, and 5 more times after
//         * improving other factors, 11 different precision values are
//         * potentially tested. Thus, the number of tests is 11 * (X*Y)^(5/2)
//         * whose logarithm value is log10(11) + 5/2 * (log10(X) + log10(Y)).
//         */
//        logNT = 5.0 * (Math.log10((double) xsize) + Math.log10((double) ysize))
//                / 2.0 + Math.log10(11.0);
//        min_reg_size = (int) (-logNT / Math.log10(p)); /*
//         * minimal number of
//         * points in region that
//         * can give a meaningful
//         * event
//         */
//
//        /* initialize some structures */
//        // if( reg_img != NULL && reg_x != NULL && reg_y != NULL ) /* save
//        // region data */
//        region = new image_int(angles.xsize, angles.ysize, 0);
//        used = new image_char(xsize, ysize, NOTUSED);
//        reg = new Point[xsize * ysize];
//
//        for (int z = 0; z < xsize * ysize; z++) {
//            reg[z] = new Point();
//        }
//        // if( reg == NULL ) error("not enough memory!");
//
//        //System.out.println("'ere1" + list_p + " " + list_p.next);
//        /* search for line segments */
//        for (; list_p != null; list_p = list_p.next) {
//            // System.out.println("'ere1.5       "+(used.data[list_p.x +
//            // list_p.y * used.xsize] == NOTUSED)+" "+(angles.data[list_p.x +
//            // list_p.y * angles.xsize] != NOTDEF));
//            if (used.data[list_p.x + list_p.y * used.xsize] == NOTUSED
//                    && angles.data[list_p.x + list_p.y * angles.xsize] != NOTDEF)
//                /*
//                 * there is no risk of double comparison problems here because we
//                 * are only interested in the exact NOTDEF value
//                 */
//            {
//                // System.out.println("'ere2");
//
//                /* find the region of connected point and ~equal angle */
//                //System.out.println("attempting to grow " + list_p.x + " "
//                //		+ list_p.y);
//                region_grow(list_p.x, list_p.y, angles, reg, /*
//                 * &reg_size,
//                 * &reg_angle,
//                 */used, prec);
//
//                /* reject small regions */
//                if (reg_size < min_reg_size) {
//                    //System.out.println("regsize " + reg_size + "   min: "
//                    //		+ min_reg_size);
//                    continue;
//                }
//
//                //System.out.println("LINE FOUND HERE");
//
//                /* construct rectangular approximation for the region */
//                region2rect(reg, reg_size, modgrad, /*reg_angle,*/ prec, p, rec);
//
//                /*
//                 * Check if the rectangle exceeds the minimal density of region
//                 * points. If not, try to improve the region. The rectangle will
//                 * be rejected if the final one does not fulfill the minimal
//                 * density condition. This is an addition to the original LSD
//                 * algorithm published in
//                 * "LSD: A Fast Line Segment Detector with a False Detection Control"
//                 * by R. Grompone von Gioi, J. Jakubowicz, J.M. Morel, and G.
//                 * Randall. The original algorithm is obtained with density_th =
//                 * 0.0.
//                 */
//                if (!refine(reg, /*reg_size,*/ modgrad, /*reg_angle,*/ prec, p, rec,
//                        used, angles, density_th))
//                    continue;
//
//                /* compute NFA value */
//                log_nfa = rect_improve(rec, angles, logNT, log_eps);
//                if (log_nfa <= log_eps)
//                    continue;
//
//                /* A New Line Segment was found! */
//                ++ls_count; /* increase line segment counter */
//
//                //System.out.println("LINE FOUND");
//                /*
//                 * The gradient was computed with a 2x2 mask, its value
//                 * corresponds to points with an offset of (0.5,0.5), that
//                 * should be added to output. The coordinates origin is at the
//                 * center of pixel (0,0).
//                 */
//                rec.x1 += 0.5;
//                rec.y1 += 0.5;
//                rec.x2 += 0.5;
//                rec.y2 += 0.5;
//
//                /* scale the result values if a subsampling was performed */
//                if (scale != 1.0) {
//                    rec.x1 /= scale;
//                    rec.y1 /= scale;
//                    rec.x2 /= scale;
//                    rec.y2 /= scale;
//                    rec.width /= scale;
//                }
//
//                /* add line segment found to output */
//
//                //System.out.println(">>>>>>>> "+rec.x1+" "+ rec.y1+" "+ rec.x2+" "+ rec.y2);
//
//                add_7tuple(out, rec.x1, rec.y1, rec.x2, rec.y2, rec.width,
//                        rec.p, log_nfa);
//
//                /* add region number to 'region' image if needed */
//                if (region != null)
//                    for (i = 0; i < reg_size; i++)
//                        region.data[reg[i].x + reg[i].y * region.xsize] = ls_count;
//            }
//        }
//
//        n_out = (int) (out.size);
//
//        return_value = out.values;
//
//        return return_value;
//        // }
//
//    }
//
//    /*----------------------------------------------------------------------------*/
//    /**
//     * LSD Simple Interface with Scale and Region output.
//     */
//    double sigma_scale = 0.6;
//    double quant = 2.0;
//    double ang_th = 22.5;
//    double log_eps = 0.0;
//    double density_th = 0.7;
//    int n_bins = 1024;
//
//    int n_out;
//
//    public List<Line2d> getLines(){
//        List<Line2d> lines = new ArrayList<>();
//
//        double [] out = runDetector(imageData, width, height, 1.0, sigma_scale, quant,
//                ang_th, log_eps, density_th, n_bins);
//
//        for(int i = 0; i < this.n_out; i++) {
//            for (int j = 0; j < 7; j++)
//                lines.add(new Line2d(
//                        (float)out[7 * i],
//                        (float)out[7 * i + 1],
//                        (float)out[7 * i + 2],
//                        (float)out[7 * i + 3]));
//        }
//
//        return lines;
//    }
//
//    int width, height;
//    double [] imageData;
//    /**Constructor */
//
//    public LineSegmentDetector(double[] imageData, int width, int height) {
//        this.width = width;
//        this.height = height;
//        this.imageData = imageData;
//    }
//    public LineSegmentDetector(FImage image){
//        this.width = image.width;
//        this.height = image.height;
//        imageData = image.getDoublePixelVector();
//    }
//
//    /*----------------------------------------------------------------------------*/
//
//}
