package uk.ac.soton.ecs.jsh2;
/*----------------------------------------------------------------------------

  LSD - Line Segment Detector on digital images

  This code is part of the following publication and was subject
  to peer review:

    "LSD: a Line Segment Detector" by Rafael Grompone von Gioi,
    Jeremie Jakubowicz, Jean-Michel Morel, and Gregory Randall,
    Image Processing On Line, 2012. DOI:10.5201/ipol.2012.gjmr-lsd
    http://dx.doi.org/10.5201/ipol.2012.gjmr-lsd

  Copyright (c) 2007-2011 rafael grompone von gioi <grompone@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @file lsd.c
 * <p>
 * LSD module code
 * @author rafael grompone von gioi <grompone@gmail.com>
 * <p>
 * Java port
 * @author chris anfractuosity.com
 * @mainpage LSD code documentation
 * <p>
 * This is an implementation of the Line Segment Detector described
 * in the paper:
 * <p>
 * "LSD: A Fast Line Segment Detector with a False Detection Control"
 * by Rafael Grompone von Gioi, Jeremie Jakubowicz, Jean-Michel Morel,
 * and Gregory Randall, IEEE Transactions on Pattern Analysis and
 * Machine Intelligence, vol. 32, no. 4, pp. 722-732, April, 2010.
 * <p>
 * and in more details in the CMLA Technical Report:
 * <p>
 * "LSD: A Line Segment Detector, Technical Report",
 * by Rafael Grompone von Gioi, Jeremie Jakubowicz, Jean-Michel Morel,
 * Gregory Randall, CMLA, ENS Cachan, 2010.
 * <p>
 * The version implemented here includes some further improvements
 * described in the following publication, of which this code is part:
 * <p>
 * "LSD: a Line Segment Detector" by Rafael Grompone von Gioi,
 * Jeremie Jakubowicz, Jean-Michel Morel, and Gregory Randall,
 * Image Processing On Line, 2012. DOI:10.5201/ipol.2012.gjmr-lsd
 * http://dx.doi.org/10.5201/ipol.2012.gjmr-lsd
 * <p>
 * The module's main function is lsd().
 * <p>
 * The source code is contained in two files: lsd.h and lsd.c.
 * <p>
 * HISTORY:
 * - version 1.6 - nov 2011:
 * - changes in the interface,
 * - max_grad parameter removed,
 * - the factor 11 was added to the number of test
 * to consider the different precision values
 * tested,
 * - a minor bug corrected in the gradient sorting
 * code,
 * - the algorithm now also returns p and log_nfa
 * for each detection,
 * - a minor bug was corrected in the image scaling,
 * - the angle comparison in "isaligned" changed
 * from < to <=,
 * - "eps" variable renamed "log_eps",
 * - "lsd_scale_region" interface was added,
 * - minor changes to comments.
 * - version 1.5 - dec 2010: Changes in 'refine', -W option added,
 * and more comments added.
 * - version 1.4 - jul 2010: lsd_scale interface added and doxygen doc.
 * - version 1.3 - feb 2010: Multiple bug correction and improved code.
 * - version 1.2 - dec 2009: First full Ansi C Language version.
 * - version 1.1 - sep 2009: Systematic subsampling to scale 0.8 and
 * correction to partially handle "angle problem".
 * - version 1.0 - jan 2009: First complete Megawave2 and Ansi C Language
 * version.
 * @author rafael grompone von gioi <grompone@gmail.com>
 */
/*----------------------------------------------------------------------------*/

/** @mainpage LSD code documentation

This is an implementation of the Line Segment Detector described
in the paper:

"LSD: A Fast Line Segment Detector with a False Detection Control"
by Rafael Grompone von Gioi, Jeremie Jakubowicz, Jean-Michel Morel,
and Gregory Randall, IEEE Transactions on Pattern Analysis and
Machine Intelligence, vol. 32, no. 4, pp. 722-732, April, 2010.

and in more details in the CMLA Technical Report:

"LSD: A Line Segment Detector, Technical Report",
by Rafael Grompone von Gioi, Jeremie Jakubowicz, Jean-Michel Morel,
Gregory Randall, CMLA, ENS Cachan, 2010.

The version implemented here includes some further improvements
described in the following publication, of which this code is part:

"LSD: a Line Segment Detector" by Rafael Grompone von Gioi,
Jeremie Jakubowicz, Jean-Michel Morel, and Gregory Randall,
Image Processing On Line, 2012. DOI:10.5201/ipol.2012.gjmr-lsd
http://dx.doi.org/10.5201/ipol.2012.gjmr-lsd

The module's main function is lsd().

The source code is contained in two files: lsd.h and lsd.c.

HISTORY:
- version 1.6 - nov 2011:
- changes in the interface,
- max_grad parameter removed,
- the factor 11 was added to the number of test
to consider the different precision values
tested,
- a minor bug corrected in the gradient sorting
code,
- the algorithm now also returns p and log_nfa
for each detection,
- a minor bug was corrected in the image scaling,
- the angle comparison in "isaligned" changed
from < to <=,
- "eps" variable renamed "log_eps",
- "lsd_scale_region" interface was added,
- minor changes to comments.
- version 1.5 - dec 2010: Changes in 'refine', -W option added,
and more comments added.
- version 1.4 - jul 2010: lsd_scale interface added and doxygen doc.
- version 1.3 - feb 2010: Multiple bug correction and improved code.
- version 1.2 - dec 2009: First full Ansi C Language version.
- version 1.1 - sep 2009: Systematic subsampling to scale 0.8 and
correction to partially handle "angle problem".
- version 1.0 - jan 2009: First complete Megawave2 and Ansi C Language
version.

 @author rafael grompone von gioi <grompone@gmail.com>
 */
/*----------------------------------------------------------------------------*/

import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.math.geometry.line.Line2d;

import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.awt.Point;
import java.util.Arrays;
import java.util.List;

/*----------------------------------------------------------------------------*/

/** Rectangle points iterator.

 The integer coordinates of pixels inside a rectangle are
 iteratively explored. This structure keep track of the fprocess and
 functions ri_ini(), ri_inc(), ri_end(), and ri_del() are used in
 the process. An example of how to use the iterator is as follows:
 \code

 struct rect * rec = XXX; // some rectangle
 rect_iter * i;
 for( i=ri_ini(rec); !ri_end(i); ri_inc(i) )
 {
 // your code, using 'i->x' and 'i->y' as coordinates
 }
 ri_del(i); // delete iterator

 \endcode
 The pixels are explored 'column' by 'column', where we call
 'column' a set of pixels with the same x value that are inside the
 rectangle. The following is an schematic representation of a
 rectangle, the 'column' being explored is marked by colons, and
 the current pixel being explored is 'x,y'.
 \verbatim

 vx[1],vy[1]
 *   *
 *       *
 *           *
 *               ye
 *                :  *
 vx[0],vy[0]           :     *
 *              :        *
 *          x,y          *
 *        :              *
 *     :            vx[2],vy[2]
 *  :                *
 y                     ys              *
 ^                        *           *
 |                           *       *
 |                              *   *
 +---> x                      vx[3],vy[3]

 \endverbatim
 The first 'column' to be explored is the one with the smaller x
 value. Each 'column' is explored starting from the pixel of the
 'column' (inside the rectangle) with the smallest y value.

 The four corners of the rectangle are stored in order that rotates
 around the corners at the arrays 'vx[]' and 'vy[]'. The first
 point is always the one with smaller x value.

 'x' and 'y' are the coordinates of the pixel being explored. 'ys'
 and 'ye' are the start and end values of the current column being
 explored. So, 'ys' < 'ye'.
 */
class RectIterator {
    double[] vx; /* rectangle's corner X coordinates in circular order */
    double[] vy; /* rectangle's corner Y coordinates in circular order */
    double ys, ye; /* start and end Y values of current 'column' */
    int x, y; /* coordinates of currently explored pixel */

    RectIterator() {
        vx = new double[4];
        vy = new double[4];
    }
}

class ImageDouble {
    double[] data;
    int width, height;

    ImageDouble(int width, int height) {
        this.width = width;
        this.height = height;

        data = new double[width * height];
    }

    ImageDouble(int width, int height, double[] data) {
        this.width = width;
        this.height = height;

        this.data = data;
    }

}

/*----------------------------------------------------------------------------*/

/**
 * Chained list of coordinates.
 */
class CoorList {
    int x, y;
    CoorList next;
};

public class LSD {

    private double[] inv = new double[TABSIZE]; /*

    /** ln(10) */
    double M_LN10 = 2.30258509299404568402;


    /** PI */
    double M_PI = 3.14159265358979323846;

    /** Label for pixels with undefined gradient. */
    double NOTDEF = -1024.0;

    /** 3/2 pi */
    double M_3_2_PI = 4.71238898038;

    /** 2 pi */
    double M_2__PI = 6.28318530718;

    /** Label for pixels not used in yet. */
    boolean NOTUSED = false;

    /** Label for pixels already used in detection. */
    boolean USED = true;


    void error(String msg) {
        System.err.println("LSD Error: " + msg);
        throw new RuntimeException("major error");
    }

    /*----------------------------------------------------------------------------*/
    /**
     * Doubles relative error factor
     */
    double RELATIVE_ERROR_FACTOR = 100.0;

    private double calculateMachineEpsilonDouble() {
        float machEps = 1.0f;

        do
            machEps /= 2.0f;
        while ((double) (1.0 + (machEps / 2.0)) != 1.0);

        return machEps;
    }

    /*----------------------------------------------------------------------------*/

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
    boolean double_equal(double a, double b) {
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

    /*----------------------------------------------------------------------------*/

    /**
     * Computes Euclidean distance between point (x1,y1) and point (x2,y2).
     */
    static double dist(double x1, double y1, double x2, double y2) {
        return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    }

    /*----------------------------------------------------------------------------*/
    /*----------------------- 'list of n-tuple' data type ------------------------*/
    /*----------------------------------------------------------------------------*/

    /*----------------------------------------------------------------------------*/

    /**
     * 'list of n-tuple' data type
     *
     * The i-th component of the j-th n-tuple of an n-tuple list 'ntl' is
     * accessed with:
     *
     * ntl.values[ i + j * ntl.dim ]
     *
     * The dimension of the n-tuple (n) is:
     *
     * ntl.dim
     *
     * The number of n-tuples in the list is:
     *
     * ntl.size
     *
     * The maximum number of n-tuples that can be stored in the list with the
     * allocated memory at a given time is given by:
     *
     * ntl.max_size
     */
    class NTupleList {
        int size;
        int max_size;
        int dim;
        double[] values;

        /*----------------------------------------------------------------------------*/

        /**
         * Create an n-tuple list and allocate memory for one element.
         *
         * @param dim
         *            the dimension (n) of the n-tuple.
         */
        NTupleList(int dim) {

            /* check parameters */
            if (dim == 0)
                error("new_ntuple_list: 'dim' must be positive.");

            /* initialize list */
            size = 0;
            max_size = 1;
            this.dim = dim;

            /* get memory for tuples */
            // n_tuple.values = new ArrayList(); (double *) malloc( *
            // sizeof(double)
            // );
            // if( n_tuple.values == NULL ) error("not enough memory.");

            values = new double[dim * max_size];

        }

    }

    /*----------------------------------------------------------------------------*/

    /**
     * Enlarge the allocated memory of an n-tuple list.
     */
    void enlarge_ntuple_list(NTupleList n_tuple) {
        /* check parameters */
        // if( n_tuple == null || n_tuple.values == null || n_tuple.max_size ==
        // 0 )
        // error("enlarge_ntuple_list: invalid n-tuple.");

        /* duplicate number of tuples */
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

    /*----------------------------------------------------------------------------*/

    /**
     * Add a 7-tuple to an n-tuple list.
     */
    void add_7tuple(NTupleList out, double v1, double v2, double v3,
					double v4, double v5, double v6, double v7) {
        /* check parameters */
        if (out == null)
            error("add_7tuple: invalid n-tuple input.");
        if (out.dim != 7)
            error("add_7tuple: the n-tuple must be a 7-tuple.");

        /* if needed, alloc more tuples to 'out' */
        if (out.size == out.max_size)
            enlarge_ntuple_list(out);
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

    /*----------------------------------------------------------------------------*/
    /*----------------------------- Image Data Types -----------------------------*/
    /*----------------------------------------------------------------------------*/

    /*----------------------------------------------------------------------------*/

    /**
     * char image data type
     *
     * The pixel value at (x,y) is accessed by:
     *
     * image.data[ x + y * image.xsize ]
     *
     * with x and y integer.
     */
	static class ImageChar {
        int[] data;
        int xsize, ysize;

		ImageChar(int xsize, int ysize, int fill_value) {
            data = new int[xsize * ysize]; /* create image */
            int N = xsize * ysize;
            int i;

            /* initialize */
            for (i = 0; i < N; i++)
                data[i] = fill_value;
            this.xsize = xsize;
            this.ysize = ysize;
        }
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

    ;

    /*----------------------------------------------------------------------------*/
    /*----------------------------- NFA computation ------------------------------*/
    /*----------------------------------------------------------------------------*/

    /*----------------------------------------------------------------------------*/

    /**
     * Computes the natural logarithm of the absolute value of the gamma
     * function of x using the Lanczos approximation. See
     * http://www.rskey.org/gamma.htm
     *
     * The formula used is
     *
     * @f[ \Gamma(x) = \frac{ \sum_{n=0}^{N} q_n x^n }{ \Pi_{n=0}^{N} (x+n) }
     *     (x+5.5)^{x+0.5} e^{-(x+5.5)}
     * @f] so
     * @f[ \log\Gamma(x) = \log\left( \sum_{n=0}^{N} q_n x^n \right) + (x+0.5)
     *     \log(x+5.5) - (x+5.5) - \sum_{n=0}^{N} \log(x+n)
     * @f] and q0 = 75122.6331530, q1 = 80916.6278952, q2 = 36308.2951477, q3 =
     *     8687.24529705, q4 = 1168.92649479, q5 = 83.8676043424, q6 =
     *     2.50662827511.
     */
    private double log_gamma_lanczos(double x) {
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

    /*----------------------------------------------------------------------------*/

    /**
     * Computes the natural logarithm of the absolute value of the gamma
     * function of x using Windschitl method. See http://www.rskey.org/gamma.htm
     *
     * The formula used is
     *
     * @f[ \Gamma(x) = \sqrt{\frac{2\pi}{x}} \left( \frac{x}{e} \sqrt{
     *     x\sinh(1/x) + \frac{1}{810x^6} } \right)^x
     * @f] so
     * @f[ \log\Gamma(x) = 0.5\log(2\pi) + (x-0.5)\log(x) - x + 0.5x\log\left(
     *     x\sinh(1/x) + \frac{1}{810x^6} \right).
     * @f] This formula is a good approximation when x > 15.
     */
    private double log_gamma_windschitl(double x) {
        return 0.918938533204673
                + (x - 0.5)
                * Math.log(x)
                - x
                + 0.5
                * x
                * Math.log(x * Math.sinh(1 / x) + 1
                / (810.0 * Math.pow(x, 6.0)));
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Computes the natural logarithm of the absolute value of the gamma
     * function of x. When x>15 use log_gamma_windschitl(), otherwise use
     * log_gamma_lanczos().
     */
    double log_gamma(double x) {
        return ((x) > 15.0 ? log_gamma_windschitl(x) : log_gamma_lanczos(x));
    }

    /*----------------------------------------------------------------------------*/
    /**
     * Size of the table to store already computed inverse values.
     */
    static int TABSIZE = 100000;

    /*----------------------------------------------------------------------------*/

    /**
     * Computes -log10(NFA).
     *
     * NFA stands for Number of False Alarms:
     *
     * @f[ \mathrm{NFA} = NT \cdot B(n,k,p)
     * @f]
     *
     *     - NT - number of tests - B(n,k,p) - tail of binomial distribution
     *     with parameters n,k and p:
     * @f[ B(n, k, p) = \sum_{j=k}^n \left(\begin{array}{c}n\\j\end{array}\right)
     *     p^{j} (1-p)^{n-j}
     * @f]
     *
     *     The value -log10(NFA) is equivalent but more intuitive than NFA: - -1
     *     corresponds to 10 mean false alarms - 0 corresponds to 1 mean false
     *     alarm - 1 corresponds to 0.1 mean false alarms - 2 corresponds to
     *     0.01 mean false alarms - ...
     *
     *     Used this way, the bigger the value, better the detection, and a
     *     logarithmic scale is used.
     * @param n
     *            ,k,p binomial parameters.
     * @param logNT
     *            logarithm of Number of Tests
     *
     *            The computation is based in the gamma function by the
     *            following relation:
     * @f[ \left(\begin{array}{c}n\\k\end{array}\right) = \frac{ \Gamma(n+1) }{
     *     \Gamma(k+1) \cdot \Gamma(n-k+1) }.
     * @f] We use efficient algorithms to compute the logarithm of the gamma
     *     function.
     *
     *     To make the computation faster, not all the sum is computed, part of
     *     the terms are neglected based on a bound to the error obtained (an
     *     error of 10% in the result is accepted).
     */
    double nfa(int n, int k, double p, double logNT) {

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

    private ImageDouble modgrad;
	private CoorList list_p;

    private ImageDouble ll_angle(double threshold, int n_bins) {
        ImageDouble g;
        int n, p, x, y, adr, i;
        double com1, com2, gx, gy, norm, norm2;
        /*
         * the rest of the variables are used for pseudo-ordering the gradient
         * magnitude values
         */
        int list_count = 0;
        CoorList[] list;

        CoorList[] range_l_s; /* array of pointers to start of bin list */
        CoorList[] range_l_e; /* array of pointers to end of bin list */
        CoorList start;
        CoorList end;
        double max_grad = 0.0;

        /* check parameters */
        if (imageData == null || imageData == null )
            error("ll_angle: invalid image.");
        if (threshold < 0.0)
            error("ll_angle: 'threshold' must be positive.");
        if (list_p == null) {
            error("ll_angle: null pointer 'list_p'.");
        }
        if (n_bins == 0)
            error("ll_angle: 'n_bins' must be positive.");

        /* image size shortcuts */
        n = height;
        p = width;

        list = new CoorList[n * p];
        for (int z = 0; z < n * p; z++) {
            list[z] = new CoorList();
        }

		/* allocate output image */
        g = new ImageDouble(width, height);

        /* get memory for the image of gradient modulus */
        modgrad = new ImageDouble(width, height);

        range_l_s = new CoorList[n_bins];
        range_l_e = new CoorList[n_bins];

        if (list == null || range_l_s == null || range_l_e == null)
            error("not enough memory.");

        for (i = 0; i < n_bins; i++) {
            range_l_s[i] = range_l_e[i] = null;
        }

        /* 'undefined' on the down and right boundaries */
        for (x = 0; x < p; x++)
            g.data[(n - 1) * p + x] = NOTDEF;

        for (y = 0; y < n; y++)
            g.data[p * y + p - 1] = NOTDEF;

        /* compute gradient on the remaining pixels */
        for (x = 0; x < p - 1; x++)
            for (y = 0; y < n - 1; y++) {
                adr = y * p + x;

                /*
                 * Norm 2 computation using 2x2 pixel window: A B C D and com1 =
                 * D-A, com2 = B-C. Then gx = B+D - (A+C) horizontal difference
                 * gy = C+D - (A+B) vertical difference com1 and com2 are just
                 * to avoid 2 additions.
                 */
                com1 = imageData[adr + p + 1] - imageData[adr];
                com2 = imageData[adr + 1] - imageData[adr + p];

                gx = com1 + com2; /* gradient x component */
                gy = com1 - com2; /* gradient y component */
                norm2 = gx * gx + gy * gy;
                norm = Math.sqrt(norm2 / 4.0); /* gradient norm */

                modgrad.data[adr] = norm; /* store gradient norm */

                // System.out.println("norm " + norm + " threshold " +
                // threshold);

                if (norm <= threshold) /* norm too small, gradient no defined */
                    g.data[adr] = NOTDEF; /* gradient angle not defined */
                else {
                    // System.out.println("gradient angle --------------");
                    /* gradient angle computation */
                    g.data[adr] = Math.atan2(gx, -gy);

                    /* look for the maximum of the gradient */
                    if (norm > max_grad)
                        max_grad = norm;
                }
            }

        /* compute histogram of gradient values */
        for (x = 0; x < p - 1; x++)
            for (y = 0; y < n - 1; y++) {
                norm = modgrad.data[y * p + x];

                /* store the point imageData the right bin according to its norm */
                i = (int) (norm * (double) n_bins / max_grad);
                if (i >= n_bins)
                    i = n_bins - 1;
                if (range_l_e[i] == null) {
                    // System.out.println("here1 "+list[list_count]);

                    range_l_s[i] = range_l_e[i] = list[list_count++];
                } else {
                    // System.out.println("here2");
                    range_l_e[i].next = list[list_count];
                    range_l_e[i] = list[list_count++];
                }
                // System.out.println(i + " "+range_l_e[i]);
                range_l_e[i].x = (int) x;
                range_l_e[i].y = (int) y;
                range_l_e[i].next = null;
            }

        /*
         * Make the list of pixels (almost) ordered by norm value. It starts by
         * the larger bin, so the list starts by the pixels with the highest
         * gradient value. Pixels would be ordered by norm value, up to a
         * precision given by max_grad/n_bins.
         */
        for (i = n_bins - 1; i > 0 && range_l_s[i] == null; i--) {
        	;
        }

        //	System.out.println("i val " + i);

        start = range_l_s[i];
        end = range_l_e[i];
        if (start != null) {
            //System.out.println("start not null");
            while (i > 0) {
                --i;
                if (range_l_s[i] != null) {
                    // System.out.println("range not null");

                    end.next = range_l_s[i];
                    end = range_l_e[i];
                }
            }
        }
        list_p = start;

        return g;
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Absolute value angle difference.
     */
    double angle_diff(double a, double b) {
        a -= b;
        while (a <= -M_PI)
            a += M_2__PI;
        while (a > M_PI)
            a -= M_2__PI;
        if (a < 0.0)
            a = -a;
        return a;
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Signed angle difference.
     */
    double angle_diff_signed(double a, double b) {
        a -= b;
        while (a <= -M_PI)
            a += M_2__PI;
        while (a > M_PI)
            a -= M_2__PI;
        return a;
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Is point (x,y) aligned to angle theta, up to precision 'prec'?
     */
	private boolean isAligned(int x, int y, ImageDouble angles, double theta,
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

    /*----------------------------------------------------------------------------*/
    /**
     * Build a region of pixels that share the same angle, up to a tolerance
     * 'prec', starting at point (x,y).
     */

    public int reg_size;
    public double reg_angle;

    void region_grow(int x, int y, ImageDouble angles, Point[] reg,
					 boolean[] used, double prec) {
        double sumdx, sumdy;
        int xx, yy, i;

        /* check parameters */
        if (x < 0 || y < 0 || x >= (int) angles.width
                || y >= (int) angles.height)
            error("region_grow: (x,y) out of the image.");
        if (angles == null || angles.data == null)
            error("region_grow: invalid image 'angles'.");
        if (reg == null)
            error("region_grow: invalid 'reg'.");

        /* first point of the region */
        reg_size = 1;
        reg[0].x = x;
        reg[0].y = y;
        reg_angle = angles.data[x + y * angles.width]; /* region's angle */
        sumdx = Math.cos(reg_angle);
        sumdy = Math.sin(reg_angle);
        used[x + y * width] = USED;

        /* try neighbors as new region points */
        for (i = 0; i < reg_size; i++)
            for (xx = reg[i].x - 1; xx <= reg[i].x + 1; xx++)
                for (yy = reg[i].y - 1; yy <= reg[i].y + 1; yy++) {
                    if (xx >= 0 && yy >= 0 && xx < width
                            && yy < height
                            && used[xx + yy * width] != USED
                            && isAligned(xx, yy, angles, reg_angle, prec)) {
                        /* add point */
                        used[xx + yy * width] = USED;
                        reg[reg_size].x = xx;
                        reg[reg_size].y = yy;
                        ++(reg_size);

                        /* update region's angle */
                        sumdx += Math.cos(angles.data[xx + yy * angles.width]);
                        sumdy += Math.sin(angles.data[xx + yy * angles.width]);
                        reg_angle = Math.atan2(sumdy, sumdx);
                    }
                }
        //System.out.println(">>>regsize " + reg_size);

    }

    /**
     * Compute region's angle as the principal inertia axis of the region.
     */
	private double getTheta(Point[] reg, int reg_size, double x, double y,
							ImageDouble modgrad, double reg_angle, double prec) {
        double lambda, theta, weight;
        double Ixx = 0.0;
        double Iyy = 0.0;
        double Ixy = 0.0;
        int i;

        /* check parameters */
        if (reg == null)
            error("get_theta: invalid region.");
        if (reg_size <= 1)
            error("get_theta: region size <= 1.");
        if (modgrad == null || modgrad.data == null)
            error("get_theta: invalid 'modgrad'.");
        if (prec < 0.0)
            error("get_theta: 'prec' must be positive.");

        /* compute inertia matrix */
        for (i = 0; i < reg_size; i++) {
            weight = modgrad.data[reg[i].x + reg[i].y * modgrad.width];
            Ixx += ((double) reg[i].y - y) * ((double) reg[i].y - y) * weight;
            Iyy += ((double) reg[i].x - x) * ((double) reg[i].x - x) * weight;
            Ixy -= ((double) reg[i].x - x) * ((double) reg[i].y - y) * weight;
        }
        if (double_equal(Ixx, 0.0) && double_equal(Iyy, 0.0)
                && double_equal(Ixy, 0.0))
            error("get_theta: null inertia matrix.");

        /* compute smallest eigenvalue */
        lambda = 0.5 * (Ixx + Iyy - Math.sqrt((Ixx - Iyy) * (Ixx - Iyy) + 4.0
                * Ixy * Ixy));

        /* compute angle */
        theta = Math.abs(Ixx) > Math.abs(Iyy) ? Math.atan2(lambda - Ixx, Ixy)
                : Math.atan2(Ixy, lambda - Iyy);

        /*
         * The previous procedure doesn't cares about orientation, so it could
         * be wrong by 180 degrees. Here is corrected if necessary.
         */
        if (angle_diff(theta, reg_angle) > prec)
            theta += M_PI;

        return theta;
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Reduce the region size, by elimination the points far from the starting
     * point, until that leads to rectangle with the right density of region
     * points or to discard the region if too small.
     */
    boolean reduceRegionRadius(Point[] reg, ImageDouble modgrad, double prec, double p,
							   Rect rec, boolean used[], ImageDouble angles, double density_th) {
        double density, rad1, rad2, rad, xc, yc;
        int i;

        /* check parameters */
        if (reg == null)
            error("reduce_region_radius: invalid pointer 'reg'.");
        // if( reg_size == null )
        // error("reduce_region_radius: invalid pointer 'reg_size'.");
        if (prec < 0.0)
            error("reduce_region_radius: 'prec' must be positive.");
        if (rec == null)
            error("reduce_region_radius: invalid pointer 'rec'.");
        if (used == null || used == null)
            error("reduce_region_radius: invalid image 'used'.");
        if (angles == null || angles.data == null)
            error("reduce_region_radius: invalid image 'angles'.");

        /* compute region points density */
        density = (double) reg_size
                / (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);

        /* if the density criterion is satisfied there is nothing to do */
        if (density >= density_th)
            return true;

        /* compute region's radius */
        xc = reg[0].x;
        yc = reg[0].y;
        rad1 = dist(xc, yc, rec.x1, rec.y1);
        rad2 = dist(xc, yc, rec.x2, rec.y2);
        rad = Math.max(rad1, rad2);

        /* while the density criterion is not satisfied, remove farther pixels */
        while (density < density_th) {
            rad *= 0.75; /* reduce region's radius to 75% of its value */

            /* remove points from the region and update 'used' map */
            for (i = 0; i < reg_size; i++)
                if (dist(xc, yc, (double) reg[i].x, (double) reg[i].y) > rad) {
                    /* point not kept, mark it as NOTUSED */
                    used[reg[i].x + reg[i].y * width] = NOTUSED;
                    /* remove point from the region */
                    reg[i].x = reg[reg_size - 1].x; /*
                     * if i==*reg_size-1 copy
                     * itself
                     */
                    reg[i].y = reg[reg_size - 1].y;
                    --(reg_size);
                    --i; /* to avoid skipping one point */
                }

            /*
             * reject if the region is too small. 2 is the minimal region size
             * for 'region2rect' to work.
             */
            if (reg_size < 2)
                return false;

            /* re-compute rectangle */
            regionToRect(reg, reg_size, modgrad, /*reg_angle,*/ prec, p, rec);

            /* re-compute region points density */
            density = (double) reg_size
                    / (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);
        }

        /* if this point is reached, the density criterion is satisfied */
        return true;
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Computes a rectangle that covers a region of points.
     */
	private void regionToRect(Point[] reg, int reg_size, ImageDouble modgrad,
			/*double reg_angle,*/ double prec, double p, Rect rec) {
        double x, y, dx, dy, l, w, theta, weight, sum, l_min, l_max, w_min, w_max;
        int i;

        /* check parameters */
        if (reg == null)
            error("region2rect: invalid region.");
        if (reg_size <= 1)
            error("region2rect: region size <= 1.");
        if (modgrad == null || modgrad.data == null)
            error("region2rect: invalid image 'modgrad'.");
        if (rec == null)
            error("region2rect: invalid 'rec'.");

        x = y = sum = 0.0;
        for (i = 0; i < reg_size; i++) {
            weight = modgrad.data[reg[i].x + reg[i].y * modgrad.width];
            x += (double) reg[i].x * weight;
            y += (double) reg[i].y * weight;
            sum += weight;
        }
        if (sum <= 0.0)
            error("region2rect: weights sum equal to zero.");
        x /= sum;
        y /= sum;

        /* theta */
        theta = getTheta(reg, reg_size, x, y, modgrad, reg_angle, prec);

        dx = Math.cos(theta);
        dy = Math.sin(theta);
        l_min = l_max = w_min = w_max = 0.0;
        for (i = 0; i < reg_size; i++) {
            l = ((double) reg[i].x - x) * dx + ((double) reg[i].y - y) * dy;
            w = -((double) reg[i].x - x) * dy + ((double) reg[i].y - y) * dx;

            if (l > l_max)
                l_max = l;
            if (l < l_min)
                l_min = l;
            if (w > w_max)
                w_max = w;
            if (w < w_min)
                w_min = w;
        }

        /* store values */
        rec.x1 = x + l_min * dx;
        rec.y1 = y + l_min * dy;
        rec.x2 = x + l_max * dx;
        rec.y2 = y + l_max * dy;
        rec.width = w_max - w_min;
        rec.x = x;
        rec.y = y;
        rec.theta = theta;
        rec.dx = dx;
        rec.dy = dy;
        rec.prec = prec;
        rec.p = p;


        if (rec.width < 1.0)
            rec.width = 1.0;
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Refine a rectangle.
     *
     * For that, an estimation of the angle tolerance is performed by the
     * standard deviation of the angle at points near the region's starting
     * point. Then, a new region is grown starting from the same point, but
     * using the estimated angle tolerance. If this fails to produce a rectangle
     * with the right density of region points, 'reduce_region_radius' is called
     * to try to satisfy this condition.
     */
    private boolean refine(Point[] reg, ImageDouble modgrad,
            double prec, double p, Rect rec, boolean used[],
                           ImageDouble angles, double density_th) {
        double angle, ang_d, mean_angle, tau, density, xc, yc, ang_c, sum, s_sum;
        int i, n;

        //this.reg_size = reg_size;
        //this.reg_angle = reg_angle;

        /* check parameters */
        if (reg == null)
            error("refine: invalid pointer 'reg'.");
        // if( reg_size == null ) error("refine: invalid pointer 'reg_size'.");
        if (prec < 0.0)
            error("refine: 'prec' must be positive.");
        if (rec == null)
            error("refine: invalid pointer 'rec'.");
        if (used == null || used == null)
            error("refine: invalid image 'used'.");
        if (angles == null || angles.data == null)
            error("refine: invalid image 'angles'.");

        /* compute region points density */
        density = (double) reg_size
                / (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);

        /* if the density criterion is satisfied there is nothing to do */
        if (density >= density_th)
            return true;

        /*------ First try: reduce angle tolerance ------*/

        /* compute the new mean angle and tolerance */
        xc = reg[0].x;
        yc = reg[0].y;
        ang_c = angles.data[reg[0].x + reg[0].y * angles.width];
        sum = s_sum = 0.0;
        n = 0;
        for (i = 0; i < this.reg_size; i++) {
            used[reg[i].x + reg[i].y * width] = NOTUSED;
            if (dist(xc, yc, (double) reg[i].x, (double) reg[i].y) < rec.width) {
                angle = angles.data[reg[i].x + reg[i].y * angles.width];
                ang_d = angle_diff_signed(angle, ang_c);
                sum += ang_d;
                s_sum += ang_d * ang_d;
                ++n;
            }
        }
        mean_angle = sum / (double) n;
        tau = 2.0 * Math.sqrt((s_sum - 2.0 * mean_angle * sum) / (double) n
                + mean_angle * mean_angle); /* 2 * standard deviation */

        /*
         * find a new region from the same starting point and new angle
         * tolerance
         */
        region_grow(reg[0].x, reg[0].y, angles, reg /* ,reg_size,reg_angle, */,
                used, tau);
        /* if the region is too small, reject */
        if (reg_size < 2)
            return false;

        /* re-compute rectangle */
        regionToRect(reg, reg_size, modgrad, /*reg_angle,*/ prec, p, rec);

        /* re-compute region points density */
        density = (double) reg_size
                / (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);

        /*------ Second try: reduce region radius ------*/
        if (density < density_th)
            return reduceRegionRadius(reg, modgrad,
                    prec, p, rec, used, angles, density_th);

        /* if this point is reached, the density criterion is satisfied */
        return true;
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Interpolate y value corresponding to 'x' value given, in the line 'x1,y1'
     * to 'x2,y2'; if 'x1=x2' return the smaller of 'y1' and 'y2'.
     *
     * The following restrictions are required: - x1 <= x2 - x1 <= x - x <= x2
     */
    private double inter_low(double x, double x1, double y1, double x2, double y2) {
        /* check parameters */
        if (x1 > x2 || x < x1 || x > x2)
            error("inter_low: unsuitable input, 'x1>x2' or 'x<x1' or 'x>x2'.");

        /* interpolation */
        if (double_equal(x1, x2) && y1 < y2)
            return y1;
        if (double_equal(x1, x2) && y1 > y2)
            return y2;
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Interpolate y value corresponding to 'x' value given, in the line 'x1,y1'
     * to 'x2,y2'; if 'x1=x2' return the larger of 'y1' and 'y2'.
     *
     * The following restrictions are required: - x1 <= x2 - x1 <= x - x <= x2
     */
    double inter_hi(double x, double x1, double y1, double x2, double y2) {
        /* check parameters */
        if (x1 > x2 || x < x1 || x > x2)
            error("inter_hi: unsuitable input, 'x1>x2' or 'x<x1' or 'x>x2'.");

        /* interpolation */
        if (double_equal(x1, x2) && y1 < y2)
            return y2;
        if (double_equal(x1, x2) && y1 > y2)
            return y1;
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Free memory used by a rectangle iterator.
     */
    void ri_del(RectIterator iter) {
        if (iter == null)
            error("ri_del: NULL iterator.");
        // free( (void *) iter );
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Check if the iterator finished the full iteration.
     *
     * See details in \ref rect_iter
     */
    boolean ri_end(RectIterator i) {
        /* check input */
        if (i == null)
            error("ri_end: NULL iterator.");

        /*
         * if the current x value is larger than the largest x value in the
         * rectangle (vx[2]), we know the full exploration of the rectangle is
         * finished.
         */
        return (double) (i.x) > i.vx[2];
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Increment a rectangle iterator.
     *
     * See details in \ref rect_iter
     */
    void ri_inc(RectIterator i) {
        /* check input */
        if (i == null)
            error("ri_inc: NULL iterator.");

        /*
         * if not at end of exploration, increase y value for next pixel in the
         * 'column'
         */
        if (!ri_end(i))
            i.y++;

        /*
         * if the end of the current 'column' is reached, and it is not the end
         * of exploration, advance to the next 'column'
         */
        while ((double) (i.y) > i.ye && !ri_end(i)) {
            /* increase x, next 'column' */
            i.x++;

            /* if end of exploration, return */
            if (ri_end(i))
                return;

            if ((double) i.x < i.vx[3])
                i.ys = inter_low((double) i.x, i.vx[0], i.vy[0], i.vx[3],
                        i.vy[3]);
            else
                i.ys = inter_low((double) i.x, i.vx[3], i.vy[3], i.vx[2],
                        i.vy[2]);


            if ((double) i.x < i.vx[1])
                i.ye = inter_hi((double) i.x, i.vx[0], i.vy[0], i.vx[1],
                        i.vy[1]);
            else
                i.ye = inter_hi((double) i.x, i.vx[1], i.vy[1], i.vx[2],
                        i.vy[2]);

            /* new y */
            i.y = (int) Math.ceil(i.ys);
        }
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Create and initialize a rectangle iterator.
     *
     * See details in \ref rect_iter
     */
    private RectIterator ri_init(Rect r) {
        double[] vx = new double[4];
        double[] vy = new double[4];
        int n, offset;
        RectIterator i;

        /* check parameters */
        if (r == null)
            error("ri_ini: invalid rectangle.");

        /* get memory */
        i = new RectIterator();
        if (i == null)
            error("ri_ini: Not enough memory.");


        vx[0] = r.x1 - r.dy * r.width / 2.0;
        vy[0] = r.y1 + r.dx * r.width / 2.0;
        vx[1] = r.x2 - r.dy * r.width / 2.0;
        vy[1] = r.y2 + r.dx * r.width / 2.0;
        vx[2] = r.x2 + r.dy * r.width / 2.0;
        vy[2] = r.y2 - r.dx * r.width / 2.0;
        vx[3] = r.x1 + r.dy * r.width / 2.0;
        vy[3] = r.y1 - r.dx * r.width / 2.0;

        if (r.x1 < r.x2 && r.y1 <= r.y2)
            offset = 0;
        else if (r.x1 >= r.x2 && r.y1 < r.y2)
            offset = 1;
        else if (r.x1 > r.x2 && r.y1 >= r.y2)
            offset = 2;
        else
            offset = 3;

        /* apply rotation of index. */
        for (n = 0; n < 4; n++) {
            i.vx[n] = vx[(offset + n) % 4];
            i.vy[n] = vy[(offset + n) % 4];
        }

        i.x = (int) Math.ceil(i.vx[0]) - 1;
        i.y = (int) Math.ceil(i.vy[0]);
        i.ys = i.ye = -Double.MAX_VALUE;

        /* advance to the first pixel */
        ri_inc(i);

        return i;
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Compute a rectangle's NFA value.
     */
	private double calculateRectNFA(Rect rec, ImageDouble angles, double logNT) {
        RectIterator i;
        int pts = 0;
        int alg = 0;

        /* check parameters */
        if (rec == null)
            error("rect_nfa: invalid rectangle.");
        if (angles == null)
            error("rect_nfa: invalid 'angles'.");

        /* compute the total number of pixels and of aligned points in 'rec' */
        for (i = ri_init(rec); !ri_end(i); ri_inc(i))
            /* rectangle iterator */
            if (i.x >= 0 && i.y >= 0 && i.x < (int) angles.width
                    && i.y < (int) angles.height) {
                ++pts; /* total number of pixels counter */
                if (isAligned(i.x, i.y, angles, rec.theta, rec.prec))
                    ++alg; /* aligned points counter */
            }
        // ri_del(i); /* delete iterator */

        return nfa(pts, alg, rec.p, logNT); /* compute NFA value */
    }

    /*----------------------------------------------------------------------------*/

    /**
     * Copy one rectangle structure to another.
     */
	private void copyRect(Rect in, Rect out) {
        /* check parameters */
        if (in == null || out == null)
            error("rect_copy: invalid 'in' or 'out'.");

        /* copy values */
        out.x1 = in.x1;
        out.y1 = in.y1;
        out.x2 = in.x2;
        out.y2 = in.y2;
        out.width = in.width;
        out.x = in.x;
        out.y = in.y;
        out.theta = in.theta;
        out.dx = in.dx;
        out.dy = in.dy;
        out.prec = in.prec;
        out.p = in.p;
    }


    /**
     * Try some rectangles variations to improve NFA value. Only if the
     * rectangle is not meaningful (i.e., log_nfa <= log_eps).
     */
	private double improveRect(Rect rec, ImageDouble angles, double logNT,
							   double log_eps) {
        Rect r = new Rect();
        double log_nfa, log_nfa_new;
        double delta = 0.5;
        double delta_2 = delta / 2.0;
        int n;

        log_nfa = calculateRectNFA(rec, angles, logNT);

        if (log_nfa > log_eps)
            return log_nfa;

        /* try finer precisions */
        copyRect(rec, r);
        for (n = 0; n < 5; n++) {
            r.p /= 2.0;
            r.prec = r.p * M_PI;
            log_nfa_new = calculateRectNFA(r, angles, logNT);
            if (log_nfa_new > log_nfa) {
                log_nfa = log_nfa_new;
                copyRect(r, rec);
            }
        }

        if (log_nfa > log_eps)
            return log_nfa;

        /* try to reduce width */
        copyRect(rec, r);
        for (n = 0; n < 5; n++) {
            if ((r.width - delta) >= 0.5) {
                r.width -= delta;
                log_nfa_new = calculateRectNFA(r, angles, logNT);
                if (log_nfa_new > log_nfa) {
                    copyRect(r, rec);
                    log_nfa = log_nfa_new;
                }
            }
        }

        if (log_nfa > log_eps)
            return log_nfa;

        /* try to reduce one side of the rectangle */
        copyRect(rec, r);
        for (n = 0; n < 5; n++) {
            if ((r.width - delta) >= 0.5) {
                r.x1 += -r.dy * delta_2;
                r.y1 += r.dx * delta_2;
                r.x2 += -r.dy * delta_2;
                r.y2 += r.dx * delta_2;
                r.width -= delta;
                log_nfa_new = calculateRectNFA(r, angles, logNT);
                if (log_nfa_new > log_nfa) {
                    copyRect(r, rec);
                    log_nfa = log_nfa_new;
                }
            }
        }

        if (log_nfa > log_eps)
            return log_nfa;

        /* try to reduce the other side of the rectangle */
        copyRect(rec, r);
        for (n = 0; n < 5; n++) {
            if ((r.width - delta) >= 0.5) {
                r.x1 -= -r.dy * delta_2;
                r.y1 -= r.dx * delta_2;
                r.x2 -= -r.dy * delta_2;
                r.y2 -= r.dx * delta_2;
                r.width -= delta;
                log_nfa_new = calculateRectNFA(r, angles, logNT);
                if (log_nfa_new > log_nfa) {
                    copyRect(r, rec);
                    log_nfa = log_nfa_new;
                }
            }
        }

        if (log_nfa > log_eps)
            return log_nfa;

        /* try even finer precisions */
        copyRect(rec, r);
        for (n = 0; n < 5; n++) {
            r.p /= 2.0;
            r.prec = r.p * M_PI;
            log_nfa_new = calculateRectNFA(r, angles, logNT);
            if (log_nfa_new > log_nfa) {
                log_nfa = log_nfa_new;
                copyRect(r, rec);
            }
        }

        return log_nfa;
    }


    /**
     * LSD full interface.
     */
    private double[] LineSegmentDetection(double quant, double ang_th,
                                          double log_eps, double density_th, int n_bins) {

        NTupleList out = new NTupleList(7);

        list_p = new CoorList();

        double[] return_value;
        ImageDouble angles;
        boolean used[];

        Rect rec = new Rect();

        Point[] reg;

        int min_reg_size;

        double rho, prec, p, log_nfa, logNT;

        /* check parameters */
        if (imageData == null || width <= 0 || height <= 0)
            error("invalid image input.");
        if (quant < 0.0)
            error("'quant' value must be positive.");
        if (ang_th <= 0.0 || ang_th >= 180.0)
            error("'ang_th' value must be in the range (0,180).");
        if (density_th < 0.0 || density_th > 1.0)
            error("'density_th' value must be in the range [0,1].");
        if (n_bins <= 0)
            error("'n_bins' value must be positive.");

        /* angle tolerance */
        prec = M_PI * ang_th / 180.0;
        p = ang_th / 180.0;
        rho = quant / Math.sin(prec); /* gradient magnitude threshold */

        modgrad = new ImageDouble(width, height);

        angles = ll_angle(rho, n_bins);

        logNT = 5.0 * (Math.log10(width) + Math.log10(height))
                / 2.0 + Math.log10(11.0);
        min_reg_size = (int) (-logNT / Math.log10(p));

//        used = new ImageChar(width, height, NOTUSED);
		used = new boolean[width*height];
		Arrays.fill(used, NOTUSED);

        reg = new Point[width * height];

        for (int z = 0; z < width * height; z++) {
            reg[z] = new Point();
        }


        /* search for line segments */
        for (; list_p != null; list_p = list_p.next) {

            if (used[list_p.x + list_p.y * width] == NOTUSED
                    && angles.data[list_p.x + list_p.y * angles.width] != NOTDEF) {

                region_grow(list_p.x, list_p.y, angles, reg, used, prec);

                /* reject small regions */
                if (reg_size < min_reg_size) {
                    continue;
                }

                //System.out.println("LINE FOUND HERE");

                /* construct rectangular approximation for the region */
                regionToRect(reg, reg_size, modgrad, prec, p, rec);

                if (!refine(reg, modgrad, prec, p, rec, used, angles, density_th))
                    continue;

                /* compute NFA value */
                log_nfa = improveRect(rec, angles, logNT, log_eps);
                if (log_nfa <= log_eps)
                    continue;

                /* A New Line Segment was found! */
//                ++ls_count; /* increase line segment counter */

                rec.x1 += 0.5;
                rec.y1 += 0.5;
                rec.x2 += 0.5;
                rec.y2 += 0.5;

                /* add line segment found to output */

                //System.out.println(">>>>>>>> "+rec.x1+" "+ rec.y1+" "+ rec.x2+" "+ rec.y2);

                add_7tuple(out, rec.x1, rec.y1, rec.x2, rec.y2, rec.width,
                        rec.p, log_nfa);

            }
        }

        n_out = out.size;

        return_value = out.values;

        return return_value;

    }

    /*----------------------------------------------------------------------------*/

    /**
     * LSD Simple Interface with Scale and Region output.
     */
    double[] lsd_scale_region() {
        /* LSD parameters */
		double quant = 2.0; /*
         * Bound to the quantization error on the gradient
         * norm.
         */
        double ang_th = 22.5; /* Gradient angle tolerance in degrees. */
        double log_eps = 0.0; /* Detection threshold: -log10(NFA) > log_eps */
        double density_th = 0.7; /*
         * Minimal density of region points in
         * rectangle.
         */
        int n_bins = 1024; /*
         * Number of bins in pseudo-ordering of gradient
         * modulus.
         */

        return LineSegmentDetection(quant,
                ang_th, log_eps, density_th, n_bins);

    }

    /*----------------------------------------------------------------------------*/

    /**
     * LSD Simple Interface with Scale.
     */
    double[] lsd_scale() {
        return lsd_scale_region();

    }

    int n_out;

    /*----------------------------------------------------------------------------*/

    /**
     * LSD Simple Interface.
     */
    double[] lsd() {
        return lsd_scale();
    }

    /*----------------------------------------------------------------------------*/
    private int width, height;
    private double[] imageData;

    //	private double scale=1.0;
    public List<Line2d> getLines() {
        List<Line2d> lines = new ArrayList<>();

        double[] out = lsd();

        for (int i = 0; i < this.n_out; i++) {
            for (int j = 0; j < 7; j++) {
                Line2d l = new Line2d(
                        (float) out[7 * i],
                        (float) out[7 * i + 1],
                        (float) out[7 * i + 2],
                        (float) out[7 * i + 3]);
//                if(l.calculateLength() > 150)
                lines.add(l);
            }
        }

        return lines;
    }

    /**Constructor */
    public LSD(double[] imageData, int width, int height) {
        this.width = width;
        this.height = height;
        this.imageData = imageData;
    }

    public LSD(FImage image) {
        this.width = image.width;
        this.height = image.height;
        imageData = image.getDoublePixelVector();
    }

    public LSD(MBFImage image) {
        this.width = image.getWidth();
        this.height = image.getHeight();
//		imageData = getDoublePixelVector(image);

        BufferedImage bf = ImageUtilities.createBufferedImageForDisplay(image);
        double[] tmp = bf.getData().getPixels(0, 0, width, height, new double[width * height * 3]);

        double[] arr2 = new double[width * height];

        int c = 0;
        for (int i = 0; i < tmp.length - 3; i += 3) {
            double B = tmp[i];
            double G = tmp[i + 1];
            double R = tmp[i + 2];
            double level = R * 0.2126 + G * 0.7152 + B * 0.0722;
            arr2[c++] = level;
        }
        imageData = arr2;

    }

    public double[] getDoublePixelVector(MBFImage image) {
        int height = image.getHeight();
        int width = image.getWidth();

        double[] arr2 = new double[width * height];

        float[][] pixelsB = image.bands.get(0).pixels;
        float[][] pixelsR = image.bands.get(1).pixels;
        float[][] pixelsG = image.bands.get(2).pixels;

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                double R = pixelsR[y][x];
                double G = pixelsG[y][x];
                double B = pixelsB[y][x];
                arr2[x + y * width] = R * 0.2126 + G * 0.7152 + B * 0.0722;
            }
        }

        return arr2;
    }
}
