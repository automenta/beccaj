package becca.core_mtj;

import com.github.fommil.netlib.BLAS;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Random;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrices;
import org.apache.commons.math3.util.Precision;



public class Util  {

    //OPENBLAS_NUM_THREADS=4
    //-Dcom.github.fommil.netlib.BLAS=com.github.fommil.netlib.NativeSystemBLAS
    static {
        System.setProperty("com.github.fommil.netlib.BLAS", "com.github.fommil.netlib.NativeSystemBLAS");
        System.out.println("BLAS: " + BLAS.getInstance());
    }
        
    
    //https://raw.githubusercontent.com/fommil/matrix-toolkits-java/master/src/main/java/no/uib/cipr/matrix/Matrix.java
    //https://raw.githubusercontent.com/fommil/matrix-toolkits-java/master/src/main/java/no/uib/cipr/matrix/Matrices.java
    //https://github.com/fommil/matrix-toolkits-java/tree/master/src/test/java/no/uib/cipr/matrix

    public final static double EPSILON = 2.0 * Precision.EPSILON;   // // sys.float_info.epsilon = 2.220446049250313e-16
    public final static double BIG = Math.pow(10, 20);
    public final static int MAX_INT16 = 32767;     // np.iinfo(np.int16).max

    static int log(int x, int base) {
        return (int) (Math.log(x) / Math.log(base));
    }

    public static String m(DenseMatrix m) {
        return "[" + m.numRows() + "," + m.numColumns() + "]";
    }

    public static DenseMatrix multMatrixMatrix(DenseMatrix a, DenseMatrix b) {
        assert (a.numRows() == b.numRows());
        assert (a.numColumns() == b.numColumns());

        final DenseMatrix c = new DenseMatrix(a.numRows(), a.numColumns());
        final double[] cd = c.getData();
        final double[] ad = a.getData();
        final double[] bd = b.getData();
        for (int i = 0; i < cd.length; i++) {
            cd[i] = ad[i] * bd[i];
        }

        return c;
    }
    
    public static DenseVector sumCols(DenseMatrix m, DenseVector target) {
        if (target == null)
            target = new DenseVector(m.numColumns());
        
        final double[] d = target.getData();
        
        for (int c = 0; c < m.numColumns(); c++) {
            for (int r = 0; r < m.numColumns(); r++) {
                d[c] += m.get(r, c);                
            }
        }
        return target;
    }
    public static DenseVector sumRows(DenseMatrix m, DenseVector target) {
        if (target == null)
            target = new DenseVector(m.numRows());
        
        final double[] d = target.getData();
        
        for (int r = 0; r < m.numColumns(); r++) {
            for (int c = 0; c < m.numColumns(); c++) {
                d[r] += m.get(r, c);                
            }
        }
        return target;
    }    

    public static DenseVector getWeightedAverage(DenseMatrix values, DenseMatrix weights) {
        //""" Perform a weighted average of values, using weights """

        assert (values.numRows() == weights.numRows());

        if (values.numColumns() < weights.numColumns()) {
            values = broadcastRows(values, weights.numColumns());
        } else if (values.numColumns() > weights.numColumns()) {
            weights = broadcastRows(weights, values.numColumns());
        }

        //weighted_sum_values = np.sum(values * weights, axis=0)                 
        DenseMatrix valueWeightProduct = multMatrixMatrix(values, weights);

        final DenseVector weightedSumValues = sumCols(valueWeightProduct, null);
        
        
        //sum_of_weights = np.sum(weights, axis=0)         
        final DenseVector sumOfWeights = sumCols(weights, null);

        //return (weighted_sum_values / (sum_of_weights + EPSILON))[:,np.newaxis]
        final double[] wsd = weightedSumValues.getData();
        final double[] sowd = sumOfWeights.getData();

        assert (weightedSumValues.size() == sumOfWeights.size());
        for (int i = 0; i < wsd.length; i++) {
            wsd[i] /= (sowd[i] + EPSILON);
        }
        return weightedSumValues;
    }

    public static void printMatrixDimensions(DenseMatrix... m) {
        String s = "";
        for (DenseMatrix r : m) {
            s += m(r) + " ";
        }
        System.out.println("matrixDim: " + s);
    }

    public static DenseMatrix eleMultiply(DenseMatrix a, DenseMatrix b, DenseMatrix c) {
        /*extern void dsbmv_(const char *uplo,
                           const int *n,
                           const int *k,
                           const double *alpha,
                           const double *a,
                           const int *lda,
                           const double *b,
                           const int *incx,
                           const double *beta,
                           double *c,
                           const int *incy);

        static const int k = 0; // Just the diagonal; 0 super-diagonal bands
        static const double alpha = 1.0;
        static const int lda = 1;
        static const int incx = 1;
        static const double beta = 0.0;
        static const int incy = 1;

        dsbmv_("L", &n, &k, &alpha, a, &lda, b, &incx, &beta, c, &incy);*/
        
        //    public abstract void dsbmv(String string, int i, int i1, double d, double[] doubles, int i2, double[] doubles1, int i3, double d1, double[] doubles2, int i4);
        
        assert(a.numRows() == b.numRows());
        assert(a.numColumns() == b.numColumns());
        if (c == null)
            c = new DenseMatrix(a.numRows(), b.numColumns());
        
        int n = a.numColumns() * a.numRows();
        BLAS.getInstance().dsbmv("N", n, 0, 1.0, a.getData(), 1, b.getData(), 1, 0.0, c.getData(), 1);
        return c;
    }

    
    public static DenseMatrix addEquals(DenseMatrix m, double x) {
        final double []d = m.getData();
        for (int i = 0; i < d.length; i++)
            d[i] += x;
        return m;
    }
    public static DenseVector addEquals(DenseVector m, double x) {
        final double []d = m.getData();
        for (int i = 0; i < d.length; i++)
            d[i] += x;
        return m;
    }    
    public static DenseVector powerEquals(DenseVector m, double x) {
        final double []d = m.getData();
        for (int i = 0; i < d.length; i++)
            d[i] = Math.pow(d[i], x);
        return m;
    }    
    
    public void identity(DenseMatrix d)    {
        int N;
        d.zero();
        if (d.numRows() < d.numColumns())
        {
            N = d.numRows();
        }
        else
        {
            N = d.numColumns();
        }

        for (int i = 0; i < N; i++)
        {
            d.set(i, i, 1.0);
        }
    }
  
    public static DenseMatrix matrixVector(final DenseMatrix matrix, final DenseMatrix vector, DenseMatrix result) {
        //ex: (8, 32) * (1, 32) -> (8, 32)

        if (result==null)        
            result = new DenseMatrix(matrix.numRows(), matrix.numColumns());
        matrix.mult(vector, result);
        return result;
    }

    public static double dot(DenseMatrix a, DenseMatrix b) {
        double[] ad = a.getData();
        double[] bd = b.getData();
        assert (ad.length == bd.length);
        assert (a.numRows() == b.numRows());
        assert (a.numColumns() == b.numColumns());
        double sum = 0;
        for (int i = 0; i < ad.length; i++) {
            sum += ad[i] * bd[i];
        }
        return sum;
    }

    /**
     * the result is a 1-D matrix which can be transposed depending on the
     * context
     */
    public static DenseVector getGeneralizedMean(final DenseMatrix values, final DenseMatrix weights, final double exponent) {
        final DenseMatrix shiftedValues = values.copy();        
        addEquals(shiftedValues, 1.0);

        final DenseMatrix valuesToPower = shiftedValues;
        matrixPower(valuesToPower, exponent);

        //mean_values_to_power = weighted_average(values_to_power, weights)
        final DenseVector meanValuesToPower = getWeightedAverage(valuesToPower, weights);

        //shifted_mean = (mean_values_to_power + EPSILON) ** (1./exponent)
        final DenseVector shiftedMean = meanValuesToPower;
        addEquals(shiftedMean, EPSILON);
        powerEquals(shiftedMean, (1.0 / exponent));

        final DenseVector mean = shiftedMean;

        //mean = shifted_mean - 1.
        addEquals(mean, -1);

        //# Find means for which all weights are zero. These are undefined.
        //# Set them equal to zero.
        //sum_weights = np.sum(weights, axis=0)
        final DenseVector sumWeights = sumCols(weights, null);

        //zero_indices = np.where(np.abs(sum_weights) < EPSILON)
        //mean[zero_indices] = 0.        
        final double[] meanD = mean.getData();
        final double[] sumWeightsD = sumWeights.getData();

        assert (sumWeights.size() == mean.size());
        for (int i = 0; i < sumWeights.size(); i++) {
            if (Math.abs(sumWeightsD[i]) < EPSILON) {
                meanD[i] = 0;
            }
        }

        return mean;
    }

    //A->B, B can be null in which case it operates just on A
    public static double[] mapOneToInf(final double[] a, double[] b) {
        //""" Map values from [0, 1] onto [0, inf) and map values from [-1, 0] onto (-inf, 0]
        final double eps = 2.2204460492503131e-16; //eps = np.finfo(np.double).eps
        if (b == null) {
            b = a;
        }
        //a_prime = np.sign(a) / (1 - np.abs(a) + eps) - np.sign(a)
        for (int i = 0; i < a.length; i++) {
            final double A = a[i];
            final double asign = Math.signum(A);
            b[i] = asign / (1 - Math.abs(A) + eps) - asign;
        }
        return b;
    }

    public static double[] mapInfToOne(final double[] b, double[] a) {
        //""" Map values from [0, inf) onto [0, 1] and map values from  (-inf, 0] onto [-1, 0] """
        if (a == null) {
            a = b;
        }

        //a = np.sign(a_prime) * (1 - 1 / (np.abs(a_prime) + 1))
        for (int i = 0; i < a.length; i++) {
            final double B = b[i];
            final double bsign = Math.signum(B);
            a[i] = bsign * (1.0 - 1.0 / (Math.abs(B) + 1));
        }
        return a;
    }

    public static DenseMatrix pad(DenseMatrix a, int rows, int cols, double defaultValue) {
        /*
         Pad a matrix to the specified shape

         If any element of shape is 0, that size remains unchanged in 
         that axis. If any element of shape is < 0, the size in that
         axis is incremented by the magnitude of that value.
         Use val (default 0) to fill in the extra spaces. 
         */

        if (rows <= 0) {
            rows = a.numRows() - rows;
        }
        if (cols <= 0) {
            cols = a.numColumns() - cols;
        }

        if (rows < a.numRows()) {
            throw new RuntimeException("Padding with fewer rows: " + " " + a.numRows() + "->" + rows);
        }
        if (cols < a.numColumns()) {
            throw new RuntimeException("Padding with fewer cols: " + " " + a.numColumns() + "->" + cols);
        }

        DenseMatrix b = new DenseMatrix(rows, cols);
        if (defaultValue!=0.0)
            fill(b, defaultValue);
        
        //SubmatrixOps.setSubMatrix(a, b, 0, 0, 0, 0, a.numRows(), a.numColumns());
        setSubMatrix(a, b);
        return b;
    }
    
    public static DenseMatrix setSubMatrix(DenseMatrix small, DenseMatrix large) {        
        for (int i = 0; i < small.numRows(); i++)
            for (int j = 0; j < small.numColumns(); j++)
                large.set(i, j, small.get(i, j));
        return large;
        
    }
    
    public static DenseMatrix fill(DenseMatrix m, double v) {
        Arrays.fill(m.getData(), v);
        return m;
    }

    public static DenseVector boundedRowSum(DenseMatrix m) {
        DenseMatrix n = m.copy();
        mapOneToInf(n.getData(), null);
        
        DenseVector r = sumCols(n, null);

        mapInfToOne(r.getData(), null);
        return r;
    }

    static double[] boundedSum(final int axis, final double[]  
        ... a) {
        /* 
            Sum elements nonlinearly, such that the total is less than 1 

            To be more precise, as long as all elements in a are between -1
            and 1, their sum will also be between -1 and 1. a can be a 
            list or a numpy array. 
        */                
        final int size = a[0].length;
        final double[] total = new double[size];
        final double[] y = new double[size];
        for (final double[] x : a) {
            assert (size == x.length);

            double[] m = mapOneToInf(x, y);
            for (int i = 0; i < size; i++) {
                total[i] += m[i];
            }
        }
        return mapInfToOne(total, null);
        /*
         def bounded_sum(a, axis=0):
            
         if type(a) is list:
         total = map_one_to_inf(a[0])
         for item in a[1:]:
         total += map_one_to_inf(item)
         return map_inf_to_one(total)
         else:
         # handle the case where a is a one-dimensional array
         if len(a.shape) == 1:
         a = a[:, np.newaxis]
         bounded_total = map_inf_to_one(np.sum(map_one_to_inf(a), axis=axis))
         return bounded_total[:,np.newaxis]
        
         */
    }

    static DenseVector extractBooleanized(DenseMatrix cableActivities, DenseVector indexProjection) {
        //exracts the elements where indexProjection!=0
        // emulates: cog_cable_activities = self.cable_activities[self.ziptie.get_index_projection(cog_index).ravel().astype(bool)]

        indexProjection = indexProjection; //both row=oriented

        double[] activities = new double[indexProjection.size()];
        int activityNum = 0;
        for (int i = 0; i < indexProjection.size(); i++) {
            double v = indexProjection.get(i);
            if (v != 0) {
                activities[activityNum++] = cableActivities.get(i, 0);
            }
        }

        //DenseMatrix result = DenseMatrix.wrap(activityNum, 1, activities);
        return new DenseVector(activities);
    }

    static NumberFormat nf = new DecimalFormat("0.0000000");

    public static void printArray(double[] d) {
        for (double x : d) {
            System.out.print(nf.format(x) + " ");
        }
        System.out.println();
    }
    /*
     //INCOMPLETE BUT POSSIBLY UNNECESSARY
     static double[][] getNonZero(DenseMatrix x) {
     //returned first element: rows where it has non-zero elements
     //returned second element: cols where it has non-zero elements
     boolean[] nonemptyRows = new boolean[x.numRows()];
     for (int i = 0; i < x.numRows(); i++) {
     boolean zero = true;
     for (int j = 0; j < x.numColumns(); j++) {
     if (x.get(i, j)!=0) { zero = false; break; }
     }
     if (!zero)
     nonemptyRows[i] = true;
     }
        
     boolean[] nonemptyCols = new boolean[bundleMap.numColumns()];
        
     return new double[][] { nonemptyRows, nonemptyCols };        
     }
     */

    /**
     * zeros the resulting matrix and puts 1 where the input matrix is non-zero
     */
    public static DenseMatrix getNonZeroMask(final DenseMatrix x) {
        final DenseMatrix y = new DenseMatrix(x.numColumns(), x.numRows());
        final double[] xd = x.getData();
        final double[] yd = y.getData();
        for (int i = 0; i < yd.length; i++) {
            if (xd[i] != 0) {
                yd[i] = 1;
            }
        }
        return y;
    }

    /**
     * modifies the parameter, and returns it
     */
    public static DenseMatrix matrixBooleanize(final DenseMatrix x) {
        final double[] d = x.getData();
        for (int j = 0; j < d.length; j++) {
            d[j] = (d[j] != 0 ? 1.0 : 0.0);
        }
        return x;
    }

    public static void matrixSign(final DenseMatrix x) {
        final double[] d = x.getData();
        for (int j = 0; j < d.length; j++) {
            d[j] = Math.signum(d[j]);
        }
    }

    public static void matrixMaximum(final DenseMatrix x, final double v) {
        final double[] d = x.getData();
        for (int j = 0; j < d.length; j++) {
            d[j] = Math.max(d[j], v);
        }
    }

    public static void matrixMinimum(final DenseMatrix x, final double v) {
        final double[] d = x.getData();
        for (int j = 0; j < d.length; j++) {
            d[j] = Math.min(d[j], v);
        }
    }

    public static void matrixAbs(final DenseMatrix x) {
        final double[] d = x.getData();
        for (int j = 0; j < d.length; j++) {
            d[j] = Math.abs(d[j]);
        }
    }

    public static void matrixDivBy(final DenseMatrix x, final DenseMatrix numerator) {
        final double[] d = x.getData();
        final double[] n = numerator.getData();
        for (int j = 0; j < d.length; j++) {
            d[j] = n[j] / d[j];
        }
    }

    public static void matrixPower(final DenseMatrix m, final double exponent) {
        final double[] d = m.getData();
        for (int i = 0; i < d.length; i++) {
            if (exponent == -1) {
                d[i] = 1.0 / d[i]; //should be faster than Math.pow
            } else {
                d[i] = Math.pow(d[i], exponent);
            }
        }
    }

    public static void matrixPowerExp(final DenseMatrix m, final double base) {
        final double[] d = m.getData();
        for (int i = 0; i < d.length; i++) {
            d[i] = Math.pow(base, d[i]);
        }
    }

    
    static DenseMatrix maxCol(DenseMatrix x) {        
        final DenseMatrix projection = new DenseMatrix(x.numRows(), 1);
        final double[] pd = projection.getData();
        for (int i = 0; i < x.numColumns(); i++) {
            for (int j = 0; j < x.numRows(); j++) {
                if (i == 0) {
                    projection.set(j, 0, x.get(j, 0));
                } else {
                    double cg = x.get(j, i);
                    if (cg < pd[j]) {
                        projection.set(j, 0, cg);
                    }
                }
            }
        }
        return projection;
    }

    //TODO unify these

    static DenseMatrix minCol(DenseMatrix x) {
        final DenseMatrix projection = new DenseMatrix(x.numRows(), 1);
        final double[] pd = projection.getData();
        for (int i = 0; i < x.numColumns(); i++) {
            for (int j = 0; j < x.numRows(); j++) {
                if (i == 0) {
                    projection.set(j, 0, x.get(j, 0));
                } else {
                    double cg = x.get(j, i);
                    if (cg < pd[j]) {
                        projection.set(j, 0, cg);
                    }
                }
            }
        }
        return projection;
    }

    static DenseMatrix maxRow(final DenseMatrix x) {
        final DenseMatrix projection = new DenseMatrix(1, x.numColumns());
        final double[] pd = projection.getData();
        for (int i = 0; i < x.numRows(); i++) {
            for (int j = 0; j < x.numColumns(); j++) {
                if (i == 0) {
                    pd[j] = x.get(0, j);
                } else {
                    final double cg = x.get(i, j);
                    if (cg > pd[j]) {
                        pd[j] = cg;
                    }
                }
            }
        }
        return projection;
    }

    //TODO unify these two funcs ^v

    static DenseMatrix minRow(final DenseMatrix x) {
        final DenseMatrix projection = new DenseMatrix(1, x.numColumns());
        final double[] pd = projection.getData();
        for (int i = 0; i < x.numRows(); i++) {
            for (int j = 0; j < x.numColumns(); j++) {
                if (i == 0) {
                    pd[j] = x.get(0, j);
                } else {
                    final double cg = x.get(i, j);
                    if (cg < pd[j]) {
                        pd[j] = cg;
                    }
                }
            }
        }
        return projection;
    }

    //DEPRECATED use elemDiv or something
    static DenseMatrix matrixDivide(final DenseMatrix a, final DenseMatrix b) {
        DenseMatrix r = new DenseMatrix(a.numRows(), a.numColumns());
        final double[] ad = a.getData();
        final double[] bd = b.getData();
        final double[] rd = r.getData();
        assert (ad.length == bd.length);
        final double[] d = new double[ad.length];
        for (int i = 0; i < d.length; i++) {
            rd[i] = ad[i] / bd[i];
        }
        return r;
    }

    static DenseMatrix broadcastRows(final DenseMatrix col, final int numColumns) {
        //TODO see if can be created by arraycopy repeatedly, depends on row ordering        
        assert (col.numColumns() == 1);
        final double[] cd = col.getData();
        int numRows = col.numRows();
        final DenseMatrix r = new DenseMatrix(numRows, numColumns);
        for (int i = 0; i < numRows; i++) {
            final double C = cd[i];
            for (int j = 0; j < numColumns; j++) {
                r.set(i, j, C);
            }
        }
        return r;
    }

    static DenseMatrix broadcastCols(final DenseMatrix row, final int numRows) {
        //TODO see if can be created by arraycopy repeatedly, depends on row ordering        
        int numColumns = row.numColumns();
        final DenseMatrix r = new DenseMatrix(numRows, numColumns);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numColumns; j++) {
                r.set(i, j, row.get(0, j));
            }
        }
        return r;
    }

    //numpy.random.normal(loc=0.0, scale=1.0, size=None)
    //Extremely fast, medium quality randomness: http://www.javamex.com/tutorials/random_numbers/java_util_random_subclassing.shtml
    public static class XORShiftRandom extends Random {

        private long seed = System.nanoTime();

        public XORShiftRandom() {
        }

        @Override
        protected int next(int nbits) {
            // N.B. Not thread-safe!
            long x = this.seed;
            x ^= (x << 21);
            x ^= (x >>> 35);
            x ^= (x << 4);
            this.seed = x;
            x &= ((1L << nbits) - 1);
            return (int) x;
        }
    }

    private static final Random random = new XORShiftRandom();

    public static DenseMatrix normRandMatrix(int numRows, int numColumns, double scale, double min) {
        //TODO use normal distribution
        final DenseMatrix r = new DenseMatrix(numRows, numColumns);
        normRand(r, scale, min);
        return r;
    }

    public static void matrixAddNoise(DenseMatrix m, double scale) {
        final double[] d = m.getData();
        for (int i = 0; i < d.length; i++) {
            d[i] += random.nextDouble() * scale;
        }
    }

    public static void normRand(DenseMatrix r, double scale, double min) {
        final double[] d = r.getData();
        for (int i = 0; i < d.length; i++) {
            d[i] = random.nextGaussian() * scale + min;
        }

    }

    public static void setSinusoidal(DenseMatrix m, int col, double t, double baseFreq, double phase) {
        setSinusoidal(m, col, t, baseFreq, phase, 0.5, 0.5);
    }

    public static void setSinusoidal(DenseMatrix m, int col, double t, double baseFreq, double phase, double scale, double offset) {
        for (int i = 0; i < m.numRows(); i++) {
            m.set(i, col, Math.sin(phase + t / ((double) (1 + i)) / 3.14159 * baseFreq) * scale + offset);
        }
    }

    public static void addNoise(DenseMatrix m, double NOISE_FACTOR) {
        double[] d = m.getData();
        for (int i = 0; i < d.length; i++) {
            d[i] += Math.random() * NOISE_FACTOR;
        }
    }

    public static DenseVector getColumn(DenseMatrix m, int j) {
        DenseVector v = new DenseVector(m.numRows());
        //May be optimized by array copy if columns stored contiguously
        for (int i = 0; i < v.size(); i++) {
            v.set(i, m.get(i, j));
        }
        return v;
    }


    public static DenseVector getRow(DenseMatrix m, int j) {
        DenseVector v = new DenseVector(m.numColumns());
        //May be optimized by array copy if rows stored contiguously
        for (int i = 0; i < v.size(); i++) {
            v.set(i, m.get(j, i));
        }
        return v;
    }

    public static DenseMatrix getColumns(DenseMatrix matrix, int[] cols) {

        int[] rowIndices = new int[matrix.numRows()];
        for (int i = 0; i < rowIndices.length; i++) {
            rowIndices[i] = i;
        }

        return (DenseMatrix) (Matrices.getSubMatrix(matrix, rowIndices,cols));
    }

    public static DenseMatrix getRows(DenseMatrix matrix, int[] rows) {

        int[] columnIndices = new int[matrix.numColumns()];
        for (int i = 0; i < columnIndices.length; i++) {
            columnIndices[i] = i;
        }

        return (DenseMatrix) (Matrices.getSubMatrix(matrix, rows, columnIndices));
    }    
    
    public static void main(String[] args) {
        
    }
    
    public static class ClassScope {

        private static java.lang.reflect.Field LIBRARIES;



        public static java.util.Vector getLoadedLibraries(final ClassLoader loader) throws Exception {
            LIBRARIES = ClassLoader.class.getDeclaredField("loadedLibraryNames");
            LIBRARIES.setAccessible(true);
            
            final java.util.Vector libraries = (java.util.Vector)LIBRARIES.get(loader);
            return libraries;
        }
    }

    static {

        
        try {
            
            //Class javaBlasClass = Class.forName("org.netlib.blas.JBLAS")
            
            //Field javaBlas = javaBlasClass.getDeclaredField("INSTANCE");
            //Field jInstance = javaBlas.getClass().getDeclaredField("modifiers");
            //jInstance.setAccessible(true);
            //jInstance.setInt(javaBlas,javaBlas.getModifiers() & ~Modifier.FINAL);
            //javaBlas.setAccessible(true);
                    
            //get native blas object and make it accessible
            
            
            /*
            Class nativeBlasClass = Class.forName("org.netlib.blas.NativeBLAS");
            Field nativeBlas = nativeBlasClass.getDeclaredField("INSTANCE");
            Field nInstance = nativeBlas.getClass().getDeclaredField("modifiers");
            nInstance.setAccessible(true);
            nInstance.setInt(nativeBlas,nativeBlas.getModifiers() & ~Modifier.FINAL);
            nativeBlas.setAccessible(true);
            
            //get blas current object and make it accessible
            Field blasCurrent = Class.forName("org.netlib.blas.BLAS").getDeclaredField("current");
            Field bInstance = blasCurrent.getClass().getDeclaredField("modifiers");
            bInstance.setAccessible(true);
            bInstance.setInt(blasCurrent, blasCurrent.getModifiers() & ~Modifier.FINAL);
            blasCurrent.setAccessible(true);
            
            //SET TO NativeBLAS
            blasCurrent.set(null, nativeBlas.get(null));
            
            
            //SET TO JBLAS
            //blasCurrent.set(null, javaBlas.get(null));
            */
        }
        catch (Exception e) {
            System.err.println(e);
        }

        /*
        try {
            final java.util.Vector<String> libraries = ClassScope.getLoadedLibraries(ClassLoader.getSystemClassLoader()); //MyClassName.class.getClassLoader()        
            
            System.out.println(libraries);
        } catch (Exception ex) {
            Logger.getLogger(MTJTest.class.getName()).log(Level.SEVERE, null, ex);
        } 
        */
    }
    
}
