package becca.core;

import org.apache.commons.math3.util.Precision;
import org.ejml.alg.dense.mult.SubmatrixOps;
import org.ejml.data.BlockMatrix64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.data.ReshapeMatrix64F;
import org.ejml.ops.CommonOps;

public class Util extends CommonOps {
    public final static double EPSILON = Precision.EPSILON;   // // sys.float_info.epsilon = 2.220446049250313e-16
    public final static double BIG = Math.pow(10, 20);
    public final static int MAX_INT16 = 32767;     // np.iinfo(np.int16).max

    static int log(int x, int base) {
        return (int) (Math.log(x) / Math.log(base));
    }

    public static String m(ReshapeMatrix64F m) {
        return "[" + m.getNumRows() + "," + m.getNumCols() + "]";
    }
    
    public static DenseMatrix64F multMatrixMatrix(DenseMatrix64F a, DenseMatrix64F b) {
        assert(a.getNumRows() == b.getNumRows());
        assert(a.getNumCols() == b.getNumCols());
        
        DenseMatrix64F c = new DenseMatrix64F(a.getNumRows(), a.getNumCols());
        double[] cd = c.getData();
        double[] ad = a.getData();
        double[] bd = b.getData();
        for (int i = 0; i < cd.length; i++) {
            cd[i] = ad[i] * bd[i];
        }
        
        return c;        
    }
    
    public static DenseMatrix64F getWeightedAverage(DenseMatrix64F values, DenseMatrix64F weights) {
        //""" Perform a weighted average of values, using weights """
        
        if (values.getNumCols() < weights.getNumCols()) {
            values = broadcastRows(values, weights.getNumCols());
        }
        else if (values.getNumCols() > weights.getNumCols()) {
            weights = broadcastRows(weights, values.getNumCols());
        }
        
        //weighted_sum_values = np.sum(values * weights, axis=0)                 
        DenseMatrix64F valueWeightProduct = multMatrixMatrix(values, weights);        
        
        DenseMatrix64F weightedSumValues = new DenseMatrix64F(valueWeightProduct.getNumCols(), 1);
        sumCols(valueWeightProduct, weightedSumValues);        
        
        //sum_of_weights = np.sum(weights, axis=0)         
        DenseMatrix64F sumOfWeights = sumCols(weights, null);        
        
        //return (weighted_sum_values / (sum_of_weights + EPSILON))[:,np.newaxis]
        add(sumOfWeights, EPSILON);
        
        final double[] wsd = weightedSumValues.getData();
        final double[] sowd = sumOfWeights.getData();
        for (int i = 0; i < wsd.length; i++)
            wsd[i] /= sowd[i];
        return weightedSumValues;        
    }

    
    public static void printMatrixDimensions(ReshapeMatrix64F... m) {
        String s = "";
        for (ReshapeMatrix64F r : m) 
            s += m(r) + " ";
        System.out.println("matrixDim: " + s);
    }
    
    public static DenseMatrix64F matrixVector(final DenseMatrix64F matrix, final DenseMatrix64F vector) {
        return matrixVector(matrix, vector, true);
    }
    public static DenseMatrix64F matrixVector(final DenseMatrix64F matrix, final DenseMatrix64F vector, final boolean multiply) {
        //ex: (8, 32) * (1, 32) -> (8, 32)
        
        final DenseMatrix64F result = new DenseMatrix64F(matrix.getNumRows(), matrix.getNumCols());
        final double[] vdata = vector.getData();
        
        if (vector.getNumRows() == 1) {
            assert(matrix.getNumCols() == vector.getNumCols());
            
            for (int i = 0; i < matrix.getNumRows(); i++) {
               for (int j = 0; j < matrix.getNumCols(); j++) {
                   if (multiply)
                       result.set(i, j, matrix.get(i, j) * vector.get(0, j));
                   else
                       result.set(i, j, matrix.get(i, j) / vector.get(0, j));
               }
            }
        }
        else {
            assert(matrix.getNumRows() == vector.getNumRows());

            for (int i = 0; i < matrix.getNumRows(); i++) {
               for (int j = 0; j < matrix.getNumCols(); j++) {
                   if (multiply)
                       result.set(i, j, matrix.get(i, j) * vector.get(i, 0));
                   else
                       result.set(i, j, matrix.get(i, j) / vector.get(i, 0));
               }
            }

        }
        return result;        
    }
    
    public static double dot(DenseMatrix64F a, DenseMatrix64F b) {
        double[] ad = a.getData();
        double[] bd = b.getData();
        assert(ad.length == bd.length);
        double sum = 0;
        for (int i = 0; i < ad.length; i++)
            sum += ad[i] * bd[i];        
        return sum;        
    }
    
    public static DenseMatrix64F getGeneralizedMean(DenseMatrix64F values, DenseMatrix64F weights, double exponent) {
        DenseMatrix64F shiftedValues = values.copy();
        add(shiftedValues, 1.0);

        DenseMatrix64F valuesToPower = shiftedValues.copy();
        matrixPower(valuesToPower, exponent);


        //mean_values_to_power = weighted_average(values_to_power, weights)
        DenseMatrix64F meanValuesToPower = getWeightedAverage(valuesToPower, weights);
        
        //shifted_mean = (mean_values_to_power + EPSILON) ** (1./exponent)
        DenseMatrix64F shiftedMean = meanValuesToPower;
        add(shiftedMean, EPSILON);
        matrixPower(shiftedMean, (1.0/exponent));
        
        
        DenseMatrix64F mean = shiftedMean;
        
        //mean = shifted_mean - 1.
        add(mean, -1);
                        
        //# Find means for which all weights are zero. These are undefined.
        //# Set them equal to zero.
        //sum_weights = np.sum(weights, axis=0)
        DenseMatrix64F sumWeights = sumCols(weights, null);
        
                
        //zero_indices = np.where(np.abs(sum_weights) < EPSILON)
        //mean[zero_indices] = 0.        
        double[] meanD = mean.getData();
        double[] sumWeightsD = sumWeights.getData();
        for (int i = 0; i < sumWeights.getNumRows(); i++) {
            double a = Math.abs(sumWeightsD[i]);
            if (a < EPSILON)
                meanD[i] = 0;
        }
        
        return mean;        
    }

    public static double[] mapOneToInf(final double[] a) {
        //""" Map values from [0, 1] onto [0, inf) and map values from [-1, 0] onto (-inf, 0]
        final double eps = 2.2204460492503131e-16; //eps = np.finfo(np.double).eps
        final double[] b = new double[a.length];
        //a_prime = np.sign(a) / (1 - np.abs(a) + eps) - np.sign(a)
        for (int i = 0; i < a.length; i++) {
            final double A = a[i];
            final double asign = (A > 0 ? +1 :((A < 0) ? -1 : 0));
            b[i] = asign / (1 - Math.abs(A) + eps) - asign;
        }
        return b;
    }
    public static double[] mapInfToOne(final double[] b) {
        //""" Map values from [0, inf) onto [0, 1] and map values from  (-inf, 0] onto [-1, 0] """
        final double eps = 2.2204460492503131e-16; //eps = np.finfo(np.double).eps
        final double[] a = new double[b.length];
        //a = np.sign(a_prime) * (1 - 1 / (np.abs(a_prime) + 1))
        for (int i = 0; i < a.length; i++) {
            final double B = b[i];
            final double bsign = (B > 0 ? +1 :((B < 0) ? -1 : 0));
            a[i] = bsign * (1 - 1.0 / (Math.abs(B) + 1));
        }
        return a;        
    }

    
    public static DenseMatrix64F pad(DenseMatrix64F a, int rows, int cols, double defaultValue) {
        /*
        Pad a matrix to the specified shape

        If any element of shape is 0, that size remains unchanged in 
        that axis. If any element of shape is < 0, the size in that
        axis is incremented by the magnitude of that value.
        Use val (default 0) to fill in the extra spaces. 
        */

        if (rows <= 0)
            rows = a.getNumRows() - rows;
        if (cols <= 0)
            cols = a.getNumCols() - cols;
        
        if (rows < a.getNumRows())
            throw new RuntimeException("Padding with fewer rows: " + " "+ a.getNumRows() + "->" + rows);
        if (cols < a.getNumCols())
            throw new RuntimeException("Padding with fewer cols: " + " "+ a.getNumCols() + "->" + cols);
        
        DenseMatrix64F b = new DenseMatrix64F(rows, cols);
        fill(b, defaultValue);
        SubmatrixOps.setSubMatrix(a, b, 0, 0, 0, 0, a.getNumRows(), a.getNumCols());
        return b;
    }

    public static DenseMatrix64F boundedRowSum(DenseMatrix64F m) {
        DenseMatrix64F[] rows = new DenseMatrix64F[m.getNumRows()];
        rowsToVector(m, rows);
        double[][] drows = new double[m.getNumRows()][];
        for (int i = 0; i < rows.length; i++)
            drows[i] = rows[i].getData();
        
        double[] result = boundedSum(0, drows);
        return DenseMatrix64F.wrap(result.length, 1, result);
    }

    static double[] boundedSum(int axis, double[]... a) {
        /* 
            Sum elements nonlinearly, such that the total is less than 1 

            To be more precise, as long as all elements in a are between -1
            and 1, their sum will also be between -1 and 1. a can be a 
            list or a numpy array. 
        */
        
        //assume all 'a' have same size
        
        BlockMatrix64F total = new BlockMatrix64F(a[0].length,1);
        int size = -1;
        for (double[] x : a) {
            
            if (size == -1) size = x.length;
            else assert(size==x.length);
            
            double[] m = mapOneToInf(x);
            BlockMatrix64F  ix = BlockMatrix64F.wrap(m, size, 1, size);
            //System.out.println(size + " " + total.getNumRows() +" " +ix.getNumRows() + " " + m.length + " " + x.length);
            
            addEquals(total, ix);
            
        }
        return mapInfToOne(total.getData());
        
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


    static DenseMatrix64F extractBooleanized(DenseMatrix64F cableActivities, DenseMatrix64F indexProjection) {
        //exracts the elements where indexProjection!=0
        // emulates: cog_cable_activities = self.cable_activities[self.ziptie.get_index_projection(cog_index).ravel().astype(bool)]
        
        indexProjection = transpose(indexProjection, null); //both row=oriented
        
        
        
        double[] activities = new double[indexProjection.getNumRows()];
        int activityNum = 0;
        for (int i = 0; i < indexProjection.getNumRows(); i++) {
            double v = indexProjection.get(i, 0);
            if (v!=0) {
                activities[activityNum++] = cableActivities.get(i, 0);
            }
        }
        
        DenseMatrix64F result = DenseMatrix64F.wrap(activityNum, 1, activities);
        
        /*if (activityNum > 0) {
            System.out.println("cable projection");
            System.out.println(transpose(cableActivities, null));
            System.out.println(transpose(indexProjection, null));            
            System.out.println(transpose(result, null));        
            System.out.println();
            System.out.println();
        }*/

        return result;
    }

    /*
    //INCOMPLETE BUT POSSIBLY UNNECESSARY
    static double[][] getNonZero(DenseMatrix64F x) {
        //returned first element: rows where it has non-zero elements
        //returned second element: cols where it has non-zero elements
        boolean[] nonemptyRows = new boolean[x.getNumRows()];
        for (int i = 0; i < x.getNumRows(); i++) {
            boolean zero = true;
            for (int j = 0; j < x.getNumCols(); j++) {
                if (x.get(i, j)!=0) { zero = false; break; }
            }
            if (!zero)
                nonemptyRows[i] = true;
        }
        
        boolean[] nonemptyCols = new boolean[bundleMap.getNumCols()];
        
        return new double[][] { nonemptyRows, nonemptyCols };        
    }
    */

    public static DenseMatrix64F getNonZeroMask(final DenseMatrix64F x) {
        final DenseMatrix64F y = x.copy();
        final double[] yd = y.getData();
        for (int i = 0; i < yd.length; i++)
            if (yd[i]!=0) yd[i] = 1;
        return y;
    }
    
    /** modifies the parameter, and returns it */
    public static DenseMatrix64F matrixBooleanize(final DenseMatrix64F x) {
        final double[] d = x.getData();
        for (int j = 0; j < d.length; j++)
            d[j] = (d[j]!=0 ? 1.0 : 0.0);
        return x;
    }
    public static void matrixSign(final DenseMatrix64F x) {
        final double[] d = x.getData();
        for (int j = 0; j < d.length; j++)
            d[j] = Math.signum(d[j]); 
    }
    public static void matrixMaximum(final DenseMatrix64F x, final double v) {
        final double[] d = x.getData();
        for (int j = 0; j < d.length; j++)
            d[j] = Math.max(d[j], v); 
    }    
    public static void matrixMinimum(final DenseMatrix64F x, final double v) {
        final double[] d = x.getData();
        for (int j = 0; j < d.length; j++)
            d[j] = Math.min(d[j], v); 
    }    
    public static void matrixAbs(final DenseMatrix64F x) {
        final double[] d = x.getData();
        for (int j = 0; j < d.length; j++)
            d[j] = Math.abs(d[j]); 
    }    
    public static void matrixDivBy(final DenseMatrix64F x, final DenseMatrix64F numerator) {
        final double[] d = x.getData();
        final double[] n = numerator.getData();
        for (int j = 0; j < d.length; j++)
            d[j] = n[j] / d[j];
    }    
    
    public static void matrixPower(final DenseMatrix64F m, final double exponent) {        
        final double[] d = m.getData();
        for (int i = 0; i < d.length; i++) {
            if (exponent == -1)
                d[i] = 1.0 / d[i]; //should be faster than Math.pow
            else
                d[i] = Math.pow(d[i], exponent);
        }
    }

    static DenseMatrix64F maxCol(DenseMatrix64F x) {
        final DenseMatrix64F projection = new DenseMatrix64F(x.getNumRows(), 1);
        for (int i = 0; i < x.getNumCols(); i++) {
            for (int j = 0; j < x.getNumRows(); j++) {
                if (i == 0) {
                    projection.set(j, 0, x.get(j, 0));
                }
                else {
                    double cg = x.get(j,i);
                    if (cg > projection.get(j,0))
                        projection.set(j, 0, cg);
                }
            }
        }
        return projection;
    }
    
    static DenseMatrix64F maxRow(DenseMatrix64F x) {
        final DenseMatrix64F projection = new DenseMatrix64F(1, x.getNumCols());
        for (int i = 0; i < x.getNumRows(); i++) {
            for (int j = 0; j < x.getNumCols(); j++) {
                if (i == 0) {
                    projection.set(0, j, x.get(0, j));
                }
                else {
                    double cg = x.get(i, j);
                    if (cg > projection.get(0, j))
                        projection.set(0, j, cg);
                }
            }
        }
        return projection;
    }
    //TODO unify these two funcs ^v
    static DenseMatrix64F minRow(DenseMatrix64F x) {
        final DenseMatrix64F projection = new DenseMatrix64F(1, x.getNumCols());
        for (int i = 0; i < x.getNumRows(); i++) {
            for (int j = 0; j < x.getNumCols(); j++) {
                if (i == 0) {
                    projection.set(0, j, x.get(0, j));
                }
                else {
                    double cg = x.get(i, j);
                    if (cg < projection.get(0, j))
                        projection.set(0, j, cg);
                }
            }
        }
        return projection;
    }
    
    
    
    //DEPRECATED use elemDiv or something
    static DenseMatrix64F matrixDivide(final DenseMatrix64F a, final DenseMatrix64F b) {
        final double[] ad = a.getData();
        final double[] bd = b.getData();
        assert(ad.length == bd.length);
        final double[] d = new double[ad.length];
        for (int i = 0; i < d.length; i++)
            d[i] = ad[i] / bd[i];
        return DenseMatrix64F.wrap(ad.length, 1, d);
    }

    static DenseMatrix64F broadcastRows(final DenseMatrix64F col, final int numCols) {
        //TODO see if can be created by arraycopy repeatedly, depends on row ordering        
        int numRows = col.getNumRows();
        final DenseMatrix64F r = new DenseMatrix64F(numRows, numCols);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                r.set(i, j, col.get(i, 0));
            }
        }
        return r;
    }

    static DenseMatrix64F broadcastCols(final DenseMatrix64F row, final int numRows) {
        //TODO see if can be created by arraycopy repeatedly, depends on row ordering        
        int numCols = row.getNumCols();
        final DenseMatrix64F r = new DenseMatrix64F(numRows, numCols);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                r.set(i, j, row.get(0, j));
            }
        }
        return r;
    }

    public static DenseMatrix64F randMatrix(int numRows, int numCols, double range) {
        //TODO use normal distribution
        final DenseMatrix64F r = new DenseMatrix64F(numRows, numCols);
        final double[] d = r.getData();
        for (int i = 0; i < d.length; i++)
            d[i] = Math.random() * range;
        return r;
    }
    
    
}
