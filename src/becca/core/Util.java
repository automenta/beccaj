package becca.core;

import org.apache.commons.math3.util.Precision;
import org.ejml.alg.dense.mult.SubmatrixOps;
import org.ejml.data.BlockMatrix64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

public class Util {
    public final static double EPSILON = Precision.EPSILON;   // // sys.float_info.epsilon = 2.220446049250313e-16
    public final static double BIG = Math.pow(10, 20);
    public final static int MAX_INT16 = 32767;     // np.iinfo(np.int16).max

    static int log(int x, int base) {
        return (int) (Math.log(x) / Math.log(base));
    }

    /*
    def pad(a, shape, val=0.):
    """
    Pad a numpy array to the specified shape
    
    If any element of shape is 0, that size remains unchanged in 
    that axis. If any element of shape is < 0, the size in that
    axis is incremented by the magnitude of that value.
    Use val (default 0) to fill in the extra spaces. 
    """
    if shape[0] <= 0:
        rows = a.shape[0] - shape[0]
    else:
        rows = shape[0]
        # debug
        if rows < a.shape[0]:
            print ' '.join(['a.shape[0] is', str(a.shape[0]), ' but trying to',
                            ' pad to ', str(rows), 'rows.'])
    if shape[1] <= 0:
        cols = a.shape[1] - shape[1]
    else:
        cols = shape[1]
        # debug
        if cols < a.shape[1]:
            print ' '.join(['a.shape[1] is', str(a.shape[1]), ' but trying to',
                            ' pad to ', str(cols), 'cols.'])
    padded = np.ones((rows,cols)) * val
    padded[:a.shape[0], :a.shape[1]] = a
    return padded
    */

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
        CommonOps.fill(b, defaultValue);
        SubmatrixOps.setSubMatrix(a, b, 0, 0, 0, 0, a.getNumRows(), a.getNumCols());
        return b;
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
            
            BlockMatrix64F  ix = BlockMatrix64F.wrap(mapOneToInf(x), size, 1, size);
            System.out.println(size + " " + total.getNumElements() +" " +ix.getNumElements());
            CommonOps.addEquals(total, ix);
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

    static double sum(BlockMatrix64F s) {
        return CommonOps.elementSum(s);
    }
    
}
