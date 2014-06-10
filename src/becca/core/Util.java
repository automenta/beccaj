package becca.core;

import org.apache.commons.math3.util.Precision;
import org.encog.mathutil.matrices.Matrix;

public class Util {
    public final static double EPSILON = Precision.EPSILON;   // // sys.float_info.epsilon = 2.220446049250313e-16
    public final static double BIG = Math.pow(10, 20);
    public final static int MAX_INT16 = 32767;     // np.iinfo(np.int16).max

    static int log(int x, int base) {
        return (int) (Math.log(x) / Math.log(base));
    }

    public static Matrix pad(Matrix m, int rows, int cols, double value) {
        return m;
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

    static Matrix boundedSum(int axis, Matrix... a) {
        /*
        def bounded_sum(a, axis=0):
            """ 
            Sum elements nonlinearly, such that the total is less than 1 

            To be more precise, as long as all elements in a are between -1
            and 1, their sum will also be between -1 and 1. a can be a 
            list or a numpy array. 
            """ 
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
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
