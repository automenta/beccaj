package becca.core_mtj;

import no.uib.cipr.matrix.Matrix;

/**
 * Faster DenseMatrix class without stupid check() calls.  Original author didnt need to make 'data' private so we work around it
 * 
 * TODO the same for DenseVector
 * 
 * @author me
 */
public class DenseMatrix extends no.uib.cipr.matrix.DenseMatrix {
    private final double[] d;

    public DenseMatrix(int numRows, int numColumns) {
        super(numRows, numColumns);        
        this.d = getData();
    }
    
    public DenseMatrix(Matrix A) {
        super(A);
        this.d = getData();        
    }
     
     @Override
    public DenseMatrix copy() {
        return new DenseMatrix(this);
    }

    public Matrix addEquals(DenseMatrix B) {        
        final double[] ad = this.getData();
        final double[] bd = B.getData();
        assert(ad.length == bd.length);
        
        for (int i = 0; i < ad.length; i++)
            ad[i] += bd[i];
        
        return this;
    }
    
    
   @Override
    public void add(int row, int column, double value) {
        d[row + column * numRows] += value;
    }

    @Override
    public void set(int row, int column, double value) {
        d[row + column * numRows] = value;
    }

    @Override
    public double get(int row, int column) {
        return d[row + column * numRows];
    }

    /**
* Checks the row and column indices, and returns the linear data index
*/
    int getIndex(int row, int column) {
        //check(row, column);
        return row + column * numRows;
    }    

    
    
}
