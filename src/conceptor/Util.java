/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package conceptor;

import java.util.Random;
import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.ops.CommonOps;
import static org.ejml.ops.CommonOps.fill;

/**
 *
 * @author me
 */
public class Util extends becca.core.Util {
    
    public static DenseMatrix64F newInternalWeights(int n, double connectivity) {
        //    # Local Variables: disp, success, connectivity, nInternalUnits, internalWeights, specRad, opts
        //    # Function calls: generate_internal_weights, not, abs, sprandn, eigs
        //    #% Create a random sparse reservoir for an ESN. Nonzero weights are normal
        //    #% distributed.
        //    #%
        //    #% inputs:
        //    #% nInternalUnits = the number of internal units in the ESN
        //    #% connectivity: a real in [0,1], the (rough) proportion of nonzero weights
        //    #%
        //    #% output:
        //    #% internalWeights = matrix of size nInternalUnits x nInternalUnits
        //    success = False
        //    while not(success):
        //        #sprandn = Sparse normally distributed random matrix
        //        #R = sprandn(m,n,density) is a random, m-by-n, sparse matrix with approximately density*m*n normally distributed nonzero entries ((0 <= density <= 1).
        //        
        //        internalWeights = sprandn(nInternalUnits, nInternalUnits, connectivity)
        //        opts.disp = 0.
        //        specRad = np.abs(eigs(internalWeights, 1., 'lm', opts))
        //        internalWeights = matdiv(internalWeights, specRad)
        //        success = True
        //
        //                
        //        
        //        
        //    return [internalWeights]
        
        DenseMatrix64F d = new DenseMatrix64F(n, n);
        boolean finished = false;
        EigenDecomposition<DenseMatrix64F> eig = DecompositionFactory.eig(n, false);
        
        while (true) {
                        
            matrixSprandN(d, connectivity);

            //eigs(A,k,sigma) and eigs(A,B,k,sigma) return k eigenvalues based on sigma, which can take any of the following values:


            if (!eig.decompose(d)) {
                continue;
            }
            double max = 0;
            for (int i = 0; i < eig.getNumberOfEigenvalues(); i++) {
                final Complex64F ev = eig.getEigenvalue(i);
                if (i == 0)
                    max = ev.getMagnitude();
                else
                    max = Math.max(max, ev.getMagnitude());                    
            }                

            CommonOps.scale(1.0/max, d);
            break;

        }

        
        return d;
    }
    
    final static Random random = new Random();
    
    public static DenseMatrix64F matrixSprandN(DenseMatrix64F d, double sparsity) {
        fill(d, 0);
        int nr = d.getNumRows();
        int nc = d.getNumCols();
        int entries = (int)Math.round(nr * nc * sparsity);
        for (int i = 0; i < entries; i++) {
            int rr = (int)(Math.random()*nr);
            int rc = (int)(Math.random()*nc);
            d.set(rr, rc, random.nextGaussian());
        }
        return d;
    }

    static void normalizeColumn(DenseMatrix64F d, int c) {
        double min = 0, max = 0;
        
        for (int i = 0; i < d.numRows; i++) {
            final double v = d.get(i, c);
            if (i == 0) {
                min = max = v;
            }
            else {
                min = Math.min(min, v);
                max = Math.max(max, v);
            }
        }
        if (min==max) return;
        for (int i = 0; i < d.numRows; i++) {
            double v = d.get(i, c);
            v = (v - min) / (max - min);
            d.set(i, c, v);
        }
    }

    static void tanh(DenseMatrix64F x) {
        double[] d = x.data;
        for (int i = 0; i < x.elements; i++) {
            d[i] = Math.tanh(d[i]);
        }
            
    }



}
