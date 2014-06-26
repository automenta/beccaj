/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package becca.test.util;

import becca.core_mtj.DenseMatrix;
import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLException.MemObjectAllocationFailure;
import com.nativelibs4java.opencl.JavaCL;
import com.nativelibs4java.opencl.blas.ujmp.CLDenseFloatMatrix2D;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ujmp.core.MatrixFactory;
import static org.ujmp.core.calculation.Calculation.ORIG;
import org.ujmp.core.doublematrix.DoubleMatrix;
import org.ujmp.core.enums.ValueType;
import org.ujmp.core.mapper.MatrixMapper;

/**
 *
 * @author me
 */
public class MTJTest  {
    static final CLContext context = JavaCL.createBestContext();    
        
    public void testSolve() {
        Vector bv = Matrices.random(50000);
        Matrix am = Matrices.random(50000, 500);
        Vector xv = new DenseVector(am.numColumns());
        for (int x = 0; x < am.numColumns(); x++) {
            xv.set(x, 1);
        }
        xv = Matrices.random(xv.size());
        xv = am.solve(bv, xv);
        System.out.println(xv);        
    }
    
    public static void testMult(char op, int size) throws Exception {
        int runs = 32;
                
        
        System.gc();

        double totalTime3 = 0;
        
        try {
            org.ujmp.core.Matrix a = MatrixFactory.dense(ValueType.FLOAT, size, size);
            CLDenseFloatMatrix2D b = new CLDenseFloatMatrix2D(size, size);
            CLDenseFloatMatrix2D v = new CLDenseFloatMatrix2D(size, 1);

            a.zeros(ORIG);
            a.plus(ORIG, true, 1.0);

            b.zeros(ORIG);
            b.plus(ORIG, true, 2.0);

            v.zeros(ORIG);
            v.plus(ORIG, true, 2.0);


            for (int i = 0; i < runs; i++) {

                long start = System.nanoTime();

                if (op == '*') {
                    org.ujmp.core.Matrix c = a.times(b);
                    DoubleMatrix d = c.toDoubleMatrix(); 
                }
                else if (op == '+') {
                    org.ujmp.core.Matrix c = a.plus(b);
                    DoubleMatrix d = c.toDoubleMatrix();
                    /*System.out.println(d.getClass() + " " + d.getRowCount() + " " + d.getColumnCount());*/
                }
                else if (op == '|') {
                    org.ujmp.core.Matrix c = a.times(v);
                    DoubleMatrix d = c.toDoubleMatrix();
                }


                totalTime3 += System.nanoTime() - start;
            }
            totalTime3 = totalTime3 / runs / 1000000; //result in ms        
        }
        catch (MemObjectAllocationFailure me) {
            System.err.println(me);
        }
        
        
        
        boolean mtj = true; double totalTime = 0;
        boolean ejml = true; double totalTime2 = 0;  
        if (size > 4096) {
            mtj = false;
            ejml = false;
        }
        
        
        if (mtj) {
            //MTJ --------

            //        Matrix am = Matrices.random(size, size);
            //        Matrix bm = Matrices.random(size, size);
            //        Matrix cm = Matrices.random(size, size);
            DenseMatrix am = becca.core_mtj.Util.normRandMatrix(size, size, 1.0, 0.0);
            DenseMatrix bm = becca.core_mtj.Util.normRandMatrix(size, size, 1.0, 0.0);
            DenseMatrix cm = new DenseMatrix(size, size);
            /*DenseVector vm = new DenseVector(size);
            DenseVector um = new DenseVector(size);*/
            DenseMatrix vm = new DenseMatrix(size, 1);
            DenseMatrix um = new DenseMatrix(size, 1);


            for (int i = 0; i < runs; i++) {
                System.gc();

                long start = System.nanoTime();

                if (op == '*')
                    am.mult(bm, cm);
                else if (op == '+') {
                    am.addEquals(bm);
                    //am.multAdd(1.0, bm, cm); //SLOW
                }
                else if (op == '|') {
                    //am.mult(vm, um);
                    am.mult(vm, um);
                }

                totalTime += System.nanoTime() - start;
            }
            totalTime = totalTime / runs / 1000000; //result in ms

        }
        
        if (ejml) {

            DenseMatrix64F em = becca.core.Util.normRandMatrix(size, size, 1.0, 0.0);
            DenseMatrix64F fm = becca.core.Util.normRandMatrix(size, size, 1.0, 0.0);
            DenseMatrix64F gm = new DenseMatrix64F(size, size);               
            DenseMatrix64F wm = new DenseMatrix64F(size, 1);
            DenseMatrix64F xm = new DenseMatrix64F(size, 1);
         
            for (int i = 0; i < runs; i++) {
                System.gc();

                long start = System.nanoTime();

                if (op == '*')
                    CommonOps.mult(em, fm, gm);
                else if (op == '+')
                    CommonOps.addEquals(em, fm);
                else if (op == '|')
                    CommonOps.mult(em, wm, xm);

                totalTime2 += System.nanoTime() - start;
            }
            totalTime2 = totalTime2 / runs / 1000000; //result in ms
        }
        
        //System.out.println(cm);
        System.out.println(op + "_" + size + " , " + totalTime + " , " + totalTime2 + ", " + totalTime3);
    }

    public static void main(String[] args) throws Exception {
        //JavaCL + UJMP -----

        System.out.println("JavaCL Context: " + context);
        System.out.println("JavaCL Max Memory: " + context.getMaxMemAllocSize());
        
        MatrixMapper.getInstance().setDenseFloatMatrix2DClassName(CLDenseFloatMatrix2D.class.getName());
        

        testMult('+', 16);
        testMult('+', 32);
        testMult('+', 64);
        testMult('+', 128);
        testMult('+', 256);
        testMult('+', 300);
        testMult('+', 512);
        testMult('+', 768);
        testMult('+', 1024);
        testMult('+', 1500);
        testMult('+', 2048);
        testMult('+', 4096);
        testMult('+', 8192);
        testMult('+', 16384);
        testMult('+', 32767);

        testMult('|', 16);
        testMult('|', 32);
        testMult('|', 64);
        testMult('|', 128);
        testMult('|', 256);
        testMult('|', 300);
        testMult('|', 512);
        testMult('|', 768);
        testMult('|', 1024);
        testMult('|', 1500);
        testMult('|', 2048);
        testMult('|', 4096);
        testMult('|', 8192);
        
        testMult('*', 16);
        testMult('*', 32);
        testMult('*', 64);
        testMult('*', 128);
        testMult('*', 256);
        testMult('*', 512);
        testMult('*', 768);
        testMult('*', 1024);
        testMult('*', 1500);
        testMult('*', 2048);
        testMult('*', 4096);
        testMult('*', 8192);
    }
}
