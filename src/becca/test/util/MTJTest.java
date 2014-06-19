/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package becca.test.util;

import becca.core_mtj.DenseMatrix;
import becca.core_mtj.Util;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 *
 * @author me
 */
public class MTJTest extends Util {


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
    
    public static void testMult(char op, int size) {
        int runs = 16;
        
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
        
        double totalTime = 0;
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
        
        //System.out.println(cm);
        
        
        DenseMatrix64F em = becca.core.Util.normRandMatrix(size, size, 1.0, 0.0);
        DenseMatrix64F fm = becca.core.Util.normRandMatrix(size, size, 1.0, 0.0);
        DenseMatrix64F gm = new DenseMatrix64F(size, size);               
        DenseMatrix64F wm = new DenseMatrix64F(size, 1);
        DenseMatrix64F xm = new DenseMatrix64F(size, 1);
        
        double totalTime2 = 0;
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
        
        //System.out.println(cm);
        System.out.println(op + "_" + size + " , " + totalTime + " , " + totalTime2);
    }
    
    public static void main(String[] args) {
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
    }
}
