/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package javacl;


import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.JavaCL;
import com.nativelibs4java.opencl.blas.ujmp.CLDenseDoubleMatrix2D;
import com.nativelibs4java.opencl.blas.ujmp.CLDenseFloatMatrix2D;
import org.ujmp.core.Matrix;
import static org.ujmp.core.calculation.Calculation.ORIG;
import org.ujmp.core.mapper.MatrixMapper;

//http://nativelibs4java.sourceforge.net/javacl/api/stable/com/nativelibs4java/opencl/blas/ujmp/CLDenseDoubleMatrix2D.html
public class JavaCLMatrix {
    public static void main(String[] args) throws Exception  {
        CLContext context = JavaCL.createBestContext();
        System.out.println("Context: " + context);
        
        MatrixMapper.getInstance().setDenseFloatMatrix2DClassName(CLDenseFloatMatrix2D.class.getName());
        
        int dim = 100;
        
        CLDenseFloatMatrix2D a = new CLDenseFloatMatrix2D(dim, dim);
        a.zeros(ORIG);
        a.plus(ORIG, true, 1.0);
                
        CLDenseFloatMatrix2D b = new CLDenseFloatMatrix2D(dim, dim);
        b.zeros(ORIG);
        b.plus(ORIG, true, 2.0);

        
        Matrix c = a.times(b);
        
        System.out.println(c);
        
        
    }
}
