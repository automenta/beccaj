/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package javacl;

import com.nativelibs4java.opencl.*;
import com.nativelibs4java.opencl.CLMem.Usage;
import com.nativelibs4java.util.*;
import org.bridj.Pointer;
import static org.bridj.Pointer.*;
import static java.lang.Math.*;
import java.nio.ByteOrder;

public class JavaCLTutorial {
    public static void main(String[] args) throws Exception  {
        CLContext context = JavaCL.createBestContext();
        System.out.println("Context: " + context);
        
        CLQueue queue = context.createDefaultQueue();
        ByteOrder byteOrder = context.getByteOrder();
        
        int n = 1024;
        Pointer<Float>
            aPtr = allocateFloats(n).order(byteOrder),
            bPtr = allocateFloats(n).order(byteOrder);

        for (int i = 0; i < n; i++) {
            aPtr.set(i, (float)cos(i));
            bPtr.set(i, (float)sin(i));
        }

        // Create OpenCL input buffers (using the native memory pointers aPtr and bPtr) :
        CLBuffer<Float> 
            a = context.createBuffer(Usage.Input, aPtr),
            b = context.createBuffer(Usage.Input, bPtr);

        // Create an OpenCL output buffer :
        CLBuffer<Float> out = context.createBuffer(Usage.Output, Float.class, n);

        
        // Read the program sources and compile them :
        String src = IOUtils.readText(JavaCLTutorial.class.getResource("TutorialKernels.cl"));
        CLProgram program = context.createProgram(src);

        // Get and call the kernel :
        CLKernel addFloatsKernel = program.createKernel("add_floats");
        addFloatsKernel.setArgs(a, b, out, n);
        CLEvent addEvt = addFloatsKernel.enqueueNDRange(queue, new int[] { n });
        
        Pointer<Float> outPtr = out.read(queue, addEvt); // blocks until add_floats finished

        // Print the first 10 output values :
        for (int i = 0; i < 10 && i < n; i++)
            System.out.println("out[" + i + "] = " + outPtr.get(i));
        
    }
}
