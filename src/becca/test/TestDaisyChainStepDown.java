package becca.test;

import becca.core.Util;
import org.ejml.data.DenseMatrix64F;


public class TestDaisyChainStepDown extends TestDaisyChain {

    double JITTER_FACTOR = 0.5;
    double NOISE_FACTOR = 0.1;
    
    public TestDaisyChainStepDown(int numCables) {
        super(numCables);                
        
    }

    public void update(double t) {
        double jitter = Math.random() * JITTER_FACTOR;
        t+=jitter;
        
        Util.setSinusoidal(cableActivities, 0, t, 8, 0);
        
        double phaseDelta = 0.1;
        for (int r = 0; r < cableGoalsIn.getNumCols(); r++) {
            Util.setSinusoidal(cableGoalsIn, r, t, 32.0, r * phaseDelta);
        }
        if (NOISE_FACTOR > 0)
            Util.addNoise(cableGoalsIn, NOISE_FACTOR);
        
        DenseMatrix64F chainActivities = d.stepUp(cableActivities);
        DenseMatrix64F cableGoalsOut = d.stepDown(cableGoalsIn);
        p.update(chainActivities, cableGoalsOut);
    }

    public static void main(String[] args) {
        new TestDaisyChainStepDown(64);
    }
    
}
