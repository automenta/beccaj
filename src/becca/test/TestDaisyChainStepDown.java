package becca.test;

import becca.core.Util;
import org.ejml.data.DenseMatrix64F;
import static org.ejml.ops.CommonOps.scale;


public class TestDaisyChainStepDown extends TestDaisyChain {

    public TestDaisyChainStepDown(int numCables) {
        super(numCables);                
        
    }

    public void update(double t) {
        Util.setSinusoidal(cableActivities, 0, t, 8, 0);
        
        double phaseDelta = 0.1;
        for (int r = 0; r < cableGoalsIn.getNumCols(); r++) {
            Util.setSinusoidal(cableGoalsIn, r, t, 32.0, r * phaseDelta);
        }
        
        DenseMatrix64F chainActivities = d.stepUp(cableActivities);
        DenseMatrix64F cableGoalsOut = d.stepDown(cableGoalsIn);
        p.update(chainActivities, cableGoalsOut);
    }

    public static void main(String[] args) {
        new TestDaisyChainStepDown(64);
    }
    
}
