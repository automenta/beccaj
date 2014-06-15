package becca.test;

import org.ejml.data.DenseMatrix64F;
import static org.ejml.ops.CommonOps.scale;


public class TestDaisyChainStepDown extends TestDaisyChain {

    public TestDaisyChainStepDown(int numCables) {
        super(numCables);                
        
        scale(0, cableActivities);
    }

    public void update(double t) {
        setSinusoidal(0, t, cableActivities, 8, 0);
        
        double phaseDelta = 0.1;
        for (int r = 0; r < cableGoalsIn.getNumCols(); r++) {
            setSinusoidal(r, t, cableGoalsIn, 32.0, r * phaseDelta);
        }
        
        System.out.println(cableActivities);
        DenseMatrix64F chainActivities = d.stepUp(cableActivities);
        DenseMatrix64F cableGoalsOut = d.stepDown(cableGoalsIn);
        p.update(chainActivities, cableGoalsOut);
    }

    public static void main(String[] args) {
        new TestDaisyChainStepDown(64);
    }
    
}
