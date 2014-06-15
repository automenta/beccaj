/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.test;

import org.ejml.data.DenseMatrix64F;

/**
 *
 * @author me
 */
public class TestDaisyChainStepUp extends TestDaisyChain {

    public TestDaisyChainStepUp(int numCables) {
        super(numCables);
    }

    public void update(double t) {
        setSinusoidal(0, t, cableActivities, 32.0, 0);
        DenseMatrix64F chainActivities = d.stepUp(cableActivities);
        p.update(chainActivities, null);
    }

    public static void main(String[] args) {
        new TestDaisyChainStepUp(64);
    }
    
}
