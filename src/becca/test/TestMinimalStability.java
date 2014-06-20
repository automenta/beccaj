/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.test;

import becca.core.BeccaAgent;
import becca.world.Grid1DDiscrete;
import becca.world.Grid1DSimple;

/**
 *
 * @author me
 */
public class TestMinimalStability {
    private final Grid1DDiscrete g;
    private final Simulation s;

    public TestMinimalStability() throws Exception {
        Class<? extends Agent> a = BeccaAgent.class;
        
        g = new Grid1DDiscrete(2, 11990000);
        
        
        s = new Simulation(a, g, 200) {

            @Override
            public void init(Agent a) {
                super.init(a);
                this.displayPeriodMS = 50;
            }
          
        };
        
    
    }
    
    public static void main(String[] args) throws Exception {
        new TestMinimalStability();
    }
}
