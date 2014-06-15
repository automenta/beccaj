/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.test;

import becca.core.DaisyChain;
import static becca.core.Util.normRand;
import becca.gui.MatrixPanel;
import javax.swing.BoxLayout;
import org.ejml.data.DenseMatrix64F;


/**
 *
 * @author me
 */
abstract public class TestDaisyChain {
    protected final DaisyChainPanel p;
    protected final DaisyChain d;
    protected final DenseMatrix64F cableActivities;
    protected final DenseMatrix64F cableGoalsIn;

    public class DaisyChainPanel extends MatrixPanel {
        private final DaisyChain d;

        public DaisyChainPanel(DaisyChain d) {
            super();
            
            setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
            
            this.d = d;
    
            update((DenseMatrix64F)null, (DenseMatrix64F)null);
        }

        public void update(DenseMatrix64F chainActivities, DenseMatrix64F cableGoalsOut) {
            removeAll();
            
            if (cableActivities!=null)
                addMatrix("cableActivities (stepup input)", cableActivities);
            if (cableGoalsIn!=null)
                addMatrix("cableGoals (stepdown input)", cableGoalsIn);
            
            if (chainActivities!=null)
                addMatrix("chainActivities (stepup output)", chainActivities);
            if (cableGoalsOut!=null)
                addMatrix("cableGoals (stepdown output)", cableGoalsOut);

            
            addMatrix("count", d.getCount());

            addMatrix("expectedCableActivations", d.getExpectedCableActivities());
            addMatrix("surprise", d.getSurprise());
            addMatrix("reaction", d.getReaction());            
            
            addMatrix("pre", d.getPre());
            addMatrix("preCount", d.getPreCount());
            addMatrix("post", d.getPost());
            addMatrix("postUncertainty", d.getPostUncertainty());
            
            if (getParent()!=null)
                getParent().validate();
            validate();
            repaint();            
        }
    }
    
    public TestDaisyChain(int numCables) {
        d = new DaisyChain(numCables);
        
        cableActivities = new DenseMatrix64F(numCables, 1);
        cableGoalsIn = new DenseMatrix64F(numCables, numCables);
        
        p = new DaisyChainPanel(d);
        
        MatrixPanel.window(p, true);
        
        init();
        
        new Thread(new Runnable() {

            @Override
            public void run() {
                long startTime = System.currentTimeMillis();
                while (true) {

                    update((System.currentTimeMillis() - startTime)/1000.0); 

                    try {
                        Thread.sleep(50);
                    } catch (InterruptedException ex) { }
                    
                }
            }
            
        }).start();
    }
   
    public void init() {
        cableActivities.set(0,0,1.0);
    }
    
    protected void setRandom(double t) {
        normRand(cableActivities, 1.0, 0.0);        
    }
    
    protected void setSinusoidal(int col, double t, DenseMatrix64F m, double baseFreq, double phase) {
        
        for (int i = 0; i < m.getNumRows(); i++) {
            m.set(i, col, Math.sin(phase + t/((double)(1+i))/3.14159*baseFreq)/2.0+1.0);                
        }
    }
    
    abstract public void update(double t);
    
    
}
