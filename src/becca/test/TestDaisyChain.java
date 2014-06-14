/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.test;

import becca.core.DaisyChain;
import becca.gui.MatrixPanel;
import javax.swing.BoxLayout;
import org.ejml.data.DenseMatrix64F;


/**
 *
 * @author me
 */
public class TestDaisyChain {
    private final DenseMatrix64F cableActivities;

    public class DaisyChainPanel extends MatrixPanel {
        private final DaisyChain d;
        private final DenseMatrix64F cableActivities;

        public DaisyChainPanel(DaisyChain d, DenseMatrix64F cableActivities) {
            super();
            
            setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
            
            this.d = d;
            this.cableActivities = cableActivities;
    
            update((DenseMatrix64F)null);
        }

        public void update(DenseMatrix64F chainActivities) {
            removeAll();
            
            addMatrix("cableActivities (stepup input)", cableActivities);
            if (chainActivities!=null)
                addMatrix("chainActivities (stepup output)", chainActivities);

            
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
        final DaisyChain d = new DaisyChain(numCables);
        
        cableActivities = new DenseMatrix64F(numCables, 1);
        
        final DaisyChainPanel p = new DaisyChainPanel(d, cableActivities);
        MatrixPanel.window(p, true);
        
        init();
        
        new Thread(new Runnable() {

            @Override
            public void run() {
                long startTime = System.currentTimeMillis();
                while (true) {
        
                    DenseMatrix64F chainActivities = d.stepUp(cableActivities);
                    
                    p.update(chainActivities);
                    
                    
                    try {
                        Thread.sleep(50);
                    } catch (InterruptedException ex) { }
                            
                    update((System.currentTimeMillis() - startTime)/1000.0); 
                    
                }
            }
            
        }).start();
    }
   
    public void init() {
        cableActivities.set(0,0,1.0);
    }
    
    public void update(double t) {
        //normRand(cableActivities, 1.0, 0.0);
        cableActivities.set(0, 0, Math.sin(t)/2.0+1.0);
        cableActivities.set(1, 0, Math.sin(t*2.0)/2.0+1.0);
        cableActivities.set(2, 0, Math.sin(t*4.0)/2.0+1.0);
        cableActivities.set(3, 0, Math.sin(t*8.0)/2.0+1.0);        
    }
    
    public static void main(String[] args) {
        new TestDaisyChain(4);
    }
}
