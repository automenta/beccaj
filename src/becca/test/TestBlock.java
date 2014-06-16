/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.test;

import becca.core.Block;
import static becca.core.Util.setSinusoidal;
import becca.gui.MatrixPanel;
import javax.swing.BoxLayout;
import javax.swing.JPanel;
import org.ejml.data.DenseMatrix64F;


/**
 *
 * @author me
 */
public class TestBlock {
    protected final BlockPanel p;
    protected final Block b;
    
    
    protected DenseMatrix64F cableActivitiesIn;
    protected DenseMatrix64F bundleGoalsIn;
    protected DenseMatrix64F bundleActivitiesOut;
    protected DenseMatrix64F cableGoalsOut;

    public class BlockPanel extends MatrixPanel {
        private final Block b;
    

        public BlockPanel(Block b) {
            this(b, null);
        }
        public BlockPanel(Block b, JPanel target) {
            super(target);
                        
            this.b = b;
            
            setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
            
            update();
        }

        public void update() {
            removeAll();
            
            
            if (cableActivitiesIn!=null)
                addMatrix("cableActivities (stepup input)", cableActivitiesIn);
            if (bundleGoalsIn!=null)
                addMatrix("bundleGoals (stepdown input)", bundleGoalsIn);         
            if (bundleActivitiesOut!=null)
                addMatrix("bundleActivities (stepup output)",bundleActivitiesOut);             if (cableGoalsOut!=null)
                addMatrix("cableGoals (stepdown output)", cableGoalsOut);
            
            
            //addMatrix("bundleActivities", b.getBundleActivities()); //Same as above
            addMatrix("cableActivities", b.getCableActivities());
            //addMatrix("hubCableGoals", b.getHubCableGoals()); //same as above
            addMatrix("surprise", b.getSurprise());
            
            

            

            if (target==this) {
                if (getParent()!=null)
                    getParent().validate();
                validate();
                repaint();            
            }
        }
    }
    
    public TestBlock(int minCables, int maxCables, int level) {
        b = new Block(minCables, maxCables, level);
                
        p = new BlockPanel(b);

        cableActivitiesIn = new DenseMatrix64F(maxCables, 1);
        bundleGoalsIn = new DenseMatrix64F(maxCables, 1);

        MatrixPanel.window(p, true);
        
        init();
        
        new Thread(new Runnable() {

            @Override
            public void run() {
                long startTime = System.currentTimeMillis();
                
                while (true) {

                    update((System.currentTimeMillis() - startTime)/1000.0); 

                    try {
                        Thread.sleep(10);
                    } catch (InterruptedException ex) { } 
                    
                }
            }
            
        }).start();
    }
   
    public void init() {
    }
    
    
    public void update(double t) {
        
        setSinusoidal(cableActivitiesIn, 0, t, 220.0, 0);
        setSinusoidal(bundleGoalsIn, 0, t, 10.0, 0.5, 1.0, 0.0);
        
        
        bundleActivitiesOut = b.stepUp(cableActivitiesIn);
        cableGoalsOut = b.stepDown(bundleGoalsIn);
        
        p.update();
        
    }
    
    public static void main(String[] args) {
        new TestBlock(16,32,0);
    }
}
