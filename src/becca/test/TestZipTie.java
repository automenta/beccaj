/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.test;

import static becca.core.Util.setSinusoidal;
import becca.core.ZipTie;
import becca.gui.MatrixPanel;
import javax.swing.BoxLayout;
import javax.swing.JPanel;
import org.ejml.data.DenseMatrix64F;


/**
 *
 * @author me
 */
public class TestZipTie {
    protected final ZipTiePanel p;
    protected final ZipTie z;
    
    protected final DenseMatrix64F cableActivitiesIn;
    protected final DenseMatrix64F bundleGoalsIn;
    protected DenseMatrix64F bundleActivitiesOut;
    protected DenseMatrix64F cableGoalsOut;

    public static class ZipTiePanel extends MatrixPanel {
        private final ZipTie z;
    

        public ZipTiePanel(ZipTie z) {
            this(z, null);
        }
        public ZipTiePanel(ZipTie z, JPanel target) {
            super(target);
            
            
            this.z = z;
            
            setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
            
            update(null,null,null,null);
        }

        public void update(
                DenseMatrix64F cableActivitiesIn,
                DenseMatrix64F bundleGoalsIn,
                DenseMatrix64F bundleActivitiesOut,
                DenseMatrix64F cableGoalsOut
                ) {
            removeAll();
            
            if (cableActivitiesIn!=null)
                addMatrix("chainActivities (stepup input)", cableActivitiesIn);
            if (bundleGoalsIn!=null)
                addMatrix("bundleGoals (stepdown input)", bundleGoalsIn);
            
            if (bundleActivitiesOut!=null)
                addMatrix("bundleActivities (stepup output)", bundleActivitiesOut);
            
            if (cableGoalsOut!=null)
                addMatrix("cableGoals (stepdown output)", cableGoalsOut);

            
            addMatrix("bundleMap", z.getBundleMap());
            addMatrix("agglomerationEnergy", z.getAgglomerationEnergy());
            addMatrix("nucleationEnergy", z.getNucleationEnergy());
            addMatrix("bundleActivities", z.getNonBundleActivities());
            addMatrix("nonBundleActivities", z.getBundleActivities());         
            
            

            

            if (target==this) {
                if (getParent()!=null)
                    getParent().validate();
                validate();
                repaint();            
            }
        }
    }
    
    public TestZipTie(boolean inBlock, int maxCables, int maxBundles, int maxBundlesPerCable) {
        z = new ZipTie(inBlock, maxCables, maxBundles, maxBundlesPerCable);
        
        cableActivitiesIn = new DenseMatrix64F(maxCables, 1);
        bundleGoalsIn = new DenseMatrix64F(maxBundles, 1);
        
        p = new ZipTiePanel(z);
        
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
    
    
    //abstract public void update(double t);
    
    int cycle = 0;
    public void update(double t) {
        
        setSinusoidal(cableActivitiesIn, 0, t, 220.0, 0);
        setSinusoidal(bundleGoalsIn, 0, t, 10.0, 0.5, 1.0, 0.0);

        
        bundleActivitiesOut = z.stepUp(cableActivitiesIn);
        
        cableGoalsOut = z.stepDown(bundleGoalsIn);
        
        p.update(cableActivitiesIn, bundleGoalsIn, 
                bundleActivitiesOut, cableGoalsOut);
        
        /*if (cycle % 2500 == 0)
            p.update();
        cycle++;*/
    }
    
    public static void main(String[] args) {
        new TestZipTie(false, 32,16,64);
    }
}
