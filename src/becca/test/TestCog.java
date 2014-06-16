package becca.test;

import becca.core.Cog;
import static becca.core.Util.setSinusoidal;
import becca.gui.MatrixPanel;
import becca.test.TestZipTie.ZipTiePanel;
import javax.swing.BoxLayout;
import org.ejml.data.DenseMatrix64F;


/**
 *
 * @author me
 */
public class TestCog {
    protected final CogPanel p;
    protected final Cog c;
    protected DenseMatrix64F cableActivitiesIn;
    protected DenseMatrix64F cableGoalsOut;
    protected DenseMatrix64F bundleGoalsIn;
    protected DenseMatrix64F bundleActivitiesOut;

    public class CogPanel extends MatrixPanel {        

        public CogPanel() {
            super();
            
            setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
                
            update();
        }

        public void update() {
            removeAll();
            
            if (cableActivitiesIn!=null)
                addMatrix("cableActivities (stepup input)", cableActivitiesIn);
            if (bundleGoalsIn!=null)
                addMatrix("cableGoals (stepdown input)", bundleGoalsIn);
            
            if (bundleActivitiesOut!=null)
                addMatrix("chainActivities (stepup output)", bundleActivitiesOut);
            if (cableGoalsOut!=null)
                addMatrix("cableGoals (stepdown output)", cableGoalsOut);

            
            addMatrix("surprise", c.getSurprise());


            addMatrix("daisychain.count", c.daisychain.getCount());
            
            
            ZipTiePanel zp = new TestZipTie.ZipTiePanel(c.ziptie, this);

                                    
            if (getParent()!=null)
                getParent().validate();
            validate();
            repaint();            
        }
    }
    
    public TestCog(int maxCables, int maxBundles, int level) {
        c = new Cog(maxCables, maxBundles, 0, level);
                
        p = new CogPanel();
        
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
        cableActivitiesIn = new DenseMatrix64F(c.maxCables, 1);
        bundleGoalsIn = new DenseMatrix64F(c.maxBundles, 1);        
    }
    
    
    public void update(double t) {
        setSinusoidal(cableActivitiesIn, 0, t, 10.0, 0);
        setSinusoidal(bundleGoalsIn, 0, t, 2, 0.5, 1.0, 0.5);
        
        bundleActivitiesOut = c.stepUp(cableActivitiesIn, true);
        
        cableGoalsOut = c.stepDown(bundleGoalsIn);
        
        
        p.update();
    }
    
    public static void main(String[] args) {
        
        new TestCog(16,8,1);
        
    }
    
    
}
