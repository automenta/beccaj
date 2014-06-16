package becca.test;

import becca.core.Hub;
import static becca.core.Util.setSinusoidal;
import becca.gui.MatrixPanel;
import java.util.LinkedList;
import javax.swing.BoxLayout;
import org.ejml.data.DenseMatrix64F;
import static org.ejml.ops.CommonOps.scale;

/**
 *
 * @author me
 */
public class TestHubDirect {
    
    protected final HubPanel p;
    protected final Hub h;
    protected final DenseMatrix64F cableActivities; //in
    protected final DenseMatrix64F hubCableGoals; //out

    public class HubPanel extends MatrixPanel {
        private final Hub h;

        public HubPanel(Hub h) {
            super();
            
            setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
            
            this.h = h;
    
            update();
            
        }

        public void update() {
            removeAll();
            
            /*if (cableActivities!=null)
                addMatrix("cableActivities (in)", cableActivities);
            if (hubCableGoals!=null)
                addMatrix("hubCableGoals (out)", h.hubCableGoals);*/
            
            
            double[] rt = h.getRewardArray();
            addMatrix("rewardTrace", DenseMatrix64F.wrap(rt.length, 1, rt));
            
            addMatrix("cableActivities", h.getCableActivities());
            addMatrix("output", h.getOutput());
            
            addMatrix("chainActivities", h.getChainActivities());
            addMatrix("expectedReward", h.getExpectedReward());
            addMatrix("estimatedRewardValue", h.getEstimatedRewardValue());

            addMatrix("count", h.getCount());
            
            LinkedList<DenseMatrix64F> pres = h.getPre();
            LinkedList<DenseMatrix64F> posts = h.getPre();
            for (int i = 0; i < pres.size(); i++) {
                addMatrix("pre." + i, pres.get(i));
            }
            for (int i = 0; i < posts.size(); i++) {
                addMatrix("post." + i, posts.get(i));            
            }
            
            if (getParent()!=null)
                getParent().validate();
            validate();
            repaint();            
        }
    }
    
    public TestHubDirect(int numCables) {
        h = new Hub(numCables);

        cableActivities = new DenseMatrix64F(numCables, 1);
        hubCableGoals = new DenseMatrix64F(numCables, numCables);
                
        p = new HubPanel(h);
        
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
    }
    
   
    protected void update(double t) {
                
        double r = Math.sin(t);
        
        setSinusoidal(h.getCableActivities(), 0, t, 5.0, 1.0);        
        h.step(null, r);
        
        scale(0.1, h.getOutput());

        p.update();
        
    }    
    
    public static void main(String[] args) {        
        new TestHubDirect(8); 
    }
    
}
