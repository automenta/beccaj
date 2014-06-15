package becca.test;

import becca.core.Block;
import becca.core.Hub;
import static becca.core.Util.setSinusoidal;
import becca.gui.MatrixPanel;
import java.util.ArrayList;
import java.util.LinkedList;
import javax.swing.BoxLayout;
import org.ejml.data.DenseMatrix64F;

/**
 *
 * @author me
 */
public class TestHub {
    
    protected final HubPanel p;
    protected final Hub h;
    protected final DenseMatrix64F cableActivities; //in
    protected final DenseMatrix64F hubCableGoals; //out
    ArrayList<Block> blocks = new ArrayList();
    private final int numBlocks;

    public class HubPanel extends MatrixPanel {
        private final Hub h;
        private final ArrayList<Block> blocks;

        public HubPanel(Hub h, ArrayList<Block> blocks) {
            super();
            
            setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
            
            this.h = h;
            this.blocks = blocks;
    
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
            
            
            for (int i = 0; i < blocks.size(); i++) {
                addMatrix("block" + i + ".cableActivities (in)", 
                        blocks.get(i).getCableActivities());
                addMatrix("block" + i + ".hubCableGoals (out)", 
                        blocks.get(i).getHubCableGoals());
            }

            addMatrix("cableActivities", h.getCableActivities());
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
    
    public TestHub(int numBlocks, int numCables) {
        h = new Hub(numCables);
        
        this.numBlocks = numBlocks;

        cableActivities = new DenseMatrix64F(numCables, 1);
        hubCableGoals = new DenseMatrix64F(numCables, numCables);
        
        for (int i = 0; i < numBlocks; i++) {
            int cablesPerBlock = numCables/numBlocks;
            blocks.add(new Block(cablesPerBlock, cablesPerBlock, i));            
        }
         
        p = new HubPanel(h, blocks);
        
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
        
        int f = blocks.size();
        for (Block b : blocks) {
            setSinusoidal(b.getCableActivities(), 0, t, 5.0*f, 0);
            f--;
        }
        
        h.step(blocks, r);

        p.update();
        
        for (Block b : blocks) {
            b.goalDecay();
        }
        
        
    }    
    
    public static void main(String[] args) {
        new TestHub(1, 32);
    }
    
}
