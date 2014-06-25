/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.gui;

import becca.core.BeccaAgent;
import becca.core.Block;
import becca.core.Cog;
import becca.core.DaisyChain;
import becca.core.ZipTie;
import javax.swing.BoxLayout;
import org.ejml.data.DenseMatrix64F;



/**
 *
 * @author me
 */
public class AgentPanel extends MatrixPanel {
    private final BeccaAgent agent;

    boolean showZipTies = true;
    boolean showDaisyChains = true;
    private final DenseMatrix64F sensorMatrix;
    private final DenseMatrix64F actionMatrix;
    
    int recreatePeriod = 100;
    
    public AgentPanel(BeccaAgent a) {
        super();

        this.agent = a;
        
        sensorMatrix = new DenseMatrix64F(agent.getPercept().length, 1);

        actionMatrix = new DenseMatrix64F(agent.action.length, 1);

        reinit();
    }

    public void reinit() {
        removeAll();
        matrices.clear();
        labels.clear();
        
        BoxLayout bl = new BoxLayout(this, BoxLayout.PAGE_AXIS);        
        setLayout(bl);
        
        
        addMatrix("sensors", sensorMatrix);
        addMatrix("actions", actionMatrix);
            
        addMatrix("hub.output", agent.hub.getOutput());       
        addMatrix("hub.count", agent.hub.getCount());                
        addMatrix("hub.cableActivities", agent.hub.getCableActivities());
        addMatrix("hub.chainActivities", agent.hub.getChainActivities());
        addMatrix("hub.expectedReward", agent.hub.getExpectedReward());
        addMatrix("hub.estimatedRewardValue", agent.hub.getEstimatedRewardValue());
        
                
        
        int n = 0;
        for (Block b : agent.blocks) {
            String bPrefix = "b" + (b.level) + ".";
            addMatrix(bPrefix + "hubCableGoals", b.getHubCableGoals());
            addMatrix(bPrefix + "bundleActivities", b.getBundleActivities());
            addMatrix(bPrefix + "cableActivities", b.getCableActivities());
            addMatrix(bPrefix + "surprise", b.getSurprise());
            
            addZipTie(bPrefix, b.ziptie);
            
                                    
            for (int i = 0; i < b.cogs.size(); i++) {
                Cog c = b.cogs.get(i);
                String cPrefix = "c" + i + ".";
                addMatrix(bPrefix + cPrefix + "surprise", c.getSurprise());
                
                ZipTie z = c.ziptie;
                if ((z!=null) && (showZipTies))
                    addZipTie(cPrefix, z);
    
                if (showDaisyChains) {
                    DaisyChain d = c.daisychain;
                    String dPrefix = bPrefix + cPrefix + "daisychain.";
                    addMatrix(dPrefix + "surprise", d.getSurprise());
                    addMatrix(dPrefix + "reaction", d.getReaction());
                    addMatrix(dPrefix + "expectedCableActivations", d.getExpectedCableActivities());
          
                    addMatrix(dPrefix + "pre", d.getPre());
                    addMatrix(dPrefix + "preCount", d.getPreCount());
                    addMatrix(dPrefix + "post", d.getPost());
                    addMatrix(dPrefix + "postUncertainty", d.getPostUncertainty());
                }
                
            }
                
        }
        
        if (getParent()!=null)
            getParent().validate();
        validate();
        repaint();
        
        update();
        
    }
    
    protected void addZipTie(String bPrefix, ZipTie z) {
        addMatrix(bPrefix + "ziptie.bundleActivities" , 
                z.getBundleActivities());
        addMatrix(bPrefix + "ziptie.nonBundleActivities" , 
                z.getNonBundleActivities());            
        addMatrix(bPrefix + "ziptie.cableActivities" , 
                z.getCableActivities());
        addMatrix(bPrefix + "ziptie.bundleMap" , 
                z.getBundleMap());
        addMatrix(bPrefix + "ziptie.agglomerationEnergy" , 
                z.getAgglomerationEnergy());
        addMatrix(bPrefix + "ziptie.nucleationEnergy" , 
                z.getNucleationEnergy());        
    }
    
    int count = 1;
    public void update() {
        System.arraycopy(agent.getPercept(), 0, sensorMatrix.data, 0, agent.getPercept().length);
        System.arraycopy(agent.getAction(), 0, actionMatrix.data, 0, agent.getAction().length);
        //DenseMatrix64F.wrap(agent.getPercept().length, 1, agent.getPercept()        //DenseMatrix64F.wrap(agent.action.length, 1, agent.action)

        if (count++ % recreatePeriod == 0)
            reinit();
        
        refresh();
        repaint();
        
        
    }

    
}
