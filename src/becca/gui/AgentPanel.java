/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.gui;

import becca.core.Agent;
import becca.core.Block;
import becca.core.Cog;
import becca.core.DaisyChain;
import javax.swing.JLabel;
import javax.swing.JPanel;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.MatrixComponent;

import static becca.core.Util.*;
import becca.core.ZipTie;
import java.awt.BorderLayout;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.border.EmptyBorder;

/**
 *
 * @author me
 */
public class AgentPanel extends MatrixPanel {
    private final Agent agent;

    boolean showZipTies = true;
    boolean showDaisyChains = true;
    
    
    public AgentPanel(Agent a) {
        super();

        
        BoxLayout bl = new BoxLayout(this, BoxLayout.PAGE_AXIS);        
        setLayout(bl);
        
        this.agent = a;
        
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
    
    public void update() {
        removeAll();
        
        addMatrix("sensors", DenseMatrix64F.wrap(agent.sensor.length, 1, agent.sensor));
        addMatrix("actions", DenseMatrix64F.wrap(agent.action.length, 1, agent.action));
            
        addMatrix("hub.count", agent.hub.getCount());                
        addMatrix("hub.cableActivities", agent.hub.getCableActivities());
        addMatrix("hub.chainActivities", agent.hub.getChainActivities());

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
                    String dPrefix = "daisychain.";
                    addMatrix(bPrefix + dPrefix + "surprise", d.getSurprise());
                    addMatrix(bPrefix + dPrefix + "reaction", d.getReaction());
                }
                
            }
                
        }
        
        if (getParent()!=null)
            getParent().validate();
        validate();
        repaint();
    }
    
}
