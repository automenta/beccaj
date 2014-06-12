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
import javax.swing.border.EmptyBorder;

/**
 *
 * @author me
 */
public class AgentPanel extends JPanel {
    private final Agent agent;

    public AgentPanel(Agent a) {
        super();
        
        BoxLayout bl = new BoxLayout(this, BoxLayout.PAGE_AXIS);        
        setLayout(bl);
        
        this.agent = a;
        
        update();
    }
    
    protected void addMatrix(String id, DenseMatrix64F m) {
        JPanel x = new JPanel(new BorderLayout());
        x.setBorder(new EmptyBorder(8,8,8,8));
        x.setAlignmentX(JPanel.LEFT_ALIGNMENT);
        
        x.add(new JLabel(id + " " + (m != null ? m(m) : "") ), BorderLayout.NORTH);

        if ((m!=null) && ((m.getNumCols() > 0) && (m.getNumRows() > 0))) {
            int px = 8; //min pixels per cell
            int maxPX = 12; //max pixels per cell
            if (m.getNumCols() == 1)
                m = transpose(m, null);
            
            int w = (int)Math.max(m.getNumCols()*px, Math.log(m.getNumCols()*maxPX));
            int h = (int)Math.max(m.getNumRows()*px, Math.log(m.getNumRows()*maxPX));
            
            MatrixComponent mv = new MatrixComponent(w,h);
            mv.setMatrix(m);
            
            x.add(mv, BorderLayout.CENTER);
        }
        else {
            x.add(new JLabel("(null)"), BorderLayout.CENTER);
        }
        
        add(x);
        
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
        
        int n = 0;
        for (Block b : agent.blocks) {
            String bPrefix = "b" + (b.level) + ".";
            addMatrix(bPrefix + "cableActivities", b.getCableActivities());
            addMatrix(bPrefix + "surprise", b.getSurprise());
            
            addZipTie(bPrefix, b.ziptie);
            
                                    
            for (int i = 0; i < b.cogs.size(); i++) {
                Cog c = b.cogs.get(i);
                String cPrefix = "c" + i + ".";
                addMatrix(bPrefix + cPrefix + "surprise", c.getSurprise());
                
                ZipTie z = c.ziptie;
                if (z!=null)
                    addZipTie(cPrefix, z);
                
                DaisyChain d = c.daisychain;
                
            }
                
        }
        
        if (getParent()!=null)
            getParent().validate();
        validate();
        repaint();
    }
    
}
