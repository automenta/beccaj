/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.gui;

import becca.core.Util;
import static becca.core.Util.m;
import java.awt.BorderLayout;
import java.awt.LayoutManager;
import java.util.HashMap;
import java.util.Map;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.EmptyBorder;
import org.ejml.data.DenseMatrix64F;
import static org.ejml.ops.CommonOps.transpose;
import org.ejml.ops.MatrixComponent;

/**
 *
 * @author me
 */
public class MatrixPanel extends JPanel {
    int maxMatrixSize = 128;
    protected JPanel target = this;
    protected Map<DenseMatrix64F, MatrixComponent> matrices = new HashMap();
    protected Map<DenseMatrix64F, JLabel> labels = new HashMap();

    public MatrixPanel(LayoutManager layout) {
        super(layout);
    }

    public MatrixPanel() {
        this((JPanel)null);
    }
    public MatrixPanel(JPanel target) {
        super();
        if (target == null)
            target = this;
        this.target = target;
    }
        
            
    public static JFrame window(JPanel content, boolean closeOnExit) {
        JFrame jf = new JFrame();        
        jf.setContentPane(new JScrollPane(content));
        jf.setSize(400,1000);
        jf.setVisible(true);
        if (closeOnExit)
            jf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);        
        return jf;
    }
        
    public void addMatrix(String id, DenseMatrix64F m) {
        if (m!=null)
            if ((m.getNumCols() > maxMatrixSize) || (m.getNumRows() > maxMatrixSize))
                return;

        JPanel x = new JPanel(new BorderLayout());
        x.setBorder(new EmptyBorder(8,8,8,8));
        x.setAlignmentX(JPanel.LEFT_ALIGNMENT);
        
        
        x.add(new JLabel(id + " " + (m != null ? m(m) : "") ), BorderLayout.NORTH);
        JLabel dlabel = new JLabel();

        if ((m!=null) && ((m.getNumCols() > 0) && (m.getNumRows() > 0))) {
            int px = 4; //min pixels per cell
            int maxPX = 8; //max pixels per cell
            if (m.getNumCols() < m.getNumRows())
                m = transpose(m, null);
            
            int w = (int)Math.max(m.getNumCols()*px, Math.log(m.getNumCols()*maxPX));
            int h = (int)Math.max(m.getNumRows()*px, Math.log(m.getNumRows()*maxPX));
            
            MatrixComponent mv = new MatrixComponent(w,h);            
            mv.setMatrix(m);
            matrices.put(m, mv);
            labels.put(m, dlabel);
            
            x.add(mv, BorderLayout.CENTER);
            x.add(dlabel, BorderLayout.SOUTH);
        }
        else {
            x.add(new JLabel("(null)"), BorderLayout.CENTER);
        }
        
        target.add(x);
        
        
    }
    
    public void refresh() {
        for (DenseMatrix64F m : matrices.keySet()) {
            MatrixComponent c = matrices.get(m);
            c.setMatrix(m);
            String t;
            
            JLabel l = labels.get(m);

            if (m!=null) {
                t = "+" + Util.elementSum(m);
            }
            else
                t = "";
            
            l.setText(t);
            
            c.repaint();
        }
    }
}
