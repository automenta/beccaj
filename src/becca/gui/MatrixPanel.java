/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.gui;

import static becca.core.Util.m;
import java.awt.BorderLayout;
import java.awt.LayoutManager;
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
    int maxMatrixSize = 64;
    protected JPanel target = this;

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
        jf.setContentPane(new JScrollPane(content, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS));
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

        if ((m!=null) && ((m.getNumCols() > 0) && (m.getNumRows() > 0))) {
            int px = 6; //min pixels per cell
            int maxPX = 8; //max pixels per cell
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
        
        target.add(x);
        
    }
    
}
