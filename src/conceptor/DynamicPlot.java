/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package conceptor;

import info.monitorenter.gui.chart.Chart2D;
import info.monitorenter.gui.chart.traces.Trace2DLtd;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 *
 * @author me
 */
public class DynamicPlot extends Chart2D {
    private Trace2DLtd[] trace;

    int history = 10000;
    
    public DynamicPlot(String... plotNames) {
        super();
        
        setMinimumSize(new Dimension(400, 200));
        setPreferredSize(new Dimension(400, 200));


        int plots = plotNames.length;
        
        trace = new Trace2DLtd[plots];
        for (int i = 0; i < plots; i++) {
            trace[i] = new Trace2DLtd(history);
            trace[i].setColor(Color.getHSBColor(((float)i)/((float)plots), 0.8f, 0.8f));
            trace[i].setName(plotNames[i]);
            addTrace(trace[i]);
        }
        


    }
    
    public void update(int x, double... values) {
        for (int i = 0; i < values.length; i++) {
            trace[i].addPoint(x, values[i]);
        }
        
    }
    

}
