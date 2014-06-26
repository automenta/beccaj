/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rlpark;

import info.monitorenter.gui.chart.Chart2D;
import info.monitorenter.gui.chart.ITrace2D;
import info.monitorenter.gui.chart.traces.Trace2DLtd;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.Timer;
import java.util.TimerTask;
import javax.swing.JFrame;
import javax.swing.JPanel;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.MatrixComponent;

/**
 *
 * @author me
 */
abstract public class DynamicChart {

    public DynamicChart() {
        this(1500);
    }
    
    public DynamicChart(long delayMS) {
        super();
        
        // Create a chart:  
        Chart2D chart = new Chart2D();
        chart.setMinimumSize(new Dimension(400, 200));

    // Create an ITrace: 
        // Note that dynamic charts need limited amount of values!!! 
        final ITrace2D trace = new Trace2DLtd(1000);
        trace.setColor(Color.BLUE);

        // Add the trace to the chart. This has to be done before adding points (deadlock prevention): 
        chart.addTrace(trace);

        final JPanel p = new JPanel(new BorderLayout());
        p.add(chart, BorderLayout.CENTER);

        final JPanel s = new JPanel(new FlowLayout());
        s.setMinimumSize(new Dimension(400, 50));
        p.add(s, BorderLayout.SOUTH);

    // Make it visible:
        // Create a frame. 
        JFrame frame = new JFrame("MinimalDynamicChart");
        //frame.setTitle(name);
        // add the chart to the frame: 
        frame.getContentPane().add(p);
        frame.setSize(700, 800);
        // Enable the termination button [cross on the upper right edge]: 
        frame.addWindowListener(
                new WindowAdapter() {
                    public void windowClosing(WindowEvent e) {
                        System.exit(0);
                    }
                }
        );
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);

        /* 
         * Now the dynamic adding of points. This is just a demo! 
         * 
         * Use a separate thread to simulate dynamic adding of date. 
         * Note that you do not have to copy this code. Dynamic charting is just about 
         * adding points to traces at runtime from another thread. Whenever you hook on 
         * to a serial port or some other data source with a polling Thread (or an event 
         * notification pattern) you will have your own thread that just has to add points 
         * to a trace. 
         */
        Timer timer = new Timer(true);
        TimerTask task = new TimerTask() {

            private double m_y = 0;
            private long m_starttime = System.currentTimeMillis();
            private double m_x;

            /**
             * @see java.util.TimerTask#run()
             */
            @Override
            public void run() {

                this.m_y = getReward();
                this.m_x = getTime();
                // This is the important thing: Point is added from separate Thread.
                /*trace.addPoint(((double) System.currentTimeMillis() - this.m_starttime),
                        this.m_y);*/
                trace.addPoint(this.m_x,
                        this.m_y);

                s.removeAll();

                MatrixComponent ma = new MatrixComponent(300, 50);
                double[] S = getSensor();
                double[] A = getAction();
                double[] SA = new double[S.length + A.length];
                System.arraycopy(S, 0, SA, 0, S.length);
                System.arraycopy(A, 0, SA, S.length, A.length);
                ma.setMatrix(DenseMatrix64F.wrap(1, SA.length, SA));
                s.add(ma);

                p.doLayout();
                p.validate();
                p.repaint();
            }


        };
        // Every 20 milliseconds a new value is collected.
        timer.schedule(task, 1000, delayMS);
    }

    abstract public double[] getSensor();

    abstract public double[] getAction();

    abstract public double getReward();

    abstract public double getTime();
}
