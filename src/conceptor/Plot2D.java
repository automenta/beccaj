package conceptor;

import info.monitorenter.gui.chart.Chart2D;
import info.monitorenter.gui.chart.ITrace2D;
import info.monitorenter.gui.chart.traces.Trace2DSimple;
import java.awt.Dimension;
import org.ejml.data.DenseMatrix64F;

/**
 *
 * @author me
 */
public class Plot2D extends Chart2D {

    private final DenseMatrix64F matrix;

    public Plot2D(DenseMatrix64F d, int xCol, int yCol, int xSize, int ySize) {
        super();
        this.matrix = d;
        // Create an ITrace: 
        ITrace2D trace = new Trace2DSimple();
        // Add the trace to the chart. This has to be done before adding points (deadlock prevention): 
        addTrace(trace);
        
        for (int i = 0; i < matrix.numRows; i++) {
            trace.addPoint(matrix.get(i, xCol), matrix.get(i, yCol));
        }

        setSize(xSize, ySize);
        setPreferredSize(new Dimension(xSize, ySize));
    }

}
