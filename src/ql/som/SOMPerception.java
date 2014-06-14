
package ql.som;

import org.encog.ml.data.basic.BasicMLData;
import org.encog.ml.data.basic.BasicMLDataSet;
import org.encog.neural.som.SOM;
import org.encog.neural.som.training.basic.BasicTrainSOM;
import org.encog.neural.som.training.basic.neighborhood.NeighborhoodSingle;
import ql.Perception;

/**
 *
 * @author me
 */
abstract public class SOMPerception extends Perception {

    double[] somoutput;
    
    public SOMPerception(int inputs, int outputs) {
        super();

        somoutput = new double[outputs];        
        
            
        
    }


    public void updateSOM() {
        if (output == null) return;
        
        
        
    }

    
    
    
}
