
package ql.som;

import static becca.core.Util.printArray;
import deeplearning.SdA.src.SdA;
import java.util.Random;
import ql.Perception;

/**
 *
 * @author me
 */
abstract public class SDAPerception extends Perception {

    private final SdA sda;

    double pretrain_lr = 0.25;
    double corruption_level = 0;
    int pretraining_epochs = 100;
    private final int n_outs;
    private final int n_ins;
    private final double[] encoded;
    private final double[] sensor;
    
    public SDAPerception(double[] sensor, int outputs) {
        super();
        
        this.sensor = sensor;
        
        Random rng = new Random(123);

        
        this.n_ins = sensor.length;
        this.n_outs = outputs;
        
        encoded = new double[n_outs];
        
        int[] hidden_layer_sizes = { (n_ins+outputs)*2,  outputs};
        int n_layers = hidden_layer_sizes.length;            

        int N = 1;
        
        // construct SdA
        sda = new SdA(N, n_ins, hidden_layer_sizes, n_outs, n_layers, rng);
        
    }

    @Override
    protected void updateInputValues() {
        updateSOM();
        
        for (int i = 0; i < encoded.length; i++)
            setNextValue(encoded[i]);        
    }

    
    public void updateSOM() {
        if (output == null) return;
        
        sda.pretrain(sensor, pretrain_lr, corruption_level, pretraining_epochs);
                                
        sda.encode(sensor, encoded);
        
        
        //printArray(output);
        //printArray(encoded);
    }

    
    
    
}
