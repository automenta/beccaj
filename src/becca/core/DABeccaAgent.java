/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package becca.core;

import static becca.core.Util.printArray;
import becca.gui.AgentPanel;
import becca.gui.MatrixPanel;
import becca.test.World;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import org.ejml.data.DenseMatrix64F;

/**
 * BECCA with Denoising Autoencoder Layer
 *
 * @author me
 */
public class DABeccaAgent extends BeccaAgent {

    public dA da;


    boolean display = false;
    int updatePeriodCycles = 256;

    double learning_rate = 0.1;
    double corruption_level = 0;

    private MatrixPanel daMatrixPanel;
    private JFrame window;
    private DenseMatrix64F dam;

    double[] reconstructedInput;

    int time = 0;
    private double[] encodedInput;

    public DABeccaAgent() {
        super();


    }

    public int getReducedSensors(int worldSensors) {
        //return worldSensors / 8;
        return (int)Math.ceil(Math.sqrt(worldSensors));
    }
    
    @Override
    public void init(World world) {

        
        int reducedSensors = numSensors = getReducedSensors(world.getNumSensors());
        
        System.out.println("World sensors: " + world.getNumSensors() + " -> Autoencoded Sensors: " + reducedSensors);
        
        this.da = new dA(world.getNumSensors(), reducedSensors);

        if (reconstructedInput == null) {
            reconstructedInput = new double[world.getNumSensors()];
            encodedInput = new double[reducedSensors];
        }

        if (display) {
            daMatrixPanel = new MatrixPanel();
            daMatrixPanel.setLayout(new BoxLayout(daMatrixPanel, BoxLayout.PAGE_AXIS));
            window = AgentPanel.window(daMatrixPanel, true);
            window.setSize(300, 100);
        }
        
        super.init(world);

    }

    public void updateDADisplay(double[] encodedInput) {
        dam = DenseMatrix64F.wrap(encodedInput.length, 1, encodedInput);
        daMatrixPanel.removeAll();
        daMatrixPanel.addMatrix("dA Encoded", dam);
        daMatrixPanel.validate();
        daMatrixPanel.repaint();
    }

    @Override
    public double[] getPercept() {
        return encodedInput;
    }

    
    @Override
    public int step(double reward) {

        
        da.train(sensor, learning_rate, corruption_level);


        da.getEncoded(sensor, encodedInput, false, true);        
        
            //printArray(sensor);
        //printArray(encodedInput);        
        da.reconstruct(sensor, reconstructedInput);
        //printArray(reconstructedInput);

        double diff = 0;
        for (int i = 0; i < reconstructedInput.length; i++) {
            diff += Math.abs(reconstructedInput[i] - sensor[i]);
        }

            //System.out.println("avg delta: " + (diff / ((double)reconstructedInput.length) ));
        if (display) {
            if (time % updatePeriodCycles == 0) {
                updateDADisplay(encodedInput);
            }
            time++;
        }

        return super.step(reward);
    }

}
