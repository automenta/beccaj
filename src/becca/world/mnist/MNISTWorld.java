/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.world.mnist;

import becca.core.BeccaAgent;
import becca.core.DABeccaAgent;
import becca.test.Simulation;
import becca.test.World;
import java.io.IOException;

/**
 *
 * @author me
 */
public class MNISTWorld extends MNIST implements World {

    int currentImage = 0, currentFrame = -1;
    int cycle = 0;
    
    
    public MNISTWorld(String path, int maxImages) throws IOException {
        super("/home/me/Downloads", maxImages);
    }

    
    @Override
    public String getName() {
        return "MNIST";
    }

    @Override
    public int getNumSensors() {
        return 28*28;
    }

    @Override
    public int getNumActions() {
        return 4; //4 bits
    }

    @Override
    public boolean isActive() {
        return true;
    }

    int bits[] = new int[4];    
    
    MNISTImage i;
    MNISTImage blank = new MNISTImage(28,28);
    
    MNISTImage retina = new MNISTImage(28,28);
    int nextColumn = 0;
    
    @Override
    public double step(double[] action, double[] sensor) {
        
        if (cycle % trainingCyclesPerImage == 0) {
            currentFrame++;
            if (currentFrame % 10 == 0) {
                i = blank;
                if (trainingCyclesPerImage < maxTrainingCyclesPerImage)
                    trainingCyclesPerImage+=5;
            }
            else {
                i = images.get(currentImage++);                
                currentImage %= images.size();
            }
            nextColumn = 0;
        }
        
        if (nextColumn<28) {
            if (cycle % scrollCycles == 0) {
                retina.scrollRight(i, nextColumn);
                nextColumn++;
            }
        }
        else {
        }
        
        
        retina.toArray(sensor, noise);
        
        
        
        double threshold = 0.75;
        
        int a = 0;
        int factor = 1;
        for (int x = 0; x < action.length; x++) {
            boolean active = (action[x] > threshold);
            if (active) {
                a += factor;
            }
            factor *= 2;
        }
        a = a - 1;
        
        double r;
        
        if (i.label == -1) {
            if (a == -1) r = 1.0;
            else r = -1.0;
        }
        else {
            if ((a < 0) || (a > 9)) r = -1.0; 
            else {


                //if (a == i.label) r = 1.0;
                //else r = 0;
                

                r = 1.0 - (Math.abs(a - i.label)/10.0)*2.0;
                r *= r;
                
            }
        }
      
        /*
        System.out.print(cycle + " " + currentFrame + " " + currentImage + " label=" + i.label + ": " + a + " " + r + " [");
        printArray(action);
        */
        
        cycle++;
        
        return r;
    }
    
    int scrollCycles = 2;
    int maxTrainingCyclesPerImage = 1000, trainingCyclesPerImage = 0+scrollCycles*28*2;
    final static double noise = 0.05;
    
    public static void main(String[] args) throws IOException, Exception {
        
        MNISTWorld m = new MNISTWorld("/home/me/Downloads", 200);
        
        
        
        DABeccaAgent a = new DABeccaAgent() {            
            @Override
            public int getReducedSensors(int worldSensors) {
                //return (int)Math.sqrt(worldSensors);
                return 32;
            }            

            @Override
            public void init(World world) {
                super.init(world); //To change body of generated methods, choose Tools | Templates.
                pretrain(m.getImageVectors(), 5, 2, 0.1, noise);
            }


            @Override
            public void update(double lastReward, int time) {
                super.update(lastReward, time);
                
                double e = 0.15 + 1/(1.0 + time/1000.0)*0.25;
                if (time%1000 == 0)
                    System.out.println(time + " exploration=" + e);
                getHub().setEXPLORATION(e);
                
            }
        };
        
        
        BeccaAgent b = new BeccaAgent() {            


            @Override
            public void update(double lastReward, int time) {
                super.update(lastReward, time);
                
                double e = 0.15 + 1/(1.0 + time/1000.0)*0.7;
                if (time%1000 == 0)
                    System.out.println(time + " exploration=" + e);
                getHub().setEXPLORATION(e);
                
            }
        };
        
                
        
        
        new Simulation(b, m, 0);
        
    }

}
