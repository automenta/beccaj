/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package ql;

import becca.test.Agent;
import becca.test.World;
import java.util.Arrays;

/**
 *
 * @author me
 */
public class QLAgent implements Agent {
    private Action[] qaction;
    private Brain brain;
    double nextReward;
    double[] sensor;
    double[] action;

    @Override
    public void init(World world) {
        int sensors = world.getNumSensors();
        int actions = world.getNumActions();
        
        qaction = new Action[actions];
        for (int i = 0; i < actions; i++) {
            final int I = i;
            qaction[i] = new Action() {                
                @Override public int execute() { return I; }
            };
        }
        
        sensor = new double[sensors];
        action = new double[actions];
        
        brain = new Brain(new Perception() {

            @Override
            public boolean isUnipolar() {
                return true;
            }

            @Override
            public double getReward() {
                return nextReward;
            }

            @Override
            protected void updateInputValues() {
                for (int i = 0; i < sensor.length; i++)
                    setNextValue(sensor[i]);        
            }
        }, qaction, new int[] { 20,16 } );
        
        /*brain = new Brain(new SDAPerception(sensor, 4) {

            @Override
            public boolean isUnipolar() {
                return true;
            }

            @Override
            public double getReward() {
                return nextReward;
            }
            
        }, qaction ); */
        
                
        brain.reset();
    }

    
    @Override
    public int step(double reward) {
        this.nextReward = reward;

        brain.getPerception().perceive();
        brain.count();
        
        Arrays.fill(action, 0.0);
        int a = brain.getAction();
        action[a] = 1.0;
        
        //brain.printStats();
        //System.out.println(reward + " " + a);
        //Util.printArray(brain.getInput());
        //Util.printArray(brain.getOutput());
        //Util.printArray(action);
        
        return 0;
    }

    @Override
    public double[] getSensor() {
        return sensor;
    }

    @Override
    public double[] getAction() {
        return action;
    }
    
}
