/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.test;



/**
 *
 * @author me
 */
public interface Agent {
    
    public int step(double reward);
    public double[] getSensor();
    public double[] getAction();

    public void init(World world);
    
}
