/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.test;

import becca.core.Agent;

/**
 *
 * @author me
 */
public class TestAgent {

    public TestAgent() {
        Agent a = new Agent("a", 8, 8);
        System.out.println(a);
    }
    
    public static void main(String[] args) {
        new TestAgent();
    }
    
}
