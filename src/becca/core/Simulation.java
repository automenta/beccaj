/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.core;

import becca.gui.AgentPanel;
import java.util.Date;
import javax.swing.JFrame;
import javax.swing.JScrollPane;



/**
 *
 * @author me
 */
public class Simulation {

    public final Agent agent;
    private AgentPanel ap;
    private JFrame jf;
    
    /*
    Run BECCA with world.  

    If restore=True, this method loads a saved agent if it can find one.
    Otherwise it creates a new one. It connects the agent and
    the world together and runs them for as long as the 
    world dictates.

    To profile BECCA's performance with world, manually set
    profile_flag in the top level script environment to True.

    */
    public void displayAgent(Agent a) {
        jf = new JFrame();
        ap = new AgentPanel(a);
        jf.setContentPane(new JScrollPane(ap, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS));
        jf.setSize(700,1000);
        jf.setVisible(true);
        jf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
    }
    private double reward;
    
    long displayPeriodMS = 50;
    long lastDisplay = -1;
    
    public Simulation(World world) {
        
        /*
        if agent_name is None:
            agent_name = '_'.join((world.name, 'agent'))
        */
        this.agent = new Agent(world.getName() + " Agent", 
                            world.getNumSensors(), world.getNumActions());
        
        displayAgent(agent);
        
        /*if restore:
            agent = agent.restore()*/
        
        /*
        # If configured to do so, the world sets some BECCA parameters to 
        # modify its behavior. This is a development hack, and 
        # should eventually be removed as BECCA matures and settles on 
        # a good, general purpose set of parameters.
        world.set_agent_parameters(agent)        

        action = np.zeros((world.num_actions,1))
        */
        
        while (world.isActive()) {
            /*    
            # Repeat the loop through the duration of the existence of the world 
                sensor, reward = world.step(action)
                world.visualize(agent)
                action = agent.step(sensor, reward)
            return agent.report_performance()
            */

            reward = world.step(agent.action, agent.sensor);

            agent.step(reward);

            long now = System.currentTimeMillis();
            if ((lastDisplay == -1) || (now - lastDisplay > displayPeriodMS)) {
                ap.update();
                jf.setTitle("Reward: " + reward);
                lastDisplay = System.currentTimeMillis();
            }
        }
        
        System.out.println("Simulation Finished.");
                          
    }
}
