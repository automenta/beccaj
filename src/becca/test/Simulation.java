/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.test;

import becca.core.BeccaAgent;
import javax.swing.JFrame;
import rlpark.DynamicChart;
import becca.gui.AgentPanel;

/**
 *
 * @author me
 */
public class Simulation {

    public final Agent agent;
    private AgentPanel ap;
    private JFrame jf;
    private int time;
    
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
        if (agent instanceof BeccaAgent) {        
            ap = new AgentPanel((BeccaAgent)a);
            jf = AgentPanel.window(ap, true);
        }

        new DynamicChart(a.getClass().getName()) {

            @Override
            public double getReward() {
                double r = rewardTotal;
                rewardTotal = 0;
                return r;
            }

            @Override
            public double[] getAction() {
                return agent.getAction();
            }

            @Override
            public double[] getSensor() {
                return agent.getSensor();
            }

            @Override
            public double getTime() {
                return time;
            }
            
            
        };
        
    }
    private double reward, rewardTotal;
    
    long cycleDelayMS = 50;
    long displayPeriodMS = 10;
    long lastDisplay = -1;
    
    public Simulation(Class<? extends Agent> agentClass, World world) throws Exception {
        
        this.agent = agentClass.newInstance();
        agent.init(world);

        
        
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
        
        time = 0;
        
        while (world.isActive()) {
            /*    
            # Repeat the loop through the duration of the existence of the world 
                sensor, reward = world.step(action)
                world.visualize(agent)
                action = agent.step(sensor, reward)
            return agent.report_performance()
            */

            reward = world.step(agent.getAction(), agent.getSensor());
            rewardTotal+=reward;
            
            agent.step(reward);

            long now = System.currentTimeMillis();
            if ((lastDisplay == -1) || (now - lastDisplay > displayPeriodMS)) {
                if (ap!=null) ap.update();
                if (jf!=null) jf.setTitle("Reward: " + reward + ", @" + time);
                lastDisplay = System.currentTimeMillis();
            }
            
            if (cycleDelayMS > 0) {
                try {
                    Thread.sleep(cycleDelayMS);
                } catch (InterruptedException ex) {
                }
            }
            if (time % 1000==0) System.out.println(time);
            
            time++;
        }
        
        System.out.println("Simulation Finished.");
                          
    }
}
