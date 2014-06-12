/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.core;



/**
 *
 * @author me
 */
public class Simulation {

    public final Agent agent;
    
    /*
    Run BECCA with world.  

    If restore=True, this method loads a saved agent if it can find one.
    Otherwise it creates a new one. It connects the agent and
    the world together and runs them for as long as the 
    world dictates.

    To profile BECCA's performance with world, manually set
    profile_flag in the top level script environment to True.

    */
    public Simulation(World world) {
        
        /*
        if agent_name is None:
            agent_name = '_'.join((world.name, 'agent'))
        */
        this.agent = new Agent(world.getName() + " Agent", 
                            world.getNumSensors(), world.getNumActions());
        
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

            double reward = world.step(agent.action, agent.sensor);

            agent.step(reward);

        }
                          
    }
}
