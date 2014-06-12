package becca.world;

import becca.core.Simulation;
import becca.core.World;
import java.util.Arrays;

/*
    One-dimensional grid task
    
    In this task, the agent steps forward and backward along 
    a nine-position line. The fourth position is rewarded and 
    the ninth position is punished. There is also a slight 
    punishment for effort expended in trying to move, 
    i.e. taking actions. This is intended to be a simple-as-possible 
    task for troubleshooting BECCA. See Chapter 4 of the 
    Users Guide for details.
    Optimal performance is a reward of about 90 per time step.
*/
public class Grid1D implements World {
    private final int size;
    private final double VISUALIZE_PERIOD;
    private final double REWARD_MAGNITUDE;
    private final double ENERGY_COST;
    private final double JUMP_FRACTION;
    private double worldState;
    private int simpleState;
    private double energy;
    private double[] action;
    private double[] sensor;
    private final int totalTime;
    
    private int time;
    private final double noise;

    public Grid1D(int size, int totalTime, double noise) {
        this.time = 0;
        this.size = size;
        this.VISUALIZE_PERIOD = Math.pow(10, 4);
        this.REWARD_MAGNITUDE = 100.0;
        this.ENERGY_COST =  this.REWARD_MAGNITUDE / 100.0;
        this.JUMP_FRACTION = 0.1;        
        this.noise = noise;
        
        //this.name_long = 'one dimensional grid world'
        
        
        this.worldState = 0.0;
        this.simpleState = 0;
        this.totalTime = totalTime;
        /*
        this.display_state = True;
        */
    }
    
    @Override    public String getName()    {     return "Grid1D";    }
    @Override    public int getNumSensors() {     return size;    }
    @Override    public int getNumActions() {     return size;    }
    @Override    public boolean isActive()  {     return time < totalTime;   }

    @Override
    public double step(double[] action, double[] sensor) {

        time++;
        
        this.action = action;
        this.sensor = sensor;
        
        //# Find the step size as combinations of the action commands
        double stepSize = (    action[0] + 
                            2 * action[1] + 
                            3 * action[2] + 
                            4 * action[3] - 
                                action[4] - 
                            2 * action[5] - 
                            3 * action[6] - 
                            4 * action[7]);

        //# Action cost is an approximation of metabolic energy
        energy=(action[0] + 
                 2 * action[1] + 
                 3 * action[2] + 
                 4 * action[3] + 
                     action[4] + 
                 2 * action[5] + 
                 3 * action[6] + 
                 4 * action[7]);
        
        worldState = worldState + stepSize;  
                
        //# At random intervals, jump to a random position in the world
        if (Math.random() < JUMP_FRACTION) {
            worldState = size * Math.random();
        }
        
        //# Ensure that the world state falls between 0 and 9
        worldState -= size * Math.floor( ((double)worldState) / ((double)size) );
        simpleState = (int)Math.floor(worldState);
        if (simpleState == 9) simpleState = 0;
        
        /*        
        # Assign basic_feature_input elements as binary. 
        # Represent the presence or absence of the current position in the bin.
        */
        //Arrays.fill(sensor, 0);
        for (int i = 0; i < sensor.length; i++) {
            sensor[i] = (Math.random()*noise);
        }            
        sensor[simpleState] = 1 - (Math.random()*noise);
        
        System.out.println(toString());
        
        return getReward(sensor);
    }
    
    public double getReward(double[] sensor) {
        
        double reward = 0.;
        reward -= sensor[8] * REWARD_MAGNITUDE;
        reward += sensor[3] * REWARD_MAGNITUDE;
        
        //# Punish actions just a little
        reward -= energy  * ENERGY_COST;
        reward = Math.max(reward, -REWARD_MAGNITUDE);
        
        return reward;    
    }
    

    @Override
    public String toString() {
        String s = "";
        for (int i = 0; i < size; i++) {
            char c;
            if (i == simpleState)
                c = 'O';
            else
                c = '.';
            s += c;
        }
        s += "\n";
        for (int i = 0; i < size; i++) {
            char c;
            if (action[i] > 0)
                c = 'X';
            else
                c = '.';
            s += c;
        }
        s += "\n";
        return s;
    }
    
    /*
    def visualize(self, agent):
        """ Show what's going on in the world """
        if (this.display_state):
            state_image = ['.'] * (this.num_sensors + this.num_actions + 2)
            state_image[this.simple_state] = 'O'
            state_image[this.num_sensors:this.num_sensors + 2] = '||'
            action_index = np.where(this.action > 0.1)[0]
            if action_index.size > 0:
                for i in range(action_index.size):
                    state_image[this.num_sensors + 2 + action_index[i]] = 'x'
            print(''.join(state_image))
            
        if (this.timestep % this.VISUALIZE_PERIOD) == 0:
            print("world age is %s timesteps " % this.timestep)    
    */

    /*
    
    def set_agent_parameters(self, agent):
    """ Turn a few of the knobs to adjust BECCA for this world """
    # Prevent the agent from forming any groups
    #agent.reward_min = -100.
    #agent.reward_max = 100.
    pass
     */
    
    public static void main(String[] args) {
        new Simulation(new Grid1D(9, 5000, 0.1));
    }
}
