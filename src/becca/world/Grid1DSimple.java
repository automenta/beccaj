package becca.world;

import becca.core.Simulation;
import becca.core.World;

/*
    Simplest one-dimensional grid task
    
    In this task, the agent's goal is to activate the same pattern as
    it sees in its sensor field.
*/
public class Grid1DSimple implements World {
    private final int size;
    
    private final double VISUALIZE_PERIOD;
    private final double REWARD_MAGNITUDE;
    private final double JUMP_FRACTION;

    private double focusPosition;
    private double focusVelocity;
        
    private double[] action;

    private final int totalTime;    
    private int time;
    
    private final double noise;

    public Grid1DSimple(int size, int totalTime, double noise, double focusVelocity) {
        this.time = 0;
        this.size = size;
        this.VISUALIZE_PERIOD = Math.pow(10, 4);
        this.REWARD_MAGNITUDE = 100.0;
        this.JUMP_FRACTION = 0.0;        
        this.noise = noise;
        this.focusVelocity = focusVelocity;
        
        this.focusPosition = 0.0;
        this.totalTime = totalTime;
    }
    
    @Override    public String getName()    {     return "Grid1D";    }
    @Override    public int getNumSensors() {     return size;    }
    @Override    public int getNumActions() {     return size;    }
    @Override    public boolean isActive()  {     return time < totalTime;   }

    @Override
    public double step(double[] action, double[] sensor) {

        time++;
        
        this.action = action;
        

        
        focusPosition = focusPosition + focusVelocity;  
                
        //# At random intervals, jump to a random position in the world
        if (Math.random() < JUMP_FRACTION) {
            focusPosition = size * Math.random();
        }
        else {
            focusPosition += focusVelocity;
        }
        
        //# Ensure that the world state falls between 0 and 9
        focusPosition -= size * Math.floor( ((double)focusPosition) / ((double)size) );
        if (focusPosition > size) focusPosition = focusPosition - size;
        if (focusPosition < 0) focusPosition = size + focusPosition;
        
        /*        
        # Assign basic_feature_input elements as binary. 
        # Represent the presence or absence of the current position in the bin.
        */
        //Arrays.fill(sensor, 0);
        for (int i = 0; i < sensor.length; i++) {
            final double exp = 3.0; //sharpen
            sensor[i] = Math.pow(1.0 / (1.0 + Math.abs( ((double)i)-focusPosition)),exp) + (Math.random()*noise);
            if (sensor[i] < 0.2)
                sensor[i] = 0;
        }            
        
        double distance = 0;        
        for (int i = 0; i < sensor.length; i++)
            distance += Math.abs(sensor[i] - action[i]);
        
        double reward = REWARD_MAGNITUDE * (1.0 / (1.0 + distance)) - REWARD_MAGNITUDE/2.0;
        
        return reward;        
    }
        

    @Override
    public String toString() {
        String s = "";
        for (int i = 0; i < size; i++) {
            char c;
            if (i == (int)focusPosition)
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
        new Simulation(new Grid1DSimple(7, 50000, 0.002, 0.001));
    }
}
