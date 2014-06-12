package becca.core;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.ejml.data.BlockMatrix64F;
import org.ejml.data.DenseMatrix64F;

import static becca.core.Util.*;
       
/**
 * A general reinforcement learning agent
 *
 * Takes in a time series of sensory input vectors and a scalar reward and puts
 * out a time series of action commands.
 */
public class Agent implements Serializable {

    public final String name;
    public final ArrayList<Block> blocks = new ArrayList();
    private int time; //current time step

    private LinkedList<Double> recentSurpriseHistory;
    private final int RECENT_SURPRISE_HISTORY_SIZE = 100;

    private final int numSensors;
    private final int numActions;
    private int timeSinceRewardLog;
    private double cumulativeReward;
    private final ArrayList<Double> rewardHistory;
    private final LinkedList<Double> surpriseHistory;
    private final ArrayList<Object> rewardSteps;
    public final Hub hub;
    public final double[] sensor;
    public final double[] action;
    private double reward;
    private double typicalSurprise;

    /*
         Configure the Agent

         num_sensors and num_actions are the only absolutely necessary
         arguments. They define the number of elements in the 
         sensor and action arrays that the agent and the world use to
         communicate with each other. 
    */
    public Agent(String name, int numSensors, int numActions) {
        
        //self.BACKUP_PERIOD = 10 ** 4

        this.name = name;

        this.time = 0;

        //TODO: Automatically adapt to the number of sensor pass in
        this.numSensors = numSensors;
        this.sensor = new double[numSensors];
        this.numActions = numActions;
        this.action = new double[numActions];

        //first_block_name = ''.join(('block_', str(self.num_blocks - 1)))
        blocks.add(new Block(numActions + numSensors));
        this.hub = new Hub(blocks.get(0).maxCables);

        this.cumulativeReward = 0;
        this.timeSinceRewardLog = 0;
        this.rewardHistory = new ArrayList();
        this.rewardSteps = new ArrayList();
        this.surpriseHistory = new LinkedList();
        this.recentSurpriseHistory = new LinkedList();

        this.reward = 0;
    }

    
    public int step(double reward) {
        //""" Step through one time interval of the agent's operation """
        //reads current 'sensor' value and sets action in 'action'
        
        this.time++;
        this.reward = reward;

        //cable_activities = np.vstack((self.action, sensor))
        DenseMatrix64F cableActivities = new DenseMatrix64F(numActions + numSensors, 1);
        System.arraycopy(action, 0, cableActivities.getData(), 0, numActions);
        System.arraycopy(sensor, 0, cableActivities.getData(), numActions, numSensors);
        
        
        /*for (int i = 0; i < numActions; i++)
            cableActivitiesData[i][0] = action[i];
        for (int i = 0; i < numSensors; i++)
            cableActivitiesData[i][1] = sensor[i];*/
                        
        //# Propogate the new sensor inputs up through the blocks
        for (final Block b : blocks) {
            b.stepUp(cableActivities);
        }
        
        //# Create a new block if the top block has had enough bundles assigned
        Block topBlock = blocks.get(blocks.size()-1); //top block        
        double blockBundlesFull = ((double)(topBlock.getBundlesCreated())) / ((double)(topBlock.maxBundles));
        
        double blockInitializationThreshold = 0.5;
        
        if (blockBundlesFull > blockInitializationThreshold) {
            Block b = new Block(numActions + numSensors, blocks.size());
            blocks.add(b);
            cableActivities = b.stepUp(cableActivities);
            hub.addCables(b.maxCables);
            //print "Added block", self.num_blocks - 1
        }
        
        hub.step(blocks, reward);

        //# Propogate the deliberation_goal_votes down through the blocks
        double agentSurprise = 0.0;
        BlockMatrix64F cableGoals = new BlockMatrix64F(cableActivities.getNumCols() * cableActivities.getNumRows(), 1);

        System.out.println("cableGoals:" + cableGoals);

        //blocks in reverse
        for (int i = blocks.size()-1; i >=0; i--) {
            Block b = blocks.get(i);
            cableGoals = b.stepDown(cableGoals);
            
            /*
            if np.nonzero(block.surprise)[0].size > 0:
                agent_surprise = np.sum(block.surprise)            
            */
            BlockMatrix64F s = b.getSurprise();
            double[] sd = s.getData();
            int nonzeros = 0;
            for (int g = 0; g < sd.length; g++)
                nonzeros += sd[g] > 0 ? 1 : 0;            
            if (nonzeros > 0)
                agentSurprise = sum(s);
        }
        
        
        recentSurpriseHistory.pop();    //remove first element
        recentSurpriseHistory.add(agentSurprise);

        //self.typical_surprise = np.median(np.array(self.recent_surprise_history))

        Median m = new Median(); 
        double[] rsh = new double[recentSurpriseHistory.size()]; int i = 0; for (Double d : recentSurpriseHistory) rsh[i++] = d;
        this.typicalSurprise = m.evaluate(rsh);
        
        double modSurprise = agentSurprise - typicalSurprise;
        surpriseHistory.add(modSurprise);
        
        /*
        # Strip the actions off the deliberation_goal_votes to make 
        # the current set of actions.
        # For actions, each goal is a probability threshold. If a roll of
        # dice comes up lower than the goal value, the action is taken
        # with a magnitude of 1.
        */
        System.arraycopy(cableGoals.getData(), 0, action, 0, numActions);
        
        
        //backup?
        /*
        if (self.timestep % self.BACKUP_PERIOD) == 0:
                self._save()    
        */        
        
        //# Log reward
        this.cumulativeReward += reward;
        this.timeSinceRewardLog += 1;
        
        return 0;

    }

    public double[][] getIndexProjections() {
        /*  
         Get representations of all the bundles in each block 
        
         Every feature is projected down through its own block and
         the blocks below it until its cable_contributions on sensor inputs 
         and action is obtained. This is a way to represent the
         receptive field of each feature.

         Returns a list containing the cable_contributions for each feature 
         in each block.
         """
         all_projections = []
         all_bundle_activities = []
         for block_index in range(len(self.blocks)):
         block_projections = []
         block_bundle_activities = []
         num_bundles = self.blocks[block_index].max_bundles
         for bundle_index in range(num_bundles):    
         bundles = np.zeros((num_bundles, 1))
         bundles[bundle_index, 0] = 1.
         cable_contributions = self._get_index_projection(
         block_index,bundles)
         if np.nonzero(cable_contributions)[0].size > 0:
         block_projections.append(cable_contributions)
         block_bundle_activities.append(self.blocks[block_index].
         bundle_activities[bundle_index])
         # Display the cable_contributions in text form if desired
         if to_screen:
         print 'cable_contributions', \
         self.blocks[block_index].name, \
         'feature', bundle_index
         for i in range(cable_contributions.shape[1]):
         print np.nonzero(cable_contributions)[0][
         np.where(np.nonzero(
         cable_contributions)[1] == i)]
         if len(block_projections) > 0:
         all_projections.append(block_projections)
         all_bundle_activities.append(block_bundle_activities)
         return (all_projections, all_bundle_activities)
         */
        return null;
    }

    public double[][] getIndexProjection(int blockIndex, double[][] bundles) {
        /*
    
         """
         Get the cable_contributions for bundles
        
         Recursively project bundles down through blocks
         until the bottom block is reached. Feature values is a 
         two-dimensional array and can contain
         several columns. Each column represents a state, and their
         order represents a temporal progression. During cable_contributions
         to the next lowest block, the number of states
         increases by one. 
        
         Return the cable_contributions in terms of basic sensor 
         inputs and action. 
         """
         */
        if (blockIndex == -1) {
            return bundles;
        }

        /*
         cable_contributions = np.zeros((self.blocks[block_index].max_cables, 
         bundles.shape[1] + 1))
         for bundle_index in range(bundles.shape[0]):
         for time_index in range(bundles.shape[1]):
         if bundles[bundle_index, time_index] > 0:
         new_contribution = self.blocks[
         block_index].get_index_projection(bundle_index)
         cable_contributions[:,time_index:time_index + 2] = (
         np.maximum(
         cable_contributions[:,time_index:time_index + 2], 
         new_contribution))
         cable_contributions = self._get_index_projection(block_index - 1, 
         cable_contributions)
         return cable_contributions
         */
        return null;
    }

    public String toString() {
        StringBuffer x = new StringBuffer();
        x.append("Agent " + this.name + " (sensors=" + this.numSensors + ", actions=" + this.numActions + ") @ time=" + this.time + ":\n");
        x.append("  " + this.blocks + "\n");
        x.append("  " + this.hub);
        
        //x.append("Reward history: " +this.rewardHistory + "\n");
        
        return x.toString();
    }
    
    /*
     def visualize(self):
     """ Show the current state and some history of the agent """
     print ' '.join([self.name, 'is', str(self.timestep), 'time steps old'])
     self.reward_history.append(float(self.cumulative_reward) / 
     (self.time_since_reward_log + 1))
     self.cumulative_reward = 0    
     self.time_since_reward_log = 0
     self.reward_steps.append(self.timestep)
     self._show_reward_history()
     for block in self.blocks:
     block.visualize()
     pass
     return
 
     def report_performance(self):
     """ Report on the reward amassed by the agent """
     performance = np.mean(self.reward_history)
     print("Final performance is %f" % performance)
     self._show_reward_history(hold_plot=self.show)
     return performance
    
     def _show_reward_history(self, hold_plot=False, 
     filename='log/reward_history.png'):
     """ Show the agent's reward history and save it to a file """
     if self.graphing:
     fig = plt.figure(1)
     plt.plot(self.reward_steps, self.reward_history)
     plt.xlabel("time step")
     plt.ylabel("average reward")
     plt.title(''.join(('Reward history for ', self.name)))
     fig.show()
     fig.canvas.draw()
     plt.savefig(filename, format='png')
     if hold_plot:
     plt.show()
     return  
  
     */
}
