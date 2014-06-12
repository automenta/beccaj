package becca.core;


/*
The hub is the central action selection mechanism 
    
    The analogy of the hub and spoke stucture is suggested by 
    the fact that the hub has a connection to each
    of the blocks. In the course of each timestep it 
    1) reads in a copy of the input cable activities to 
        each of the blocks
    2) updates its reward distribution estimate for all 
        the goals it could select
    3) selects a goal and
    4) declares that goal in the appropriate block.
*/
import java.util.ArrayList;
import org.ejml.data.BlockMatrix64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.data.Matrix64F;
import org.ejml.ops.CommonOps;

import static becca.core.Util.*;

public class Hub {
    
    private int numCables;
    private final double INITIAL_REWARD;
    private final double UPDATE_RATE;
    private final double REWARD_DECAY_RATE;
    private final double FORGETTING_RATE;
    private final int TRACE_LENGTH;
    private final double EXPLORATION;
    private double rewardMin;
    private double rewardMax;
    
    private final double oldReward;
    private final double[] rewardTrace;
    
    private DenseMatrix64F count;
    private DenseMatrix64F expectedReward;
    private DenseMatrix64F cableActivities;
    private DenseMatrix64F[] pre;
    private DenseMatrix64F[] post;

    public Hub(int initialNumCables) {
        this.numCables = initialNumCables;
        
        // Set constants that adjust the behavior of the hub
        this.INITIAL_REWARD = 1.0;
        this.UPDATE_RATE = Math.pow(10, -2);
        this.REWARD_DECAY_RATE = .3;
        this.FORGETTING_RATE = Math.pow(10, -5);
        this.TRACE_LENGTH = 10;
        this.EXPLORATION = .1;
        
        
        //# Initialize variables for later use
        this.rewardMin = Util.BIG;        
        this.rewardMax = -Util.BIG;
        this.oldReward = 0.0;
                
        this.count = new DenseMatrix64F(this.numCables, this.numCables); // np.zeros((this.num_cables, this.num_cables))
        this.rewardTrace = new double[this.TRACE_LENGTH];
        
        this.expectedReward = new DenseMatrix64F(this.numCables, this.numCables);
        CommonOps.fill(expectedReward, INITIAL_REWARD);
        
        this.cableActivities = new DenseMatrix64F(this.numCables, 1);
        
        /*# pre represents the feature and sensor activities at a given
          # time step.
          # post represents the goal or action that was taken following. */
        this.pre = new DenseMatrix64F[this.TRACE_LENGTH];
        this.post = new DenseMatrix64F[this.TRACE_LENGTH];
        for (int i = 0; i < this.TRACE_LENGTH; i++) {
            /*this.pre = [np.zeros((this.num_cables, 1))] * (this.TRACE_LENGTH) 
            this.post = [np.zeros((this.num_cables, 1))] * (this.TRACE_LENGTH)*/
            pre[i] = new DenseMatrix64F(this.numCables, 1);
            post[i] = new DenseMatrix64F(this.numCables, 1);
        }
                
    }
    
    @Override
    public String toString() {
        String s = super.toString();
        return s;
    }


    void addCables(int numNewCables) {    
        //""" Add new cables to the hub when new blocks are created """ 
        numCables += numNewCables;
        
        expectedReward = pad(expectedReward, numCables, numCables, INITIAL_REWARD);
        cableActivities = pad(cableActivities, numCables, 1, 0.0);
        count = pad(count, numCables, numCables, 0.0);

        
        //# All the cable activities from all the blocks, at the current time
        for (int i = 0; i < pre.length; i++) {
            pre[i] = pad(pre[i], numCables, 1, 0.0);
            post[i] = pad(post[i], numCables, 1, 0.0);            
        }
    
    }
        
    void step(ArrayList<Block> blocks, double unscaledReward) {
        /*
        """ Advance the hub one step:
        1. Comb tower of blocks, collecting cable activities from each
        2. Update all-to-all reward model
        3. Select a goal
        4. Modify the goal in the block
        """
        */
        /*
        # Adapt the reward so that it falls between -1 and 1 
        this.reward_min = np.minimum(unscaled_reward, this.reward_min)
        this.reward_max = np.maximum(unscaled_reward, this.reward_max)
        spread = this.reward_max - this.reward_min
        new_reward = ((unscaled_reward - this.reward_min) / 
                       (spread + tools.EPSILON))
        this.reward_min += spread * this.FORGETTING_RATE
        this.reward_max -= spread * this.FORGETTING_RATE

        # Use change in reward, rather than absolute reward
        delta_reward = new_reward - this.old_reward
        this.old_reward = new_reward
        # Update the reward trace, a brief history of reward
        this.reward_trace.append(delta_reward)
        this.reward_trace.pop(0)
        
        # Gather the cable activities from all the blocks
        cable_index = 0
        block_index = 0
        for block in blocks:
            block_size =  block.cable_activities.size
            this.cable_activities[cable_index: cable_index + block_size] = \
                    block.cable_activities.copy()
            cable_index += block_size 
            block_index += 1

        # Update the reward model.
        # It has a structure similar to the chain transtion model 
        # in daisychain.
        # pre is composed of all the cable activities.
        # post is the selected goal that followed.
        this.chain_activities = this.pre[0] * this.post[0].T
        # Update the count of how often each feature has been active
        this.count = this.count + this.chain_activities
        # Decay the count gradually to encourage occasional re-exploration 
        this.count *= 1. - this.FORGETTING_RATE
        this.count = np.maximum(this.count, 0)
        # Calculate the rate at which to update the reward estimate
        update_rate_raw = (this.chain_activities * ((1 - this.UPDATE_RATE) / 
                                               (this.count + tools.EPSILON) + 
		                                       this.UPDATE_RATE)) 
        update_rate = np.minimum(0.5, update_rate_raw)
        # Collapse the reward history into a single value for this time step
        reward_array = np.array(this.reward_trace)
        # TODO: substitute np.arange in this statement
        decay_exponents = (1. - this.REWARD_DECAY_RATE) ** (
                np.cumsum(np.ones(this.TRACE_LENGTH)) - 1.)
        decayed_array = reward_array.ravel() * decay_exponents
        reward = np.sum(decayed_array.ravel())
        reward_difference = reward - this.expected_reward 
        this.expected_reward += reward_difference * update_rate
        # Decay the reward value gradually to encourage re-exploration 
        this.expected_reward *= 1. - this.FORGETTING_RATE
        # Use the count to estimate the uncertainty in the expected 
        # value of the reward estimate.
        # Use this to scale additive random noise to the reward estimate,
        # encouraging exploration.
        reward_uncertainty = (np.random.normal(size=this.count.shape) *
                              this.EXPLORATION / (this.count + 1.))
        this.estimated_reward_value = this.expected_reward + reward_uncertainty

        # Select a goal cable.
        # First find the estimated reward associated with each chain.   
        chain_votes = (this.cable_activities * this.estimated_reward_value + 
                       tools.EPSILON)
        # Find the maximum estimated reward associated with each potential goal
        hi_end = np.max(chain_votes, axis=0)
        # And the minimum estimated reward associated with each potential goal
        lo_end = np.min(chain_votes, axis=0)
        # Sum the maxes and mins to find the goal with the highest mid-range  
        goal_votes = hi_end + lo_end
        potential_winners = np.where(goal_votes == np.max(goal_votes))[0] 
        # Break any ties by lottery
        winner = potential_winners[np.random.randint(potential_winners.size)]
        # Figure out which block the goal cable belongs to 
        goal_cable = np.remainder(winner, this.cable_activities.size)
        cable_index = goal_cable
        for block in blocks:
            block_size =  block.hub_cable_goals.size
            if cable_index >= block_size:
                cable_index -= block_size
                continue
            else:
                # Activate the goal
                block.hub_cable_goals[cable_index] = 1.
                new_post  = np.zeros(this.post[0].shape)
                new_post[goal_cable] = 1.
                # Remove deliberate goals and actions from pre
                new_pre = np.maximum(this.cable_activities.copy() - 
                                     this.post[-1].copy(), 0.)
                # Update pre and post
                this.pre.append(new_pre)
                this.pre.pop(0)
                this.post.append(new_post)
                this.post.pop(0)
                this._display()
                return
        print 'No goal chosen'
        return 
        */
    }
        

    /*
    def _display(self):
        """ Give a visual update of the internal workings of the hub """
        DISPLAY_PERIOD = 1000.
        #if np.random.random_sample() < 1. / DISPLAY_PERIOD:
        if False:

            # Plot reward value
            fig311 = plt.figure(311)
            plt.gray()
            plt.imshow(this.expected_reward, interpolation='nearest')
            plt.title('reward')
            fig311.show()
            fig311.canvas.draw()
            plt.savefig('log/reward_image.png', bbox_inches=0.)
            
            # Plot weighted chain votes
            fig313 = plt.figure(313)
            plt.gray()
            plt.imshow(np.maximum(this.cable_activities * 
                                  this.estimated_reward_value, 0.), 
                       interpolation='nearest')
            plt.title('cable activities * reward')
            fig313.show()
            fig313.canvas.draw()
            
            # Plot the count 
            fig314 = plt.figure(314)
            plt.gray()
            plt.imshow(1. / (this.count + 1.), interpolation='nearest')
            plt.title('1 / count')
            fig314.show()
            fig314.canvas.draw()
            
            # Plot the reward value plus exploration
            fig315 = plt.figure(315)
            plt.gray()
            plt.imshow(np.maximum(this.estimated_reward_value, 0.), 
                       interpolation='nearest')
            plt.title('estimated_reward_value')
            fig315.show()
            fig315.canvas.draw()
            plt.savefig('log/estimated_reward_image.png', bbox_inches=0.)
            
            # Plot the chain activities 
            fig316 = plt.figure(316)
            plt.gray()
            plt.imshow(this.chain_activities, interpolation='nearest')
            plt.title('chain_activities')
            fig316.show()
            fig316.canvas.draw()
    
    */


}
