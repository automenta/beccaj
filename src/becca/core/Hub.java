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
import org.ejml.data.DenseMatrix64F;

import static becca.core.Util.*;
import java.util.LinkedList;

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

    private double oldReward;
    public final LinkedList<Double> rewardTrace;

    private DenseMatrix64F count;
    private DenseMatrix64F expectedReward;
    private DenseMatrix64F cableActivities;
    private final LinkedList<DenseMatrix64F> pre;
    private final LinkedList<DenseMatrix64F> post;
    private DenseMatrix64F chainActivities;
    private DenseMatrix64F estimatedRewardValue;
    private final double[] rewardTraceArray;
    final private ArrayList<Integer> potentialWinners = new ArrayList();
    private DenseMatrix64F rewardUncertainty;
    private DenseMatrix64F activatedHubCableGoal;

    public Hub(int initialNumCables) {
        this.numCables = initialNumCables;

        // Set constants that adjust the behavior of the hub
        this.INITIAL_REWARD = BeccaParams.hubInitialReward;
        this.UPDATE_RATE = BeccaParams.hubUpdateRate;
        this.REWARD_DECAY_RATE = BeccaParams.hubRewardDecayRate;
        this.FORGETTING_RATE = BeccaParams.hubForgettingRate;
        this.TRACE_LENGTH = BeccaParams.hubTraceLength;
        this.EXPLORATION = BeccaParams.hubExploration;

        //# Initialize variables for later use
        this.rewardMin = Util.BIG;
        this.rewardMax = -Util.BIG;
        this.oldReward = 0.0;

        this.count = new DenseMatrix64F(this.numCables, this.numCables); // np.zeros((this.num_cables, this.num_cables))
        //this.count = new DenseMatrix(this.numCables, this.numCables);

        this.rewardTrace = new LinkedList();
        for (int i = 0; i < TRACE_LENGTH; i++) {
            rewardTrace.add(0.0);
        }
        rewardTraceArray = new double[rewardTrace.size()];

        this.expectedReward = new DenseMatrix64F(this.numCables, this.numCables);
        fill(expectedReward, INITIAL_REWARD);

        this.cableActivities = new DenseMatrix64F(numCables, 1);
        this.chainActivities = new DenseMatrix64F(numCables, numCables);
        this.estimatedRewardValue = new DenseMatrix64F(numCables, numCables);
        this.rewardUncertainty = new DenseMatrix64F(numCables, numCables);

        /*# pre represents the feature and sensor activities at a given
         # time step.
         # post represents the goal or action that was taken following. */
        this.pre = new LinkedList();
        this.post = new LinkedList();
        for (int i = 0; i < this.TRACE_LENGTH; i++) {
            /*this.pre = [np.zeros((this.num_cables, 1))] * (this.TRACE_LENGTH) 
             this.post = [np.zeros((this.num_cables, 1))] * (this.TRACE_LENGTH)*/
            pre.add(new DenseMatrix64F(this.numCables, 1));
            post.add(new DenseMatrix64F(this.numCables, 1));
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

        //reallocate matrices
        expectedReward = pad(expectedReward, numCables, numCables, INITIAL_REWARD);
        cableActivities = pad(cableActivities, numCables, 1, 0.0);
        chainActivities = pad(chainActivities, numCables, numCables, 0.0);
        estimatedRewardValue = pad(estimatedRewardValue, numCables, numCables, 0.0);
        
        count = pad(count, numCables, numCables, 0.0);
        rewardUncertainty = pad(rewardUncertainty, numCables, numCables, 0.0);

        //# All the cable activities from all the blocks, at the current time
        for (int i = 0; i < pre.size(); i++) {
            pre.set(i, pad(pre.get(i), numCables, 1, 0.0));
            post.set(i, pad(post.get(i), numCables, 1, 0.0));
        }
    }

    public void step(ArrayList<Block> blocks, double unscaledReward) {
        /*
         """ Advance the hub one step:
         1. Comb tower of blocks, collecting cable activities from each
         2. Update all-to-all reward model
         3. Select a goal
         4. Modify the goal in the block
         """
         */

        //# Adapt the reward so that it falls between -1 and 1 
        rewardMin = Math.min(unscaledReward, rewardMin);
        rewardMax = Math.max(unscaledReward, rewardMax);
        double spread = rewardMax - rewardMin;
        double newReward = ((unscaledReward - rewardMin) / (spread + EPSILON));
        rewardMin += spread * FORGETTING_RATE;
        rewardMax -= spread * FORGETTING_RATE;

        //# Use change in reward, rather than absolute reward
        double deltaReward = newReward - oldReward;
        oldReward = newReward;

        //# Update the reward trace, a brief history of reward
        rewardTrace.addLast(deltaReward);
        rewardTrace.removeFirst();

        //# Gather the cable activities from all the blocks
        int cableIndex = 0;

        if (blocks != null) {
            for (Block b : blocks) {
                int blockSize = b.getCableActivities().getNumElements();
                final double[] cad = cableActivities.getData();
                System.arraycopy(b.getCableActivities().getData(), 0, cableActivities.getData(), cableIndex, blockSize);
                cableIndex += blockSize;
            }
        } else {
            //use the data in cableActivities directly
        }

        /*
         # Update the reward model.
         # It has a structure similar to the chain transtion model 
         # in daisychain.
         # pre is composed of all the cable activities.
         # post is the selected goal that followed.
         */
        //this.chain_activities = this.pre[0] * this.post[0].T
        DenseMatrix64F firstPre = pre.peekFirst();
        DenseMatrix64F firstPost = post.peekFirst();

        assert (firstPre.getNumRows() == count.getNumRows());
        assert (firstPre.getNumRows() == firstPost.getNumRows());
        assert (firstPre.getNumCols() == 1);
        assert (firstPost.getNumCols() == 1);
        mult(firstPre, transpose(firstPost, null), chainActivities);

        //self.count = self.count + self.chain_activities        
        addEquals(count, chainActivities);

        //# Decay the count gradually to encourage occasional re-exploration 
        //this.count *= 1. - this.FORGETTING_RATE
        scale(1.0 - FORGETTING_RATE, count);

        //this.count = np.maximum(this.count, 0)
        matrixMaximum(count, 0);

        //# Calculate the rate at which to update the reward estimate
        //update_rate_raw = (this.chain_activities * ((1 - this.UPDATE_RATE) / (this.count + tools.EPSILON) + this.UPDATE_RATE))
        final DenseMatrix64F updateRateRawFactor = count.copy();
        add(updateRateRawFactor, EPSILON);
        matrixPower(updateRateRawFactor, -1.0);
        scale((1.0 - UPDATE_RATE), updateRateRawFactor);
        add(updateRateRawFactor, UPDATE_RATE);
        final DenseMatrix64F updateRate = new DenseMatrix64F(chainActivities.getNumRows(), updateRateRawFactor.getNumCols());
        elementMult(chainActivities, updateRateRawFactor, updateRate);

        //update_rate = np.minimum(0.5, update_rate_raw)
        matrixMinimum(updateRate, 0.5);

        /*
         # Collapse the reward history into a single value for this time step
         reward_array = np.array(this.reward_trace)
         */
        int ra = 0;
        for (Double d : rewardTrace) {
            rewardTraceArray[ra++] = d;
        }

        /*
         # TODO: substitute np.arange in this statement
         decay_exponents = (1. - this.REWARD_DECAY_RATE) ** (np.cumsum(np.ones(this.TRACE_LENGTH)) - 1.)
         */
        double[] decayExponents = new double[TRACE_LENGTH];
        for (int i = 0; i < TRACE_LENGTH; i++) {
            decayExponents[i] = Math.pow((1.0 - REWARD_DECAY_RATE), i);
        }

        //decayed_array = reward_array.ravel() * decay_exponents        
        //reward = np.sum(decayed_array.ravel())
        double reward = 0;
        for (int i = 0; i < rewardTraceArray.length; i++) {
            rewardTraceArray[i] *= decayExponents[i];
            reward += rewardTraceArray[i];
        }

        //System.out.println("ER: " + elementSum(expectedReward) );        
        //reward_difference = reward - this.expected_reward
        DenseMatrix64F rewardDifference = expectedReward.copy();
        scale(-1, rewardDifference);
        add(rewardDifference, reward);

        //System.out.println("ER: " + elementSum(expectedReward) + " RD:" + elementSum(rewardDifference) + " " + reward);
        //this.expected_reward += reward_difference * update_rate
        //DenseMatrix64F erDelta = new DenseMatrix64F(rewardDifference.getNumRows(), updateRate.getNumCols());
        //mult(rewardDifference, updateRate, erDelta);
        DenseMatrix64F erDelta = multMatrixMatrix(rewardDifference, updateRate);
        addEquals(expectedReward, erDelta);

        //System.out.println("ERDelta: " + elementSum(erDelta));
        //# Decay the reward value gradually to encourage re-exploration 
        //this.expected_reward *= 1. - this.FORGETTING_RATE
        scale(1.0 - FORGETTING_RATE, expectedReward);

        /*
         # Use the count to estimate the uncertainty in the expected 
         # value of the reward estimate.
         # Use this to scale additive random noise to the reward estimate,
         # encouraging exploration.
         reward_uncertainty = (np.random.normal(size=this.count.shape) * this.EXPLORATION / (this.count + 1.))
         */
        DenseMatrix64F uncertainNumerator = normRandMatrix(count.getNumRows(), count.getNumCols(), EXPLORATION, 0.0);
        rewardUncertainty.set(count);
        add(rewardUncertainty, 1.0);
        matrixDivBy(rewardUncertainty, uncertainNumerator);

        //this.estimated_reward_value = this.expected_reward + reward_uncertainty
        this.estimatedRewardValue.set(expectedReward);
        addEquals(estimatedRewardValue, rewardUncertainty);

        /*        
         # Select a goal cable.
         # First find the estimated reward associated with each chain.
         chain_votes = (this.cable_activities * this.estimated_reward_value + tools.EPSILON)
         */
        DenseMatrix64F chainVotes = matrixVector(estimatedRewardValue, cableActivities);
        add(chainVotes, EPSILON);

        /*
         # Find the maximum estimated reward associated with each potential goal
         hi_end = np.max(chain_votes, axis=0)
         # And the minimum estimated reward associated with each potential goal
         lo_end = np.min(chain_votes, axis=0)
         */
        DenseMatrix64F hiEnd = maxRow(chainVotes);
        DenseMatrix64F loEnd = minRow(chainVotes);

        /*
         # Sum the maxes and mins to find the goal with the highest mid-range  
         goal_votes = hi_end + lo_end
         */
        DenseMatrix64F goalVotes = loEnd; //loEnd.copy();
        addEquals(goalVotes, hiEnd);

        //potential_winners = np.where(goal_votes == np.max(goal_votes))[0] 
        double maxGoalVote = elementMax(goalVotes);
        potentialWinners.clear();
        potentialWinners.ensureCapacity(goalVotes.getNumElements());
        double[] gvd = goalVotes.getData();
        for (int i = 0; i < gvd.length; i++) {
            if (gvd[i] == maxGoalVote) {
                potentialWinners.add(i);
            }
        }

        if ((goalVotes.getNumElements() == 0) || (potentialWinners.size() == 0)) {
            /*System.err.println("No goals");
             System.err.println("Chain Votes: " + chainVotes);
             System.err.println("Estimated Reward Value: " + estimatedRewardValue);
             System.err.println("Reward Uncertainty: " + rewardUncertainty);
            
             System.err.println("Update Rate: " + updateRate);
             System.err.println("Reward: " + reward);
             System.err.println("Reward Trace: " + rewardTrace);
             System.err.println("Unscaled Reward: " + unscaledReward);
             System.err.println("Reward Difference: " + rewardDifference);
             System.err.println("Expected Reward: " + expectedReward);
            
             System.err.println("ERDelta: " + erDelta);*/

            //System.exit(1);
            return;
        }

        //# Break any ties by lottery
        /*
         if potential_winners.size < 1:
         print 'npw', potential_winners.size
         print 'max', np.max(goal_votes)
         winner = 0
         else:
         winner = potential_winners[np.random.randint(
         potential_winners.size)]        
         */
        final int winner;
        if (potentialWinners.isEmpty()) {
            winner = 0;
        } else {
            int pwi = (int) (Math.random() * potentialWinners.size());
            winner = potentialWinners.get(pwi);
        }

        //# Figure out which block the goal cable belongs to 
        //goal_cable = np.remainder(winner, this.cable_activities.size)
        //cable_index = goal_cable
        int goalCable = winner % cableActivities.getNumElements();
        cableIndex = goalCable;

        if (blocks != null) {
            for (Block b : blocks) {
                DenseMatrix64F h = b.getHubCableGoals();
                activatedHubCableGoal = h;
                int block_size = h.getNumElements();
                if (cableIndex >= block_size) {
                    cableIndex -= block_size;
                    continue;
                } else {
                    activateGoal(h, cableIndex, goalCable);
                    return;
                }
            }
        } else {
            if (activatedHubCableGoal == null)
                activatedHubCableGoal = new DenseMatrix64F(numCables, 1);
            activateGoal(activatedHubCableGoal, cableIndex, goalCable);            
        }

        //System.err.println("No goal chosen");
    }

    public DenseMatrix64F getOutput() {
        return activatedHubCableGoal;
    }
    
    protected void activateGoal(DenseMatrix64F h, int cableIndex, int goalCable) {
                    //# Activate the goal
        //block.hub_cable_goals[cable_index] = 1.
        assert (h.getNumCols() == 1);
        h.set(cableIndex, 0, 1.0);

        //new_post  = np.zeros(self.post[0].shape)
        DenseMatrix64F newPost = new DenseMatrix64F(post.peekFirst().getNumRows(), post.peekFirst().getNumCols());
        assert (newPost.getNumCols() == 1);
        newPost.set(goalCable, 0, 1);

                    //# Remove deliberate goals and actions from pre
        //new_pre = np.maximum(self.cable_activities.copy() - self.post[-1].copy(), 0.)
        DenseMatrix64F newPre = cableActivities.copy();
        DenseMatrix64F newPreDiff = post.peekLast().copy();
        subEquals(newPre, newPreDiff);
        matrixMaximum(newPre, 0);

        /*                        
         //# Update pre and post
         self.pre.append(new_pre)
         self.pre.pop(0)
         self.post.append(new_post)
         self.post.pop(0)
         */
        pre.removeFirst();
        pre.addLast(newPre);
        post.removeFirst();
        post.addLast(newPost);

        //System.err.println("Goal");
        //self._display()
    }

    public DenseMatrix64F getCount() {
        return count;
    }

    public DenseMatrix64F getCableActivities() {
        return cableActivities;
    }

    public DenseMatrix64F getChainActivities() {
        return chainActivities;
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
    public double[] getRewardArray() {
        return rewardTraceArray;
    }

    public DenseMatrix64F getExpectedReward() {
        return expectedReward;
    }

    public DenseMatrix64F getEstimatedRewardValue() {
        return estimatedRewardValue;
    }

    public LinkedList<DenseMatrix64F> getPre() {
        return pre;
    }

    public LinkedList<DenseMatrix64F> getPost() {
        return post;
    }

}
