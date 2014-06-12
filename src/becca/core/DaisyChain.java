package becca.core;

import org.ejml.data.DenseMatrix64F;

import static becca.core.Util.*;

/**
    An incremental model-based reinforcement learning algorithm
    
    Input channels are chained together in two-element sequences,
    and an expected reward is associated with each.
    A helpful metaphor is joining cables end-to-end in a daisy chain,
    each containing a pre cable and a post cable.
    If activity in the post cable follows activity in pre cable, 
    the activity in the chain they share is high.

    The daisychain is the basis for prediction in the agent. 
    It incrementally estimates the likelihood that a post cable 
    will be active given the set of pre cables that have recently
    been active. 

    Each pre-post chain also tracks the expected reward, the uncertainty in
    the estimates of the post activity and reward, and a count of how many
    times the chain has been active.
 */
public class DaisyChain {
    private final double AGING_TIME_CONSTANT;
    private final double CHAIN_UPDATE_RATE;
    private final int maxCables;
    private final int time;
    private final DenseMatrix64F count;
    private final DenseMatrix64F expectedCableActivities;
    private final DenseMatrix64F postUncertainty;
    private DenseMatrix64F pre;
    private final DenseMatrix64F preCount;
    private DenseMatrix64F post;
    private int numCables;
    private DenseMatrix64F surprise;
    private DenseMatrix64F reaction;

    public DaisyChain(int maxCables) {
        
        this.maxCables = maxCables;

        this.AGING_TIME_CONSTANT = Math.pow(10, 6); //# real, large
        this.CHAIN_UPDATE_RATE = Math.pow(10,-1); // # real, 0 < x < 1
                
        this.time = 0;
        
        //this.shape = (max_num_cables, max_num_cables)        
        this.count = new DenseMatrix64F(maxCables, maxCables);        
        this.expectedCableActivities = new DenseMatrix64F(maxCables, maxCables);        
        this.postUncertainty = new DenseMatrix64F(maxCables, maxCables);
        
        //state_shape = (max_num_cables,1)
        this.pre = new DenseMatrix64F(maxCables, 1);
        this.preCount = new DenseMatrix64F(maxCables, 1);
        this.post = new DenseMatrix64F(maxCables, 1);
                
        this.numCables = 0;
        
        //#this.deliberation_vote = np.zeros((max_num_cables, 1))
        
        this.surprise = new DenseMatrix64F(maxCables, 1); //np.ones((max_num_cables, 1))
        fill(surprise, 1.0);
        
    }

    public DenseMatrix64F getSurprise() { 
        /*
        def get_surprise(self):
            return self.surprise[:self.num_cables]
        */
        
        if (surprise.getNumRows() > numCables) {
            if (numCables == 0)
                return new DenseMatrix64F(0, 1);
            //try {
                return extract(surprise, 0, numCables, 0, 1);
            /*}
            catch (Exception e) {
                System.err.println("DaisyChain getSurprise() trying to extract " + numCables + " columns from " + m(surprise));
                System.exit(1);
            }*/
        }
        return surprise;
    }
    
    public int getNumCables() {
        return numCables;
    }

    public DenseMatrix64F stepUp(DenseMatrix64F cableActivities) {
        //""" Train the daisychain using the current cable_activities """
    
        
        //self.num_cables = np.maximum(self.num_cables, cable_activities.size)        
        numCables = Math.max(numCables, cableActivities.getNumElements());
        
        //# Pad the incoming cable_activities array out to its full size 
        //cable_activities = tools.pad(cable_activities, (self.max_num_cables, 0))
        cableActivities = pad(cableActivities, maxCables, 0, 0.0);
        
        
        /* self.pre = self.post.copy()
           self.post = cable_activities.copy() */
        pre = post.copy();
        post = cableActivities.copy();
                
        DenseMatrix64F postT = transpose(post, null);

        //chain_activities = self.pre * self.post.T
        DenseMatrix64F chainActivities = new DenseMatrix64F(pre.getNumRows(), postT.getNumCols());
        mult(pre, postT, chainActivities);
        

        
        
        //set the main diagonal (of size pre) of chainActivities to zero
        //chain_activities[np.nonzero(np.eye(self.pre.size))] = 0.
        DenseMatrix64F eye = identity(pre.getNumElements());
        scale(-1, eye);
        add(eye, 1);
        elementMult(chainActivities, eye);
        
        //self.count += chain_activities
        addEquals(count, chainActivities);
        
        //self.count -= 1 / (self.AGING_TIME_CONSTANT * self.count + tools.EPSILON)
        DenseMatrix64F countDelta = count.copy();
        scale(AGING_TIME_CONSTANT, countDelta);
        add(countDelta, EPSILON);
        matrixPower(countDelta, -1);
        subEquals(count, countDelta);
        
        //self.count = np.maximum(self.count, 0)
        matrixMaximum(count, 0);
                
        /*        
        update_rate_raw_post = (self.pre * ((1 - self.CHAIN_UPDATE_RATE) /
                                   (self.pre_count + tools.EPSILON) +
                                    self.CHAIN_UPDATE_RATE))         
        */        
        DenseMatrix64F updateRateRawPost = preCount.copy(); add(updateRateRawPost, EPSILON);
        matrixPower(updateRateRawPost, -1);
        scale(1 - CHAIN_UPDATE_RATE, updateRateRawPost);
        add(updateRateRawPost, CHAIN_UPDATE_RATE);
        elementMult(updateRateRawPost, pre);
        
        //update_rate_post = np.minimum(0.5, update_rate_raw_post)
        DenseMatrix64F updateRatePost = updateRateRawPost.copy();
        matrixMinimum(updateRatePost, 0.5);
        
        //self.pre_count += self.pre
        addEquals(preCount, pre);
                        
        //self.pre_count -= 1 / (self.AGING_TIME_CONSTANT * self.pre_count + tools.EPSILON)
        DenseMatrix64F preCountDelta = preCount.copy();
        scale(AGING_TIME_CONSTANT, preCountDelta);
        add(preCountDelta, EPSILON);
        matrixPower(preCountDelta, -1);
        subEquals(preCount, preCountDelta);
       
        //self.pre_count = np.maximum(self.pre_count, 0)
        matrixMaximum(preCount, 0);
        
        //post_difference = np.abs(self.pre * self.post.T - self.expected_cable_activities)
        
        DenseMatrix64F postDifference = new DenseMatrix64F(pre.getNumRows(), postT.getNumCols());
        mult(pre, postT, postDifference);
        subEquals(postDifference, expectedCableActivities);
        matrixAbs(postDifference);
        
        
        //self.expected_cable_activities += update_rate_post * (self.pre * self.post.T - self.expected_cable_activities)
        DenseMatrix64F ecDelta = new DenseMatrix64F(pre.getNumRows(), postT.getNumCols());
        mult(pre, postT, ecDelta);
        subEquals(ecDelta, expectedCableActivities);        
        matrixVector(ecDelta, updateRatePost);        
        addEquals(expectedCableActivities, ecDelta);
        

        
        //self.post_uncertainty += (post_difference - self.post_uncertainty) * update_rate_post
        DenseMatrix64F puDelta = postDifference.copy();
        subEquals(puDelta, postUncertainty);        
        matrixVector(puDelta, updateRatePost);
        addEquals(postUncertainty, puDelta);
                
        //# Reaction is the expected post, turned into a deliberation_vote
        //self.reaction = tools.weighted_average(self.expected_cable_activities, self.pre)        
            //daisychain stepup reaction param (8, 8) (8, 1)
        reaction = getWeightedAverage(expectedCableActivities, pre);
        assert(reaction.getNumRows() == expectedCableActivities.getNumRows());
        

        //# Surprise is the difference between the expected post and the actual one
        //self.surprise = tools.weighted_average(
        //                  np.abs(self.post.T - self.expected_cable_activities),
        //                  self.pre / (self.post_uncertainty + tools.EPSILON))
        DenseMatrix64F sa = broadcastCols(postT, expectedCableActivities.getNumCols());
        subEquals(sa, expectedCableActivities);
        matrixAbs(sa);

        DenseMatrix64F sb = postUncertainty.copy();
        add(sb, EPSILON);

        matrixDivBy(sb, broadcastRows(pre, sb.getNumRows()));        
        
        surprise = getWeightedAverage(sa, sb);
        
        //# Reshape chain activities into a single column
        //return chain_activities.ravel()[:,np.newaxis]
        //check the order: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ravel.html?highlight=ravel
        return DenseMatrix64F.wrap(chainActivities.getNumElements(), 1, chainActivities.getData());
    }
    
    public DenseMatrix64F stepDown(DenseMatrix64F chainGoals) {
        //""" Propogate goals down through the transition model """
        
        
        //# Reshape chain_goals back into a square array 
        //chain_goals = np.reshape(chain_goals, (self.post.size, -1))
        chainGoals = DenseMatrix64F.wrap(
                post.getNumElements(), chainGoals.getNumElements() / post.getNumElements(),
                chainGoals.getData()
        );
                        
        //# Weight chain goals by the current cable activities   
        //upstream_goals = tools.bounded_sum(self.post * chain_goals.T)
                
        DenseMatrix64F upstreamGoals = matrixVector(transpose(chainGoals, null), post);
                
        upstreamGoals = boundedRowSum(upstreamGoals);
        
        //cable_goals = tools.bounded_sum([upstream_goals, self.reaction])
        DenseMatrix64F cableGoals = 
                DenseMatrix64F.wrap(upstreamGoals.getNumRows(), 1,
                        boundedSum(0, upstreamGoals.getData(), reaction.getData()));
        
        //return cable_goals[:self.num_cables]        
        if (numCables == 0) 
            return new DenseMatrix64F(cableGoals.getNumRows(), 0);
        
        //try {
            return extract(cableGoals, 0, numCables, 0, 1);
        /*}
        catch (Exception e) {
            System.err.println("DaisyChain stepDown() trying to extract " + numCables + " columns from " + m(cableGoals));
            System.exit(1);
            return null;
        }*/
        
    }
    
    /*    


    def get_index_projection(self, map_projection):
        """ Find the projection from chain activities to cable signals """
        num_cables = np.int(map_projection.size ** .5)
        projection = np.zeros((num_cables,2))
        chains = np.reshape(map_projection, (num_cables,num_cables))
        projection[:,0] = np.sign(np.max(chains, axis=1))
        projection[:,1] = np.sign(np.max(chains, axis=0))
        return projection
    
    def visualize(self, save_eps=True):
        """ Show the internal state of the daisychain in a pictorial format """
        tools.visualize_array(self.reward_value, 
                                  label=self.name + '_reward')
        #tools.visualize_array(self.reward_uncertainty, 
        #                          label=self.name + '_reward_uncertainty')
        tools.visualize_array(np.log(self.count + 1.), 
                                  label=self.name + '_count')
        #tools.visualize_daisychain(self, self.num_primitives, 
        #                          self.num_actions, 10)
        return
    */    
    public String toString() {
        StringBuffer sb = new StringBuffer(super.toString());
        sb.append("{ count=" + this.count + "}");
        return sb.toString();
    }

    

}
