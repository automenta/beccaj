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
     
Dasiychain is an incremental algorithm that estimates
the probability of one cable being active following an-
other. It represents this as a conditional probability:
given that one cable is active what is the expected ac-
tivity of a second cable in the next time step. (Figure
7) High expected chain activities indicate sequences of
cable activities that co-occur regularly. They identify
temporal structure in the data.

Expected chain activities are similar to transition prob-
abilities in Markov models. The difference is that in a
Markov model, only one state can be occupied at each
time step. This is analogous to just one cable being
active. In a daisychain, many cables can be completely
or partially active at once. As a result, transition prob-
abilities can sum to much more than one. 
 
A temporal sequence of one cable being active, followed
by another, is a chain. The activity of a chain is given
by the product of the two cable activities involved.

A leaky accumulation of the activity on each cable and
on each chain is also maintained.

The expected chain activities are maintained and up-
dated based on the current chain activities. 

In addition, the expected deviation from the expected
chain activities are maintained and updated based on
the difference between the current and expected chain
activities.

The temporal structure captured in the expected chain
activities provide a basis for making short-term predic-
tions.The reaction is the predicted next set of cable ac-
tivities.

The most recently observed cable activities can be compared
to those that would have been predicted from the previous 
cable activities to find surprising events.

As chain goals are propagated down through the daisy-
chain, they are combined to form the cable goals. Each
cable goal is a weighted, bounded sum of all the goals
of the chains it belongs to.

 */
public class DaisyChain {
    private final double countDecayRate;
    private final double chainUpdateRate;
    private final int maxCables;
    private final DenseMatrix64F count;
    private final DenseMatrix64F expectedCableActivities;
    private DenseMatrix64F pre;
    private final DenseMatrix64F preCount;
    private DenseMatrix64F post;
    private final DenseMatrix64F postUncertainty;
    private int numCables;
    private final DenseMatrix64F surprise;
    private DenseMatrix64F reaction;
    private final boolean allowSelfTransitions;
    private DenseMatrix64F chainActivities;
    private DenseMatrix64F nullGoals;
    private double agingTimeconstant = BeccaParams.daisyAgingTimeConstant;
    private DenseMatrix64F updateRatePost;
    private DenseMatrix64F eye;
    private DenseMatrix64F postDifference;
    private DenseMatrix64F ecDeltaFactor, puDelta;
    private DenseMatrix64F sb;
    private DenseMatrix64F expectedCableDelta;


    
    public DaisyChain(int maxCables) {
        
        this.maxCables = maxCables;

        this.countDecayRate = BeccaParams.daisyCountDecayRate;
        this.chainUpdateRate = BeccaParams.daisyChainUpdateRate;; // # real, 0 < x < 1         
        this.allowSelfTransitions = BeccaParams.daisyAllowSelfTransitions;
                
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
        
        this.surprise = new DenseMatrix64F(1, maxCables); //np.ones((max_num_cables, 1))
        this.reaction = new DenseMatrix64F(1, maxCables); //np.ones((max_num_cables, 1))
        fill(surprise, 1.0);
        
    }

    public DenseMatrix64F getSurprise() { 
        /*
        def get_surprise(self):
            return self.surprise[:self.num_cables]
        */

        /*
        if (surprise.getNumRows() > numCables) {
            if (numCables == 0) {
                surprise.reshape(0,1);
                //return new DenseMatrix64F(0, 1);
            }
            else {
                surprise.reshape(numCables, 1);
                //return extract(surprise, 0, numCables, 0, 1);
            }
        }*/
        /*if (numCables < surprise.elements) {
            Arrays.fill(surprise.getData(), numCables, surprise.elements, 0);
        }*/
        return surprise;
    }
    
    public int getNumCables() {
        return numCables;
    }

    public DenseMatrix64F stepUp(DenseMatrix64F cableActivities) {
        //""" Train the daisychain using the current cable_activities """
            
        //self.num_cables = np.maximum(self.num_cables, cable_activities.size)        
        numCables = Math.max(numCables, cableActivities.elements);
        
        //# Pad the incoming cable_activities array out to its full size 
        //cable_activities = tools.pad(cable_activities, (self.max_num_cables, 0))
        cableActivities = pad(cableActivities, maxCables, 0, 0.0);
        
        
        /* self.pre = self.post.copy()
           self.post = cable_activities.copy() */       
        pre.set(post);
        post.set(cableActivities);
        
                
        DenseMatrix64F postT = transpose(post, null);

        //chain_activities = self.pre * self.post.T
        chainActivities = ensureSize(chainActivities, pre.getNumRows(), postT.getNumCols());
        assert(pre.getNumCols() == 1);
        assert(postT.getNumRows() == 1);
        mult(pre, postT, chainActivities);
        

        
        
        //set the main diagonal (of size pre) of chainActivities to zero
        //chain_activities[np.nonzero(np.eye(self.pre.size))] = 0.
        if (!allowSelfTransitions) {
            if ((eye == null) || (eye.getNumRows()!=pre.elements)) {
                eye = identity(pre.elements);            
                scale(-1, eye);
                add(eye, 1);
            }
            elementMult(chainActivities, eye);
        }
        
        //self.count += chain_activities
        addEquals(count, chainActivities);
        
        
        //self.count -= 1 / (self.AGING_TIME_CONSTANT * self.count + tools.EPSILON)
        //self.count = np.maximum(self.count, 0)
        final double[] cod = count.getData();
        for (int i = 0; i < count.elements; i++) {
            double cd = cod[i];
            cd -= 1 / (agingTimeconstant * cd + EPSILON);            
            if (cd < 0) cd = 0.0;
            cod[i] = cd;
        }   
        
                
        /*        
        update_rate_raw_post = (self.pre * ((1 - self.CHAIN_UPDATE_RATE) /
                                   (self.pre_count + tools.EPSILON) +
                                    self.CHAIN_UPDATE_RATE))         
        */        
//        final DenseMatrix64F updateRateRawPost = preCount.copy(); add(updateRateRawPost, EPSILON);
//        matrixPower(updateRateRawPost, -1);
//        scale(1 - CHAIN_UPDATE_RATE, updateRateRawPost);
//        add(updateRateRawPost, CHAIN_UPDATE_RATE);
//        elementMult(updateRateRawPost, pre);                
//        
//        //update_rate_post = np.minimum(0.5, update_rate_raw_post)
//        final DenseMatrix64F updateRatePost = updateRateRawPost;
//        matrixMinimum(updateRatePost, 0.5);
        
        //ALTERNATE CALCULATION of updateRatePost:
        updateRatePost = ensureCopy(updateRatePost, preCount);
        final double[] updateRatePostD = updateRatePost.getData();
        for (int i = 0; i < updateRatePost.elements; i++) {
            final double u = preCount.getData()[i];
            final double v = 
                    Math.min(
                            (1-chainUpdateRate) * 
                            (1 + chainUpdateRate* (u + EPSILON)) /
                            (u + EPSILON), 0.5);      
                  /*update_rate_raw_post = (self.pre * ((1 - self.CHAIN_UPDATE_RATE) /
                                   (self.pre_count + tools.EPSILON) +
                                    self.CHAIN_UPDATE_RATE))     */
            
            updateRatePostD[i] = v;
        }
        
        
        //self.pre_count += self.pre
        addEquals(preCount, pre);
  
        
        
        //self.pre_count -= 1 / (self.AGING_TIME_CONSTANT * self.pre_count + tools.EPSILON)
        //self.pre_count = np.maximum(self.pre_count, 0)        
        for (int i = 0; i < preCount.elements; i++) {
            double pcd = preCount.getData()[i];
            pcd -= 1 / (agingTimeconstant * pcd + EPSILON);            
            if (pcd < 0) pcd = 0.0;
            preCount.getData()[i] = pcd;
        }

                       
        
        
        //post_difference = np.abs(self.pre * self.post.T - self.expected_cable_activities)        
        postDifference = ensureSize(postDifference, pre.getNumRows(), postT.getNumCols());
        mult(pre, postT, postDifference);        
        subEquals(postDifference, expectedCableActivities);


        //ecDelta = self.pre * self.post.T - self.expected_cable_activities
        ecDeltaFactor = ensureCopy(ecDeltaFactor, postDifference); //used below
        matrixAbs(postDifference);
        
        
        //self.expected_cable_activities += update_rate_post * (self.pre * self.post.T - self.expected_cable_activities)                
        expectedCableDelta = matrixVector(ecDeltaFactor, updateRatePost, expectedCableDelta);
        addEquals(expectedCableActivities, expectedCableDelta);
        

        
        //self.post_uncertainty += (post_difference - self.post_uncertainty) * update_rate_post
        
        puDelta = ensureCopy(puDelta, postDifference);
        subEquals(puDelta, postUncertainty);        
        matrixVector(puDelta, updateRatePost);
        addEquals(postUncertainty, puDelta);
                
        //# Reaction is the expected post, turned into a deliberation_vote
        //self.reaction = tools.weighted_average(self.expected_cable_activities, self.pre)        
            //daisychain stepup reaction param (8, 8) (8, 1)
        getWeightedAverage(expectedCableActivities, pre, reaction);
        

        //# Surprise is the difference between the expected post and the actual one
        //self.surprise = tools.weighted_average(
        //                  np.abs(self.post.T - self.expected_cable_activities),
        //                  self.pre / (self.post_uncertainty + tools.EPSILON))
        DenseMatrix64F sa = broadcastCols(postT, expectedCableActivities.getNumRows());
        subEquals(sa, expectedCableActivities);
        matrixAbs(sa);

        sb = ensureCopy(sb, postUncertainty);
        add(sb, EPSILON);

        matrixDivBy(sb, broadcastRows(pre, sb.getNumRows()));        
        
        getWeightedAverage(sa, sb, surprise);
        
        /*if (elementSum(surprise) > 0) {
            System.out.println(postT);
            System.out.println(broadcastCols(postT, expectedCableActivities.getNumRows()));
            System.out.println(expectedCableActivities);
            System.out.print(elementSum(surprise) + " = ");
            printArray(surprise.getData());
            System.out.println("---");
        }*/
        
        //# Reshape chain activities into a single column
        //return chain_activities.ravel()[:,np.newaxis]
        //check the order: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ravel.html?highlight=ravel        
        return DenseMatrix64F.wrap(chainActivities.elements, 1, chainActivities.getData() );
    }
    
    public DenseMatrix64F stepDown(DenseMatrix64F chainGoals) {
        //""" Propogate goals down through the transition model """
        
        
        //# Reshape chain_goals back into a square array 
        //chain_goals = np.reshape(chain_goals, (self.post.size, -1))
        chainGoals = DenseMatrix64F.wrap(
                post.elements, chainGoals.elements / post.elements,
                chainGoals.getData()
        );
                        
        //# Weight chain goals by the current cable activities   
        //upstream_goals = tools.bounded_sum(self.post * chain_goals.T)
                
        DenseMatrix64F upstreamGoals = matrixVector(transpose(chainGoals, null), post);
                
        upstreamGoals = boundedRowSum(upstreamGoals);
        
        //cable_goals = tools.bounded_sum([upstream_goals, self.reaction])
        DenseMatrix64F cableGoals = 
                DenseMatrix64F.wrap(upstreamGoals.getNumRows(), 1,
                        boundedSum(0, upstreamGoals.getData(), 
                                reaction.getData()));
        
        //return cable_goals[:self.num_cables]        
        if (numCables == 0) {
            if (nullGoals == null)
                nullGoals = new DenseMatrix64F(cableGoals.getNumRows(), 0);
            else
                if (nullGoals.numRows != cableGoals.numRows)
                    nullGoals.numRows = cableGoals.numRows;
            return nullGoals;
        }
        else {
//            if (cableGoals.getData()!=null)
//                if (cableGoals.getData().length>0) {
//                    //reaction is causing the '1 . 0 0 0 0 '
//                    System.out.println(upstreamGoals.getData()[0] + " " + reaction.getData()[0]);
//                    System.out.println(expectedCableActivities);
//                    System.out.println(pre);
//
//                    new Exception().printStackTrace();
//                }
            
            return extract(cableGoals, 0, numCables, 0, 1);
        }
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

    public DenseMatrix64F getReaction() {
        return reaction;
    }

    public DenseMatrix64F getCount() {
        return count;
    }

    public DenseMatrix64F getPre() {
        return pre;
    }

    public DenseMatrix64F getPreCount() {
        return preCount;
    }

    public DenseMatrix64F getPost() {
        return post;
    }

    public DenseMatrix64F getPostUncertainty() {
        return postUncertainty;
    }

    public DenseMatrix64F getExpectedCableActivities() {
        return expectedCableActivities;
    }


    
}
