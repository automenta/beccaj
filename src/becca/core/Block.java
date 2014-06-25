/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.core;

import java.util.ArrayList;
import org.ejml.data.DenseMatrix64F;

import static becca.core.Util.*;
import static java.lang.Double.NaN;

/**
    Blocks are the building block of which the agent is composed

    Blocks are arranged hierarchically within the agent. 
    The agent begins with only one block, and creates additional
    blocks in a tower arrangement as lower ones mature. 
    
    The block's input channels (cables) are organized 
    into clusters (bundles)
    whose activities are passed up to the next block in the hierarchy.  
    Each block performs the same two functions, 
    1) a step_up 
    where cable activities are converted to bundle activities and
    passed up the tower and 
    2) a step_down where bundle activity goals are passed back down
    and converted into cable activity goals. 
    Internally, a block contains a number of cogs that work in parallel
    to convert cable activities into bundle activities and back again.
    
    */
//TODO get changes from: https://github.com/brohrer/becca/commit/4b0c6270b3974478213831652fb3895c58734daf
public class Block  {
    public final int maxCables;
    private final int maxBundlesPerCog;
    private final int maxCablesPerCog;
    private final int maxCogs;
    public final int maxBundles;
    public final int level;
    public final ZipTie ziptie;
    
    public final ArrayList<Cog> cogs;
    
    private DenseMatrix64F cableActivities;
    private final double fillFractionThreshold;
    private final double activityDecayRate;
    private final double rangeDecayRate;
    private DenseMatrix64F maxVals;
    private DenseMatrix64F minVals;
    private final DenseMatrix64F surprise;
    private DenseMatrix64F bundleActivities;
    private DenseMatrix64F hubCableGoals;
    private boolean parallelCogs = BeccaParams.cogParallel;
    int stepCounter = 0;
    private final DenseMatrix64F rawCableActivities;
    private final DenseMatrix64F prevCableActivities;
    private final int stepMultiplier;
    private final int stepCounterOffset;
    private final DenseMatrix64F clusterTrainingActivities;
    
    public Block(int minCables) {    
        this(minCables, 0);
    }
    

    public Block(int minCables, int level) {
        this(minCables, (int)Math.round( Math.pow(2, 1 + Math.ceil(Util.log(minCables, 2)))), level);
    }

   //""" Initialize the level, defining the dimensions of its cogs """
    public Block(int minCables, int maxCables, int level) {        
        this.level = level;
        this.maxCables = maxCables;        
        
        this.maxCablesPerCog = BeccaParams.blockMaxCablesPerCog;
        this.maxBundlesPerCog = BeccaParams.blockMaxBundlesPerCog;
        this.maxCogs = Math.max(1, maxCables / maxBundlesPerCog);
        this.maxBundles = this.maxCogs * this.maxBundlesPerCog;
        
        this.ziptie = new ZipTie(true, this.maxCables, this.maxCogs, this.maxCablesPerCog, -2);
            //ziptie_name = ''.join(('ziptie_', self.name))
        
        this.cogs = new ArrayList<Cog>();
        for (int i = 0; i < this.maxCogs; i++) {
            final Cog c = new Cog(this.maxCablesPerCog,
                                 this.maxBundlesPerCog,
                                  this.maxCablesPerCog,
                                   this.level);
            cogs.add(c);
        }
        
        this.cableActivities = new DenseMatrix64F(this.maxCables,1); //np.zeros((self.max_cables, 1))
        this.rawCableActivities = new DenseMatrix64F(this.maxCables,1); //np.zeros((self.max_cables, 1))
        this.prevCableActivities = new DenseMatrix64F(this.maxCables,1); //np.zeros((self.max_cables, 1))
        this.hubCableGoals = new DenseMatrix64F(this.maxCables, 1); //np.zeros((self.max_cables, 1))
        
        this.fillFractionThreshold = BeccaParams.blockFillFractionThreshold;
        this.activityDecayRate = BeccaParams.blockActivityDecayRate;
        this.rangeDecayRate = BeccaParams.blockRangeDecayRate; //# Constants for adaptively rescaling the cable activities        
        
        this.surprise = new DenseMatrix64F(this.maxCables, 1);
        
        this.maxVals = new DenseMatrix64F(this.maxCables, 1); // np.zeros((self.max_cables, 1)) 
        this.minVals = new DenseMatrix64F(this.maxCables, 1); // np.zeros((self.max_cables, 1))
        this.bundleActivities = new DenseMatrix64F(this.maxCables, 1); // np.zeros((self.max_cables, 1))
        this.clusterTrainingActivities = new DenseMatrix64F(this.maxCables, 1);  // n        
        this.stepMultiplier = (int)Math.pow(2, level);
        this.stepCounterOffset = 1000000;
        /*
        self.bundle_activities = np.zeros((self.max_bundles, 1))
         self.raw_cable_activities = np.zeros((self.max_cables, 1))
         self.previous_cable_activities = np.zeros((self.max_cables, 1))        
        */
        
    }

    
    public DenseMatrix64F stepUp(DenseMatrix64F newCableActivities) {
                
//        """ Find bundle_activities that result from new_cable_activities """
//        # Condition the cable activities to fall between 0 and 1
        
        if (newCableActivities.getNumRows() < maxCables) {           
              newCableActivities = pad(newCableActivities, maxCables, 1, 0);
        }
        
        
        /*self.raw_cable_activities *= 1. - (self.ACTIVITY_DECAY_RATE /
                                           float(self.step_multiplier))
        self.raw_cable_activities += new_cable_activities
        self.step_counter += 1
        if self.step_counter < self.step_multiplier:
            return self.bundle_activities*/
        scale( 1.0 - (activityDecayRate / ((double)stepMultiplier)),  rawCableActivities);
        addEquals(rawCableActivities, newCableActivities);
        if (this.stepCounter + stepCounterOffset < stepMultiplier)
            return bundleActivities;
                    
        double[] rcad = rawCableActivities.getData();
        double[] minD = minVals.getData();
        double[] maxD = maxVals.getData();
        for (int i = 0; i < maxCables; i++) {
            minD[i] = Math.min(minD[i], rcad[i]);
            maxD[i] = Math.max(maxD[i], rcad[i]);
        }
        
        final DenseMatrix64F spread = maxVals.copy();
        subEquals(spread, minVals);
        
//        self.cable_activities = ((self.raw_cable_activities - self.min_vals) / 
//                   (self.max_vals - self.min_vals + tools.EPSILON))
        double[] caD = this.cableActivities.getData();
        for (int i = 0; i < cableActivities.elements; i++) {            
            final double min = minD[i];
            final double max = maxD[i];
            
            caD[i] = (rcad[i]-min)/(max-min+EPSILON);        
        }
        
        DenseMatrix64F spreadDecay = spread.copy();
        scale(rangeDecayRate, spreadDecay);
        addEquals(minVals, spreadDecay);
        subEquals(maxVals, spreadDecay);

        /*cluster_training_activities = np.maximum(
                self.previous_cable_activities, self.cable_activities)        
        */
        double[] ctad = clusterTrainingActivities.getData();
        for (int i = 0; i < prevCableActivities.elements; i++) {
            ctad[i] = Math.max(prevCableActivities.getData()[i], cableActivities.getData()[i]);
        }
        
        //self.previous_cable_activities = self.cable_activities.copy() 
        prevCableActivities.set(cableActivities);
        
        

        
//        # Update cable_activities, incorporating sensing dynamics
//        self.cable_activities = tools.bounded_sum([
//                new_cable_activities, 
//                self.cable_activities * (1. - self.ACTIVITY_DECAY_RATE)])
//
//        # Update the map from self.cable_activities to cogs
//        self.ziptie.step_up(self.cable_activities)
//        
        
        
        ziptie.stepUp(clusterTrainingActivities);
        
        //# Process the upward pass of each of the cogs in the block        

        
        
        int cogIndex = 0;
        for (final Cog c : cogs) {
            /* # Pick out the cog's cable_activities, process them, 
               # and assign the results to block's bundle_activities*/
            
            //cog_cable_activities = self.cable_activities[self.ziptie.get_index_projection(cog_index).astype(bool)]
            DenseMatrix64F cogCableActivities = Util.extractBooleanized(cableActivities, ziptie.getIndexProjection(cogIndex));
            
            

            /*# Cogs are only allowed to start forming bundles once 
              # the number of cables exceeds the fill_fraction_threshold*/
            boolean enoughCables = ziptie.getCableFractionInBundle(cogIndex) > fillFractionThreshold;
            
            //System.out.print(enoughCables + " ");

            c.preStepUp(cogCableActivities, enoughCables);
            cogIndex++;
        }
        

        if (parallelCogs) {
            cogs.parallelStream().forEach(c -> c.stepUp(null, false));
        }
        else {
            for (final Cog c : cogs) {
                c.stepUp(null, false);
            }
        }        
        
        fill(this.bundleActivities, 0);
        
        int numBundleActivitiez = 0;
        cogIndex = 0;        
        for (final Cog c : cogs) {                                                   
            //self.bundle_activities = np.concatenate((self.bundle_activities,cog_bundle_activities))
            DenseMatrix64F cogBundleActivities = c.getBundleActivities();
            double[] cbaData = cogBundleActivities.getData();
            System.arraycopy(cbaData, 0, this.bundleActivities.getData(), numBundleActivitiez, cogBundleActivities.elements);
            numBundleActivitiez += cogBundleActivities.elements;
            
            cogIndex++;
        }
        
        goalDecay();
        
        stepCounter++;
        
        return bundleActivities;
    }
    
    public void goalDecay() {
        //# Goal fulfillment and decay
        /*self.hub_cable_goals -= self.cable_activities
        self.hub_cable_goals *= self.ACTIVITY_DECAY_RATE
        self.hub_cable_goals = np.maximum(self.hub_cable_goals, 0.) */
        subEquals(hubCableGoals, cableActivities);
        scale(activityDecayRate, hubCableGoals);
        matrixMaximum(hubCableGoals, 0);        
    }

    public DenseMatrix64F stepDown(DenseMatrix64F bundleGoals) {

                
        //""" Find cable_activity_goals, given a set of bundle_goals """
        
        //bundle_goals = tools.pad(bundle_goals, (self.max_bundles, 1))        
        
        bundleGoals = pad(bundleGoals, maxBundles, 1, 0.0);
        
        
        



        //# Process the downward pass of each of the cogs in the level
        for (int cogIndex = 0; cogIndex < cogs.size(); cogIndex++) {
            Cog c = cogs.get(cogIndex);

            /*
            #Gather the goal inputs for each cog
            cog_bundle_goals = bundle_goals[cog_index * self.max_bundles_per_cog:cog_index + 1 * self.max_bundles_per_cog,:]
            */
            DenseMatrix64F cogBundleGoals = extract(bundleGoals, 
                    cogIndex * maxBundlesPerCog, (cogIndex+1) * maxBundlesPerCog,
                    0, bundleGoals.getNumCols()
                    );
                        
            //# Update the downward outputs for the level 
            //cable_goals_by_cog = cog.step_down(cog_bundle_goals)
            c.preStepDown(cogBundleGoals);
        }
        
        
        if (parallelCogs) {
            cogs.parallelStream().forEach(c -> c.stepDown(null));
        }
        else {
            for (final Cog c : cogs) {                
                c.stepDown(null);
            }
        }        
        
        
        
        final DenseMatrix64F cableGoals = new DenseMatrix64F(maxCables, 1);

        //self.surprise = np.zeros((self.max_cables, 1))
        fill(surprise, 0);            

        for (int cogIndex = 0; cogIndex < cogs.size(); cogIndex++) {
            final Cog c = cogs.get(cogIndex);         
            final DenseMatrix64F cableGoalsByCog = c.getCableGoals();
            final double[] cableGoalsByCogD = cableGoalsByCog.getData();
            
            //cog_cable_indices = self.ziptie.get_index_projection(cog_index).astype(bool)
            DenseMatrix64F cogCableIndices = matrixBooleanize(ziptie.getIndexProjection(cogIndex));

            
            assert(cogCableIndices.getNumRows() == 1);
            
            //cog_cable_activities = self.cable_activities[self.ziptie.get_index_projection(cog_index).astype(bool)]            
            //DenseMatrix64F cogCableActivities = Util.extractBooleanized(cableActivities, ziptie.getIndexProjection(cogIndex));            
                       
            //paddingReaction = pad(c.getReaction(), cogCableIndices.getNumRows(), 0, 0.0);
            final DenseMatrix64F cogSurprise = c.getSurprise(); //transpose(c.getSurprise(), null);
            

        
            final double[] csd = cogSurprise.getData();
            final double[] ccid = cogCableIndices.getData();
            final double[] cgd = cableGoals.getData();
            final double[] sd = this.surprise.getData();
            
            //printArray(cogCableIndices.getData());
            
            assert(cogCableIndices.elements == cableGoals.elements);
            assert(cogCableIndices.elements == surprise.elements);
            
            assert(cableGoalsByCog.numCols < 2);
            for (int j = 0; j < cableGoalsByCog.elements; j++) {       
                if (ccid[j] > 0) {
                    cgd[j] = Math.max(cgd[j], cableGoalsByCogD[j]);
                    sd[j] = Math.max(sd[j], csd[j]);
                }
            }     
            
        }                
        
        //System.out.println(transpose(cableGoals,null));
        //System.out.println(transpose(hubCableGoals,null));

        if (BeccaParams.BlockGoalBoundedSum)
            hubCableGoals.setData( boundedSum(0, hubCableGoals.getData(), cableGoals.getData() ) );
        
        
        //System.out.println(transpose(hubCableGoals,null));
        //System.out.println();
        
        
        
       //test if action contains NaN or Inf        
        final double[] action = hubCableGoals.getData();
        boolean invalidAction = false;
        for (int ii = 0; ii < action.length; ii++) {
            if ((!Double.isFinite(action[ii])) || (action[ii] == NaN)) {
                //System.out.println("WTF is " + action[ii]);
                invalidAction = true; 
                
                action[ii] = 0;
            }
        }
        if (invalidAction) {
            System.err.println("Invalid action");
            printArray(action);
            //System.out.println(this);
            //System.exit(1);            
        }

                
        return this.hubCableGoals;
    }


    
        
/*
    def get_index_projection(self, bundle_index):
        """ Represent one of the bundles in terms of its cables """
        # Find which cog it belongs to and which output it corresponds to
        cog_index = int(bundle_index / self.max_bundles_per_cog)
        cog_bundle_index = bundle_index - cog_index * self.max_bundles_per_cog
        # Find the projection to the cog's own cables
        cog_cable_indices = self.ziptie.get_index_projection(
                cog_index).ravel().astype(bool)
        num_cables_in_cog = np.sum(cog_cable_indices)
        cog_projection = self.cogs[cog_index].get_index_projection(
                cog_bundle_index)
        # Then re-sort them to the block's cables
        projection = np.zeros((self.max_cables, 2))
        projection[cog_cable_indices,:] = cog_projection[:num_cables_in_cog,:]
        return projection

    */
    
    public DenseMatrix64F getSurprise() {
        return surprise;
    }
    
    public int getBundlesCreated() {
        int total = 0;
        for (Cog c: cogs) {
            //# Check whether all cogs have created all their bundles
            total += c.getNumBundles();
        }
                        
        /*if np.random.random_sample() < 0.01:
            print total, 'bundles in', self.name, ', max of', self.max_bundles*/
        
        return total;
        
    }
    
    /*
    def visualize(self):
        """ Show what's going on inside the level """
        self.ziptie.visualize()
        #for cog in self.cogs:
        #    cog.visualize()
        return
    */    
    public String toString() {
        String x = super.toString() + " {" + ziptie + ",\n";
        for (Cog c : cogs) {
            x += ("    " + c + "\n");
        }
        x += "  }";
        return x;
    }

    public DenseMatrix64F getBundleActivities() {
        return bundleActivities;
    }

    public DenseMatrix64F getCableActivities() {
        return cableActivities;
    }

    public DenseMatrix64F getHubCableGoals() {
        return hubCableGoals;
    }

    
}
