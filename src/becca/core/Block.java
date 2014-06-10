/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.core;

import java.util.ArrayList;
import org.encog.mathutil.matrices.Matrix;

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
public class Block {
    public final int maxCables;
    private final int maxBundlesPerCog;
    private final int maxCablesPerCog;
    private final int maxCogs;
    public final int maxBundles;
    public final int level;
    public final ZipTie ziptie;
    private final ArrayList<Cog> cogs;
    private final Matrix cableActivities;
    private Matrix hubCableGoals;
    private final double fillFractionThreshold;
    private final double activityDecayRate;
    private final double rangeDecayRate;
    private final Matrix maxVals;
    private final Matrix minVals;
    private Matrix surprise;
    private Matrix bundleActivities;

    Block(int minCables) {    
        this(minCables, 0);
    }
    
    //""" Initialize the level, defining the dimensions of its cogs """
    Block(int minCables, int level) {
        this.maxCables = (int)Math.round( Math.pow(2, Math.ceil(Util.log(minCables, 2))));
        this.level = level;
        
        this.maxCablesPerCog = 8;
        this.maxBundlesPerCog = 4;
        this.maxCogs = maxCablesPerCog / maxBundlesPerCog;
        this.maxBundles = this.maxCogs * this.maxBundlesPerCog;
        
        this.ziptie = new ZipTie(this.maxCables, this.maxCogs, this.maxCablesPerCog, -2);
            //ziptie_name = ''.join(('ziptie_', self.name))
        
        this.cogs = new ArrayList<Cog>();
        for (int i = 0; i < this.maxCogs; i++) {
            final Cog c = new Cog(this.maxCablesPerCog,
                                 this.maxBundlesPerCog,
                                  this.maxCablesPerCog,
                                   this.level);
            cogs.add(c);
        }
        
        this.cableActivities = new Matrix(this.maxCables,1); //np.zeros((self.max_cables, 1))
        this.hubCableGoals = new Matrix(this.maxCables, 1); //np.zeros((self.max_cables, 1))
        
        this.fillFractionThreshold = 0.7;
        
        this.activityDecayRate = 1.0;       //real, 0 < x < 1
        
        //# Constants for adaptively rescaling the cable activities
        this.rangeDecayRate = Math.pow(10, -5);        
        this.maxVals = new Matrix(this.maxCables, 1); // np.zeros((self.max_cables, 1)) 
        this.minVals = new Matrix(this.maxCables, 1); // np.zeros((self.max_cables, 1))

    }

    
    public Matrix stepUp(Matrix newCableActivities) {
        /*        
        """ Find bundle_activities that result from new_cable_activities """
        # Condition the cable activities to fall between 0 and 1
        if new_cable_activities.size < self.max_cables:
            new_cable_activities = tools.pad(new_cable_activities, 
                                             (self.max_cables, 1))
        self.min_vals = np.minimum(new_cable_activities, self.min_vals)
        self.max_vals = np.maximum(new_cable_activities, self.max_vals)
        spread = self.max_vals - self.min_vals
        new_cable_activities = ((new_cable_activities - self.min_vals) / 
                   (self.max_vals - self.min_vals + tools.EPSILON))
        self.min_vals += spread * self.RANGE_DECAY_RATE
        self.max_vals -= spread * self.RANGE_DECAY_RATE
        # Update cable_activities, incorporating sensing dynamics
        self.cable_activities = tools.bounded_sum([
                new_cable_activities, 
                self.cable_activities * (1. - self.ACTIVITY_DECAY_RATE)])

        # Update the map from self.cable_activities to cogs
        self.ziptie.step_up(self.cable_activities)
        */
        
        //# Process the upward pass of each of the cogs in the block        
        this.bundleActivities = new Matrix(0, 1);   //self.bundle_activities = np.zeros((0, 1))
        
        /*
        for cog_index in range(len(self.cogs)):
            # Pick out the cog's cable_activities, process them, 
            # and assign the results to block's bundle_activities
            cog_cable_activities = self.cable_activities[
                    self.ziptie.get_index_projection(
                    cog_index).ravel().astype(bool)]
            # Cogs are only allowed to start forming bundles once 
            # the number of cables exceeds the fill_fraction_threshold
            enough_cables = (self.ziptie.cable_fraction_in_bundle(cog_index)
                             > self.fill_fraction_threshold)
            cog_bundle_activities = self.cogs[cog_index].step_up(
                    cog_cable_activities, enough_cables)
            self.bundle_activities = np.concatenate((self.bundle_activities, 
                                                     cog_bundle_activities))
        # Goal fulfillment and decay
        self.hub_cable_goals -= self.cable_activities
        self.hub_cable_goals *= self.ACTIVITY_DECAY_RATE
        self.hub_cable_goals = np.maximum(self.hub_cable_goals, 0.)
        return self.bundle_activities
    */
        return bundleActivities;
    }

    public Matrix stepDown(Matrix bundleGoals) {
        /*
        """ Find cable_activity_goals, given a set of bundle_goals """
        bundle_goals = tools.pad(bundle_goals, (self.max_bundles, 1))
        */
        
        Matrix cableGoals = new Matrix(maxCables, 1);
        
        /*
        self.surprise = np.zeros((self.max_cables, 1))
        # Process the downward pass of each of the cogs in the level
        cog_index = 0
        for cog in self.cogs:
            # Gather the goal inputs for each cog
            cog_bundle_goals = bundle_goals[
                    cog_index * self.max_bundles_per_cog:
                    cog_index + 1 * self.max_bundles_per_cog,:]
            # Update the downward outputs for the level 
            cable_goals_by_cog = cog.step_down(cog_bundle_goals)
            cog_cable_indices = self.ziptie.get_index_projection(
                    cog_index).ravel().astype(bool)
            cable_goals[cog_cable_indices] = np.maximum(
                    cable_goals_by_cog, cable_goals[cog_cable_indices]) 
            #self.reaction[cog_cable_indices] = np.maximum(
            #        tools.pad(cog.reaction, (cog_cable_indices[0].size, 0)),
            #        self.reaction[cog_cable_indices]) 
            self.surprise[cog_cable_indices] = np.maximum(
                    cog.surprise, self.surprise[cog_cable_indices]) 
            cog_index += 1
        */
        hubCableGoals = Util.boundedSum(0, hubCableGoals, cableGoals);
        
        return hubCableGoals;
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
    
    public Matrix getSurprise() {
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

    
}
