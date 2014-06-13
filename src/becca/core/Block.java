/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.core;

import java.util.ArrayList;
import org.ejml.data.BlockMatrix64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static becca.core.Util.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;

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
    
    public final ArrayList<Cog> cogs;
    
    private DenseMatrix64F cableActivities;
    private final double fillFractionThreshold;
    private final double activityDecayRate;
    private final double rangeDecayRate;
    private BlockMatrix64F maxVals;
    private BlockMatrix64F minVals;
    private DenseMatrix64F surprise;
    private DenseMatrix64F bundleActivities;
    private DenseMatrix64F hubCableGoals;

    public Block(int minCables) {    
        this(minCables, 0);
    }
    
    //""" Initialize the level, defining the dimensions of its cogs """
    public Block(int minCables, int level) {
        this.maxCables = (int)Math.round( Math.pow(2, 1 + Math.ceil(Util.log(minCables, 2))));
        this.level = level;
        
        this.maxCablesPerCog = 8;
        this.maxBundlesPerCog = 4;
        this.maxCogs = maxCables / maxBundlesPerCog;
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
        
        this.cableActivities = new DenseMatrix64F(this.maxCables,1); //np.zeros((self.max_cables, 1))
        this.hubCableGoals = new DenseMatrix64F(this.maxCables, 1); //np.zeros((self.max_cables, 1))
        
        this.fillFractionThreshold = 0.7;
        
        this.activityDecayRate = 1.0;       //real, 0 < x < 1
        
        //# Constants for adaptively rescaling the cable activities
        this.rangeDecayRate = Math.pow(10, -5);        
        this.maxVals = new BlockMatrix64F(this.maxCables, 1); // np.zeros((self.max_cables, 1)) 
        this.minVals = new BlockMatrix64F(this.maxCables, 1); // np.zeros((self.max_cables, 1))

    }

    
    public DenseMatrix64F stepUp(DenseMatrix64F newCableActivities) {
                
//        """ Find bundle_activities that result from new_cable_activities """
//        # Condition the cable activities to fall between 0 and 1
        
        if (newCableActivities.getNumRows() < maxCables) {           
              newCableActivities = Util.pad(newCableActivities, maxCables, 1, 0);
        }
        
        double[] ncaD = newCableActivities.getData();
        double[] minD = minVals.getData();
        double[] maxD = maxVals.getData();
        for (int i = 0; i < maxCables; i++) {
            minD[i] = Math.min(minD[i], ncaD[i]);
            maxD[i] = Math.max(maxD[i], ncaD[i]);
        }
        
        final BlockMatrix64F spread = maxVals.copy();
        CommonOps.subEquals(spread, minVals);
        
//        new_cable_activities = ((new_cable_activities - self.min_vals) / 
//                   (self.max_vals - self.min_vals + tools.EPSILON))
        for (int i = 0; i < maxCables; i++) {
            final double n = ncaD[i];
            final double min = minD[i];
            final double max = maxD[i];
            
            final double v = (n-min)/(max-min+Util.EPSILON);
            
            ncaD[i] = v;
        
        }
        
        BlockMatrix64F spreadDecay = spread.copy();
        CommonOps.scale(rangeDecayRate, spreadDecay);
        CommonOps.addEquals(minVals, spreadDecay);
        CommonOps.subEquals(maxVals, spreadDecay);

        DenseMatrix64F decayedActivities = cableActivities.copy();
        CommonOps.scale(1.0 - activityDecayRate, decayedActivities);
                
        assert(newCableActivities.getNumRows() == decayedActivities.getNumRows() );
        cableActivities = DenseMatrix64F.wrap(newCableActivities.getNumRows(), 1,
                            boundedSum(0, newCableActivities.getData(), decayedActivities.getData() ));        
        
//        # Update cable_activities, incorporating sensing dynamics
//        self.cable_activities = tools.bounded_sum([
//                new_cable_activities, 
//                self.cable_activities * (1. - self.ACTIVITY_DECAY_RATE)])
//
//        # Update the map from self.cable_activities to cogs
//        self.ziptie.step_up(self.cable_activities)
//        
        ziptie.stepUp(cableActivities);
        
        //# Process the upward pass of each of the cogs in the block        

        ArrayList<double[]> bundleActivitiez = new ArrayList(cogs.size());
        int numBundleActivitiez = 0;
        
        int cogIndex = 0;
        for (Cog c : cogs) {
            /* # Pick out the cog's cable_activities, process them, 
               # and assign the results to block's bundle_activities*/
            
            //cog_cable_activities = self.cable_activities[self.ziptie.get_index_projection(cog_index).ravel().astype(bool)]
            DenseMatrix64F cogCableActivities = Util.extractBooleanized(cableActivities, ziptie.getIndexProjection(cogIndex));
            
            /*# Cogs are only allowed to start forming bundles once 
              # the number of cables exceeds the fill_fraction_threshold*/
            boolean enoughCables = ziptie.getCableFractionInBundle(cogIndex) > fillFractionThreshold;
               
            DenseMatrix64F cogBundleActivities = c.stepUp(cogCableActivities, enoughCables);
                    
            //self.bundle_activities = np.concatenate((self.bundle_activities,cog_bundle_activities))
            double[] cbaData = cogBundleActivities.getData();
            bundleActivitiez.add(cbaData);
            numBundleActivitiez += cbaData.length;
            
            cogIndex++;
        }
        this.bundleActivities = new DenseMatrix64F(numBundleActivitiez, 1);
        int p = 0;
        for (int i = 0; i < bundleActivitiez.size(); i++) {
            double[] c = bundleActivitiez.get(i);
            System.arraycopy(c, 0, bundleActivities.getData(), p, c.length);
            p+=c.length;
        }
        //concat finished
        
        
        //# Goal fulfillment and decay
        /*self.hub_cable_goals -= self.cable_activities
        self.hub_cable_goals *= self.ACTIVITY_DECAY_RATE
        self.hub_cable_goals = np.maximum(self.hub_cable_goals, 0.) */
        subEquals(hubCableGoals, cableActivities);
        scale(activityDecayRate, hubCableGoals);
        matrixMaximum(hubCableGoals, 0);
        
        
        return bundleActivities;
    }

    public DenseMatrix64F stepDown(DenseMatrix64F bundleGoals) {

                
        //""" Find cable_activity_goals, given a set of bundle_goals """
        
        //bundle_goals = tools.pad(bundle_goals, (self.max_bundles, 1))        
        bundleGoals = pad(bundleGoals, maxBundles, 1, 0.0);
        
        DenseMatrix64F cableGoals = new DenseMatrix64F(maxCables, 1);
        
        //self.surprise = np.zeros((self.max_cables, 1))
        surprise = new DenseMatrix64F(maxCables, 1);


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
            DenseMatrix64F cableGoalsByCog = c.stepDown(cogBundleGoals);
                        
            //cog_cable_indices = self.ziptie.get_index_projection(cog_index).ravel().astype(bool)
            DenseMatrix64F cogCableIndices = matrixBooleanize(ziptie.getIndexProjection(cogIndex));
            assert(cogCableIndices.getNumRows() == 1);
                       
            //paddingReaction = pad(c.getReaction(), cogCableIndices.getNumRows(), 0, 0.0);
            final DenseMatrix64F cs = c.getSurprise();
            
            final double[] ccid = cogCableIndices.getData();
            
            for (int i = 0; i < ccid.length; i++) {
                if (ccid[i]>0) {
                    
                    for (int j = 0; j < cableGoals.getNumCols(); j++) {
                        
                        //System.out.println(i + " " + j + " " + m(cableGoals) + " " + m(cableGoalsByCog)+ " " + m(cogBundleGoals));
                        
                        //TODO: DECIDE IF THIS IS CORRECT
                        //cable_goals[cog_cable_indices] = np.maximum(cable_goals_by_cog, cable_goals[cog_cable_indices])            
                        if (cableGoalsByCog.getNumRows() > i)
                            cableGoals.set(i, j, 
                                    Math.max(cableGoals.get(i, j), cableGoalsByCog.get(i, 0)));
                    }
                    
                        //#self.reaction[cog_cable_indices] = np.maximum(
                        //#        tools.pad(cog.reaction, (cog_cable_indices[0].size, 0)),
                        //#        self.reaction[cog_cable_indices]) 


                        //TODO: DECIDE IF THIS IS CORRECT
                        //self.surprise[cog_cable_indices] = np.maximum(cog.surprise, self.surprise[cog_cable_indices])
                        //System.out.println(j + " " + i + " " + m(cs) + " " + m(surprise));                                               
                    for (int j = 0; j < surprise.getNumCols(); j++) {
                    
                        if (cs.getNumRows() > i)
                            surprise.set(i, j, 
                                    Math.max(surprise.get(i, j), cs.get(i, 0)));

                    }
                
                
                }
                
            }
            
        }
        
        //hubCableGoals.setData( boundedSum(0, hubCableGoals.getData(), cableGoals.getData() ) );
                        
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
