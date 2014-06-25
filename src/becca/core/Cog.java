/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.core;

import org.ejml.data.DenseMatrix64F;

import static becca.core.Util.*;


/**
    The basic units of which blocks are composed

    Cogs are named for their similarity to clockwork cogwheels.
    They are simple and do the same task over and over, but by
    virtue of how they are connected to their fellows, they 
    collectively bring about interesting behavior.  

    Input channels are similar to cables in that they carry activity 
    signals that vary over time.
    Each cog contains two important parts, a daisychain and a ziptie.
    The daisychain is an object that builds cables into short sequences,
    and the ziptie is an object that takes the resulting chains
    and performs clustering on them, creating bundles.
    During upward processing, cable activities are used to train
    the daisychain and ziptie, making new bundles, maturing existing 
    bundles, and calculating the activity in bundle. 
    During downward processing, 
    the daisychain and ziptie use the bundle activity goals from 
    the next level higher to create goals for the cables. 
*/
public class Cog {
    public final int maxCables;
    public final int maxBundles;
    public final int maxChainsPerBundle;
    public final DaisyChain daisychain;
    public final ZipTie ziptie;
    private DenseMatrix64F bundleActivities;
    private DenseMatrix64F preCogCableActivities;
    private boolean preEnoughCable;
    private DenseMatrix64F cableGoals;
    private DenseMatrix64F preStepDownGoals;
    private final DenseMatrix64F emptyActivity;

    public Cog(int maxCables, int maxBundles, int maxChainsPerBundle, int level) {
        
        this.maxCables = maxCables;
        this.maxBundles = maxBundles;
        
        if (maxChainsPerBundle == 0) 
            maxChainsPerBundle = (int)(Math.pow(maxCables,2) / ((double)maxBundles) );
        
        this.maxChainsPerBundle = maxChainsPerBundle;
        
        this.daisychain = new DaisyChain(maxCables);        
        this.emptyActivity = new DenseMatrix64F(maxBundles, 1);
        
        if (maxBundles > 0)
            this.ziptie = new ZipTie(false, (int)Math.pow(maxCables, 2), maxBundles, maxChainsPerBundle);
        else
            this.ziptie = null;
        
    }

    //""" cable_activities percolate upward through daisychain and ziptie """
    public DenseMatrix64F stepUp(DenseMatrix64F cableActivities, boolean enoughCables) {                
        if (cableActivities == null) {
            cableActivities = preCogCableActivities;
            enoughCables = preEnoughCable;
        }

        /*
        # TODO: fix this so that cogs can gracefully handle more cables 
        # or else never be assigned them in the first place
        if cable_activities.size > self.max_cables:
            cable_activities = cable_activities[:self.max_cables, :]
            print '-----  Number of max cables exceeded in', self.name, \
                    '  -----'
        */
        if (cableActivities.getNumRows() > maxCables) {
            cableActivities.reshape(maxCables,1);
            System.err.println("Cog: Number of max cables exceeded in " + this);
        }
        
        final DenseMatrix64F chainActivities = daisychain.stepUp(cableActivities);
        
        if (enoughCables) {
            bundleActivities = ziptie.stepUp(chainActivities);//.copy();
        }
        else {
            bundleActivities = emptyActivity;
        }
        
        if (cableActivities.getNumRows() < maxBundles)
            bundleActivities = pad(bundleActivities, maxBundles, 1, 0.0);
                
        return bundleActivities;
    }

    public DenseMatrix64F getBundleActivities() {
        return bundleActivities;
    }

    public DenseMatrix64F getCableGoals() {
        return cableGoals;
    }
        
    
    
    //""" bundle_goals percolate downward """
    public DenseMatrix64F stepDown(DenseMatrix64F bundleGoals) {
        if (bundleGoals == null)
            bundleGoals = preStepDownGoals;
        
        final DenseMatrix64F chainGoals = ziptie.stepDown(bundleGoals);
        cableGoals = daisychain.stepDown(chainGoals);
        
        return cableGoals;
    }
    
    //""" How many bundles have been created in this cog? """
    public int getNumBundles() {            
        return ziptie.getNumBundles();
    }
    
    //""" How full is the set of cables for this cog? """
    public double getFractionFilled() {    
        return ((double)daisychain.getNumCables()) / ((double)maxCables);
    }
    
    /*           
    def get_index_projection(self, bundle_index):
        """ Project a bundle down through the ziptie and daisychain """
        chain_projection = self.ziptie.get_index_projection(bundle_index)
        cable_projection = self.daisychain.get_index_projection(
                chain_projection)
        return cable_projection
         
    */
    
    //""" Show the internal state of the daisychain and ziptie """
    public String toString() {
        String x = this.daisychain.toString();
        if (this.maxBundles > 0)
            x += " -> " + this.ziptie;
        return x;
    }

    public DenseMatrix64F getSurprise() {
        return daisychain.getSurprise();
    }

    /** used to preload parameters in case it is parallel executed */
    void preStepUp(DenseMatrix64F cogCableActivities, boolean enoughCables) {
        this.preCogCableActivities = cogCableActivities;
        this.preEnoughCable = enoughCables;
    }
    void preStepDown(DenseMatrix64F goals) {
        this.preStepDownGoals = goals;
    }

    
}
