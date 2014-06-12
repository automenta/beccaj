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
    private DenseMatrix64F surprise;

    public Cog(int maxCables, int maxBundles, int maxChainsPerBundle, int level) {
        
        this.maxCables = maxCables;
        this.maxBundles = maxBundles;
        
        if (maxChainsPerBundle == 0) 
            maxChainsPerBundle = (int)(Math.pow(maxCables,2) / ((double)maxBundles) );
        
        this.maxChainsPerBundle = maxChainsPerBundle;
        
        this.daisychain = new DaisyChain(maxCables);        
        
        if (maxBundles > 0)
            this.ziptie = new ZipTie((int)Math.pow(maxCables, 2), maxBundles, maxChainsPerBundle);
        else
            this.ziptie = null;
        
    }

    //""" cable_activities percolate upward through daisychain and ziptie """
    public DenseMatrix64F stepUp(DenseMatrix64F activities, boolean enoughCables) {

        /*
        # TODO: fix this so that cogs can gracefully handle more cables 
        # or else never be assigned them in the first place
        if cable_activities.size > self.max_cables:
            cable_activities = cable_activities[:self.max_cables, :]
            print '-----  Number of max cables exceeded in', self.name, \
                    '  -----'
        */
        if (activities.getNumRows() > maxCables) {
            activities.reshape(maxCables,1);
            System.err.println("Cog: Number of max cables exceeded in " + this);
        }
        
        activities = daisychain.stepUp(activities);
        surprise = daisychain.getSurprise();
        
        if (enoughCables) {
            activities = ziptie.stepUp(activities);
        }
        else {
            activities = new DenseMatrix64F(0, 1);
        }
        
        activities = pad(activities, maxBundles, 1, 0.0);
        
        return activities;
    }
    
    //""" bundle_goals percolate downward """
    public DenseMatrix64F stepDown(DenseMatrix64F goals) {
        goals = ziptie.stepDown(goals);
        goals = daisychain.stepDown(goals);
        return goals;
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
        return surprise;
    }

    
}
