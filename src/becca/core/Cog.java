/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.core;

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
    private final int maxCables;
    private final int maxBundles;
    private final int maxChainsPerBundle;
    private final DaisyChain daisychain;
    private ZipTie ziptie;


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
    
/*

    def step_up(self, cable_activities, enough_cables):
        """ cable_activities percolate upward through daisychain and ziptie """
        # TODO: fix this so that cogs can gracefully handle more cables 
        # or else never be assigned them in the first place
        if cable_activities.size > self.max_cables:
            cable_activities = cable_activities[:self.max_cables, :]
            print '-----  Number of max cables exceeded in', self.name, \
                    '  -----'
        chain_activities = self.daisychain.step_up(cable_activities)
        self.surprise = self.daisychain.get_surprise()
        if enough_cables is True:
            bundle_activities = self.ziptie.step_up(chain_activities)
        else:
            bundle_activities = np.zeros((0,1))
        bundle_activities = tools.pad(bundle_activities, (self.max_bundles, 0))
        return bundle_activities

    def step_down(self, bundle_goals):
        """ bundle_goals percolate downward """
        chain_goals = self.ziptie.step_down(bundle_goals) 
        cable_goals = self.daisychain.step_down(chain_goals)     
        return cable_goals

    def get_index_projection(self, bundle_index):
        """ Project a bundle down through the ziptie and daisychain """
        chain_projection = self.ziptie.get_index_projection(bundle_index)
        cable_projection = self.daisychain.get_index_projection(
                chain_projection)
        return cable_projection
         
    def fraction_filled(self):
        """ How full is the set of cables for this cog? """
        return float(self.daisychain.num_cables) / float(self.max_cables)

    def num_bundles(self):
        """ How many bundles have been created in this cog? """
        return self.ziptie.num_bundles
            
    def visualize(self):
        """ Show the internal state of the daisychain and ziptie """
        self.daisychain.visualize()
        if self.max_bundles > 0:
            self.ziptie.visualize()
        return    
    */    
    
    public String toString() {
        String x = this.daisychain.toString();
        if (this.maxBundles > 0)
            x += " -> " + this.ziptie;
        return x;
    }
}
