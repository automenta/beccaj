/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.core;

import org.ejml.data.BlockMatrix64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.data.Matrix64F;
import org.ejml.ops.CommonOps;

/**
    An incremental unsupervised learning algorithm

    Input channels are clustered together into mutually co-active sets.
    A helpful metaphor is bundling cables together with zip ties.
    Cables that carry related signals are commonly co-active, 
    and are grouped together. Cables that are co-active with existing
    bundles can be added to those bundles. A single cable may be ziptied
    into several different bundles. Co-activity is estimated 
    incrementally, that is, the algorithm updates the estimate after 
    each new set of signals is received.     
 */
public class ZipTie {
    private final int maxCables;
    private final int maxBundles;
    private final int maxCablesPerBundle;
    private final int numBundles;
    private final double AGGLOMERATION_ENERGY_RATE;
    private final double NUCLEATION_ENERGY_RATE;
    private final double ENERGY_DECAY_RATE;
    private final double JOINING_THRESHOLD;
    private final double NUCLEATION_THRESHOLD;
    private final double MEAN_EXPONENT;
    private final double ACTIVATION_WEIGHTING_EXPONENT;
    private final boolean bundlesFull;
    
    private final Matrix64F bundleActivities;
    private final int[] mapSize;
    private final Matrix64F bundleMap;
    private final Matrix64F agglomerationEnergy;
    private final Matrix64F nucleationEnergy;

    public ZipTie(int maxCables, int maxBundles, int maxCablesPerBundle) {
        this(maxCables, maxBundles, maxCablesPerBundle, -4);
    }

    public ZipTie(int maxCables, int maxBundles, int maxCablesPerBundle, double meanExponent) {
        this(maxCables, maxBundles, maxCablesPerBundle, meanExponent, 0.05, 1.0);
    }
    
    public ZipTie(int maxCables, int maxBundles, int maxCablesPerBundle, double meanExponent, double joiningThreshold, double speedup) {
        
        this.maxCables = maxCables;
        this.maxBundles = maxBundles;
        
        if (maxCablesPerBundle == 0)
            maxCablesPerBundle = (int)( ((double)maxCables) / ((double)maxBundles) );
        
        this.maxCablesPerBundle = maxCablesPerBundle;
        
        this.numBundles = 0;
                
        //User-defined constants        
        this.AGGLOMERATION_ENERGY_RATE = Math.pow(10,-2) * speedup;
        this.NUCLEATION_ENERGY_RATE = Math.pow(10,-4) * speedup;
        this.ENERGY_DECAY_RATE = Math.pow(10, -2);

                
        /*
        # Coactivity value which, if it's every exceeded, causes a 
        # cable to be added to a bundle
        # real, 0 < x < 1, small
        */
        this.JOINING_THRESHOLD = joiningThreshold;
        this.NUCLEATION_THRESHOLD = joiningThreshold;
                
        /*
        # Exponent for calculating the generalized mean of signals in 
        # order to find bundle activities
        # real, x != 0
        */
        this.MEAN_EXPONENT = meanExponent;
                
        //# Exponent controlling the strength of inhibition between bundles
        this.ACTIVATION_WEIGHTING_EXPONENT = 6.0;

        this.bundlesFull = false;
                
        this.bundleActivities = new BlockMatrix64F(this.maxBundles, 1);
        
        this.mapSize = new int[]{ maxBundles, maxCables };
               
        this.bundleMap = new BlockMatrix64F(mapSize[0], mapSize[1]);
        this.agglomerationEnergy = new BlockMatrix64F(mapSize[0], mapSize[1]);
        this.nucleationEnergy = new BlockMatrix64F(this.maxCables, 1);                
        
    }
    


    public DenseMatrix64F stepUp(DenseMatrix64F activities) {
        /*
   def step_up(self, cable_activities):
        """ Update co-activity estimates and calculate bundle activity """
        self.cable_activities = cable_activities
        # Find bundle activities by taking the generalized mean of
        # the signals with a negative exponent.
        # The negative exponent weights the lower signals more heavily.
        # At the extreme value of negative infinity, the 
        # generalized mean becomes the minimum operator.
        # Make a first pass at the bundle activation levels by 
        # multiplying across the bundle map.
        initial_bundle_activities = tools.generalized_mean(
                self.cable_activities, self.bundle_map.T, self.MEAN_EXPONENT)
        bundle_contribution_map = np.zeros(self.bundle_map.shape)
        bundle_contribution_map[np.nonzero(self.bundle_map)] = 1.
        activated_bundle_map = (initial_bundle_activities * 
                                bundle_contribution_map)
        # Find the largest bundle activity that each input contributes to
        max_activation = np.max(activated_bundle_map, axis=0) + tools.EPSILON
        # Divide the energy that each input contributes to each bundle
        input_inhibition_map = np.power(activated_bundle_map / max_activation, 
                                        self.ACTIVATION_WEIGHTING_EXPONENT)
        # Find the effective strength of each cable to each bundle 
        # after inhibition.
        inhibited_cable_activities = (input_inhibition_map * 
                                      self.cable_activities.T)
        final_bundle_activities = tools.generalized_mean(
                inhibited_cable_activities.T, self.bundle_map.T, 
                self.MEAN_EXPONENT)
        self.bundle_activities = final_bundle_activities
        # Calculate how much energy each input has left to contribute 
        # to the co-activity estimate. 
        final_activated_bundle_map = (final_bundle_activities * 
                                      bundle_contribution_map)
        combined_weights = np.sum(final_activated_bundle_map, 
                                  axis=0)[:,np.newaxis]
        self.nonbundle_activities = np.maximum(0., (cable_activities - 
                                                    combined_weights))
        # As appropriate update the co-activity estimate and 
        # create new bundles
        if not self.bundles_full:
            self._create_new_bundles()
        self._grow_bundles()
        return self.bundle_activities[:self.num_bundles,:]
        */
        //return activities.getMatrix(0,numBundles, 0,1);
        if (activities.getNumRows() > numBundles)
            return CommonOps.extract(activities, 0, 1, 0, numBundles);
        return activities;        
    }

    public BlockMatrix64F stepDown(BlockMatrix64F goals) {
        /*
    def step_down(self, bundle_goals):
        """ 
        Project the bundle goal values to the appropriate cables

        Multiply the bundle goals across the cables that contribute 
        to them, and perform a bounded sum over all bundles to get 
        the estimated activity associated with each cable.
        """
        if bundle_goals.size > 0:
            bundle_goals = tools.pad(bundle_goals, (self.max_num_bundles, 0))
            cable_activity_goals = tools.bounded_sum(self.bundle_map * 
                                                     bundle_goals, axis=0)
        else:
            cable_activity_goals = np.zeros((self.max_num_cables, 1))
        return cable_activity_goals
                */
        return goals;
    }
    
/*
    def _create_new_bundles(self):
        """ If the right conditions have been reached, create a new bundle """
        # Bundle space is a scarce resource
        # Decay the energy        
        self.nucleation_energy -= (self.cable_activities *
                                   self.nucleation_energy * 
                                   self.NUCLEATION_ENERGY_RATE * 
                                   self.ENERGY_DECAY_RATE)
        self.nucleation_energy += (self.nonbundle_activities * 
                                   (1. - self.nucleation_energy) *
                                   self.NUCLEATION_ENERGY_RATE)
        cable_indices = np.where(self.nucleation_energy > 
                                 self.NUCLEATION_THRESHOLD)
        # Add a new bundle if appropriate
        if cable_indices[0].size > 0:
            # Randomly pick a new cable from the candidates, 
            # if there is more than one
            cable_index = cable_indices[0][int(np.random.random_sample() * 
                                                  cable_indices[0].size)]
            self.bundle_map[self.num_bundles, cable_index] = 1.
            self.num_bundles += 1
            if self.num_bundles == self.max_num_bundles:
                self.bundles_full = True
            print self.name, 'ci', cable_index, 'added as a bundle nucleus'
            self.nucleation_energy[cable_index, 0] = 0.
            self.agglomeration_energy[:, cable_index] = 0.
        return 
          
    def _grow_bundles(self):
        """ Update an estimate of co-activity between all cables """
        coactivities = np.dot(self.bundle_activities, 
                              self.nonbundle_activities.T)
        # Each cable's nonbundle activity is distributed to agglomeration energy
        # with each bundle proportionally to their coactivities.
        proportions_by_bundle = (self.bundle_activities / 
                                 np.sum(self.bundle_activities + tools.EPSILON))
        proportions_by_cable = np.dot(proportions_by_bundle, 
                                      self.nonbundle_activities.T)
        # Decay the energy        
        self.agglomeration_energy -= (proportions_by_cable * 
                                      self.cable_activities.T *
                                      self.agglomeration_energy * 
                                      self.AGGLOMERATION_ENERGY_RATE * 
                                      self.ENERGY_DECAY_RATE)
        self.agglomeration_energy += (proportions_by_cable * 
                                      coactivities * 
                                      (1. - self.agglomeration_energy) *
                                      self.AGGLOMERATION_ENERGY_RATE)
        # For any bundles that are already full, don't change their coactivity
        # TODO: make this more elegant than enforcing a hard maximum count
        full_bundles = np.zeros((self.max_num_bundles, 1))
        cables_per_bundle = np.sum(self.bundle_map, axis=1)[:,np.newaxis]
        full_bundles[np.where(cables_per_bundle >= 
                              self.max_cables_per_bundle)] = 1.
        self.agglomeration_energy *= 1 - full_bundles
        new_candidates = np.where(self.agglomeration_energy >= 
                                  self.JOINING_THRESHOLD)
        num_candidates =  new_candidates[0].size 
        if num_candidates > 0:
            candidate_index = np.random.randint(num_candidates) 
            candidate_cable = new_candidates[1][candidate_index]
            candidate_bundle = new_candidates[0][candidate_index]
            self.bundle_map[candidate_bundle, candidate_cable] = 1.
            self.nucleation_energy[candidate_cable, 0] = 0.
            self.agglomeration_energy[:, candidate_cable] = 0.
            print self.name, 'cable', candidate_cable, 'added to bundle', \
                    candidate_bundle
        return
        
        
    def get_index_projection(self, bundle_index):
        """ Project bundle indices down to their cable indices """
        bundle = np.zeros((self.max_num_bundles, 1))
        bundle[bundle_index, 0] = 1.
        projection = np.sign(np.max(self.bundle_map * bundle, 
                                    axis=0))[np.newaxis, :]
        return projection
        
    def cable_fraction_in_bundle(self, bundle_index):
        cable_count = np.nonzero(self.bundle_map[bundle_index,:])[0].size
        cable_fraction = float(cable_count) / float(self.max_cables_per_bundle)
        return cable_fraction

    def visualize(self, save_eps=False):
        """ Show the internal state of the map in a pictorial format """
        print self.name, '0', np.nonzero(self.bundle_map)[0]
        print self.name, '1', np.nonzero(self.bundle_map)[1]
        print self.max_num_bundles, 'bundles maximum'
        return    
    */   

    public String toString() {
        StringBuffer b = new StringBuffer();
        b.append(super.toString() + "{" + this.bundleMap + "}");
        return b.toString();
    }

    public int getNumBundles() {
        return numBundles;
    }

    
}
