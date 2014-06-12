/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package becca.core;

import org.ejml.data.BlockMatrix64F;
import org.ejml.data.DenseMatrix64F;

import static becca.core.Util.*;
import java.util.ArrayList;
import java.util.List;


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
    private int numBundles;
    private final double AGGLOMERATION_ENERGY_RATE;
    private final double NUCLEATION_ENERGY_RATE;
    private final double ENERGY_DECAY_RATE;
    private final double JOINING_THRESHOLD;
    private final double NUCLEATION_THRESHOLD;
    private final double MEAN_EXPONENT;
    private final double ACTIVATION_WEIGHTING_EXPONENT;
    private boolean bundlesFull;
    
    private final DenseMatrix64F bundleMap;
    private DenseMatrix64F bundleActivities;
    private DenseMatrix64F agglomerationEnergy;
    private final DenseMatrix64F nucleationEnergy;
    private DenseMatrix64F cableActivities;
    private DenseMatrix64F nonBundleActivities;

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
                
        this.bundleActivities = new DenseMatrix64F(this.maxBundles, 1);
        
        
        this.bundleMap = new DenseMatrix64F(maxBundles, maxCables);
        
        
        this.agglomerationEnergy = new DenseMatrix64F(maxBundles, maxCables);
        this.nucleationEnergy = new DenseMatrix64F(this.maxCables, 1);                
        
    }
    


    public DenseMatrix64F stepUp(DenseMatrix64F cableActivities) {
        // Update co-activity estimates and calculate bundle activity """

        this.cableActivities = cableActivities;

        DenseMatrix64F cableActivitiesTranspose = cableActivities.copy();
        transpose(cableActivitiesTranspose);
        
        /*
        # Find bundle activities by taking the generalized mean of
        # the signals with a negative exponent.
        # The negative exponent weights the lower signals more heavily.
        # At the extreme value of negative infinity, the 
        # generalized mean becomes the minimum operator.
        # Make a first pass at the bundle activation levels by 
        # multiplying across the bundle map.
         */
        DenseMatrix64F bundleMapTranspose = bundleMap.copy();
        transpose(bundleMapTranspose);
        
        /*initial_bundle_activities = tools.generalized_mean(
                self.cable_activities, self.bundle_map.T, self.MEAN_EXPONENT)*/
        DenseMatrix64F initialBundleActivities = getGeneralizedMean(cableActivitiesTranspose, bundleMapTranspose, MEAN_EXPONENT);
        
                
        //bundle_contribution_map = np.zeros(self.bundle_map.shape)
        //bundle_contribution_map[np.nonzero(self.bundle_map)] = 1.
        //DenseMatrix64F bundleContributionMap = new DenseMatrix64F(bundleMap.getNumRows(), bundleMap.getNumCols());
        //DenseMatrix64F bundleMapNonZero = Util.getNonZero(bundleMap);
        DenseMatrix64F bundleContributionMap = getNonZeroMask(bundleMap);
        
        //activated_bundle_map = (initial_bundle_activities * bundle_contribution_map)        
        DenseMatrix64F activatedBundleMap = new DenseMatrix64F(initialBundleActivities.getNumRows(), bundleContributionMap.getNumCols());
        mult(initialBundleActivities, bundleContributionMap, activatedBundleMap);
        
        /*
        # Find the largest bundle activity that each input contributes to
        max_activation = np.max(activated_bundle_map, axis=0) + tools.EPSILON
         */
        DenseMatrix64F maxActivation = maxRow(activatedBundleMap);
        add(maxActivation, Util.EPSILON);
                
        /*
        # Divide the energy that each input contributes to each bundle
        input_inhibition_map = np.power(activated_bundle_map / max_activation, 
                                        self.ACTIVATION_WEIGHTING_EXPONENT)
        */
        DenseMatrix64F inputInhibitionMap = matrixDivide(activatedBundleMap, maxActivation);
        matrixPower(inputInhibitionMap, ACTIVATION_WEIGHTING_EXPONENT);

        //# Find the effective strength of each cable to each bundle after inhibition.
        //inhibited_cable_activities = (input_inhibition_map * self.cable_activities.T)
        
        
        
        DenseMatrix64F inhibitedCableActivities = new DenseMatrix64F(inputInhibitionMap.getNumRows(), cableActivitiesTranspose.getNumCols());
        mult(inputInhibitionMap, cableActivitiesTranspose, inhibitedCableActivities);
        
        DenseMatrix64F inhibitedCableActivitiesT = inhibitedCableActivities.copy(); 
        transpose(inhibitedCableActivitiesT);
        
        /*final_bundle_activities = tools.generalized_mean(inhibited_cable_activities.T, self.bundle_map.T, self.MEAN_EXPONENT)*/        
        //self.bundle_activities = final_bundle_activities
        DenseMatrix64F finalBundleActivities = this.bundleActivities = 
                Util.getGeneralizedMean(inhibitedCableActivitiesT, bundleMapTranspose, MEAN_EXPONENT);               
                
        
        //# Calculate how much energy each input has left to contribute to the co-activity estimate. 
        //final_activated_bundle_map = (final_bundle_activities * bundle_contribution_map)
        DenseMatrix64F finalActivatedBundleMap = new DenseMatrix64F(finalBundleActivities.getNumRows(), bundleMap.getNumCols());
        mult(finalBundleActivities, bundleMap, finalActivatedBundleMap);
                
        //combined_weights = np.sum(final_activated_bundle_map, axis=0)[:,np.newaxis]
        DenseMatrix64F combinedWeights = sumCols(finalActivatedBundleMap, null);
        transpose(combinedWeights);
        
        
        //self.nonbundle_activities = np.maximum(0., (cable_activities - combined_weights))
        this.nonBundleActivities = new DenseMatrix64F(cableActivities.getNumRows(), cableActivities.getNumCols());
        sub(cableActivities, combinedWeights, nonBundleActivities);
        matrixMaximum(nonBundleActivities, 0);
        
        //# As appropriate update the co-activity estimate and 
        //# create new bundles
        if (!bundlesFull) {
            createNewBundles();
        }
        growBundles();
                        
        DenseMatrix64F r = DenseMatrix64F.wrap(numBundles, 1, cableActivities.getData());
        return r;
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

    public DenseMatrix64F getIndexProjection(int bundleIndex) {
        //""" Project bundle indices down to their cable indices """
        DenseMatrix64F bundle = new DenseMatrix64F(maxBundles, 1);
        bundle.set(bundleIndex, 0, 1.0);

        //projection = np.sign(np.max(self.bundle_map * bundle,axis=0))[np.newaxis, :]        
        DenseMatrix64F bmb = transpose(matrixVector(transpose(bundleMap,null), transpose(bundle,null)),null);
        
        //TODO this may need to be iterated by rows, depending on which axis is meant        
        
        DenseMatrix64F projection = maxRow(bmb);
        matrixSign(projection);
                
        return projection;
    }
    
    protected void createNewBundles() {
        
        //""" If the right conditions have been reached, create a new bundle """
        /*
        # Bundle space is a scarce resource
        # Decay the energy        
        */
        
        /*self.nucleation_energy -= (self.cable_activities *
                                   self.nucleation_energy * 
                                   self.NUCLEATION_ENERGY_RATE * 
                                   self.ENERGY_DECAY_RATE)       */ 
        

        DenseMatrix64F nucleationEnergyDelta = new DenseMatrix64F(cableActivities.getNumRows(), nucleationEnergy.getNumCols());
        elementMult(cableActivities, nucleationEnergy, nucleationEnergyDelta);
        scale(NUCLEATION_ENERGY_RATE * ENERGY_DECAY_RATE, nucleationEnergyDelta);
        
        sub(nucleationEnergy, nucleationEnergyDelta, nucleationEnergy);
        
        //--
        /*
        self.nucleation_energy += (self.nonbundle_activities * 
                                   (1. - self.nucleation_energy) *
                                   self.NUCLEATION_ENERGY_RATE)
        */
        
        DenseMatrix64F nucleationEnergyDelta2 = new DenseMatrix64F(nonBundleActivities.getNumRows(), nucleationEnergy.getNumCols());
        
        DenseMatrix64F nucleationEnergyOnes = new DenseMatrix64F(nucleationEnergy.getNumRows(), nucleationEnergy.getNumCols());
        fill(nucleationEnergyOnes, 1.0);
        sub(nucleationEnergyOnes, nucleationEnergy, nucleationEnergyOnes);
        
        
        elementMult(nonBundleActivities, nucleationEnergyOnes, nucleationEnergyDelta2);
        scale(NUCLEATION_ENERGY_RATE, nucleationEnergyDelta2);
        
        add(nucleationEnergy, nucleationEnergyDelta2, nucleationEnergy);
        
        /* cable_indices = np.where(self.nucleation_energy > 
                                 self.NUCLEATION_THRESHOLD)  */
        List<Integer> cableIndices = new ArrayList();
        double[] ned = nucleationEnergy.getData();
        for (int i = 0; i < ned.length; i++)
            if (ned[i] > NUCLEATION_THRESHOLD)
                cableIndices.add(i);
                
        //# Add a new bundle if appropriate
        if (cableIndices.size() > 0) {
        
            //# Randomly pick a new cable from the candidates, if there is more than one
            
            /*cable_index = cable_indices[0][int(np.random.random_sample() * 
                                                  cable_indices[0].size)]*/
            int cableIndex = (int)(Math.random() * cableIndices.size());

            //self.bundle_map[self.num_bundles, cable_index] = 1.
            bundleMap.set(numBundles, cableIndex, 1.0);
            
            numBundles++;
            
            if (numBundles == maxBundles)
                bundlesFull = true;
                                
            //print self.name, 'ci', cable_index, 'added as a bundle nucleus'
            
            //self.nucleation_energy[cable_index, 0] = 0.
            nucleationEnergy.set(cableIndex, 0, 0.0);
            
            //self.agglomeration_energy[:, cable_index] = 0.
            for (int i = 0; i < agglomerationEnergy.getNumRows(); i++)
                agglomerationEnergy.set(i, cableIndex, 0.0);
            
            
        }
                
    }
    
    protected void growBundles() {
        //""" Update an estimate of co-activity between all cables """
       
        DenseMatrix64F nonBundleActivitiesT = transpose(nonBundleActivities, null);
        DenseMatrix64F bundleActivitiesT = transpose(bundleActivities, null);
        
        //coactivities = np.dot(self.bundle_activities, self.nonbundle_activities.T)       
        
        DenseMatrix64F coactivities = new DenseMatrix64F(bundleActivitiesT.getNumRows(), nonBundleActivitiesT.getNumCols());
        mult(bundleActivitiesT, nonBundleActivitiesT, coactivities);                
        
        /* Each cable's nonbundle activity is distributed to agglomeration energy
           with each bundle proportionally to their coactivities.

        proportions_by_bundle = (self.bundle_activities / 
                                 np.sum(self.bundle_activities + tools.EPSILON))   */
        DenseMatrix64F proportionsByBundle = bundleActivities.copy();
        add(proportionsByBundle, EPSILON);
        double sum = elementSum(proportionsByBundle);
        proportionsByBundle.set(bundleActivities);
        scale(1.0 / sum, proportionsByBundle);
                                
        /*
        proportions_by_cable = np.dot(proportions_by_bundle, 
                                      self.nonbundle_activities.T)        */

        DenseMatrix64F proportionsByBundleT = transpose(proportionsByBundle, null);
        DenseMatrix64F proportionsByCable = new DenseMatrix64F(proportionsByBundleT.getNumRows(), nonBundleActivitiesT.getNumCols());
        mult(proportionsByBundleT, nonBundleActivitiesT, proportionsByCable);

        
        /*
        # Decay the energy        
        self.agglomeration_energy -= (proportions_by_cable * 
                                      self.cable_activities.T *
                                      self.agglomeration_energy * 
                                      self.AGGLOMERATION_ENERGY_RATE * 
                                      self.ENERGY_DECAY_RATE)
        */
        //(8, 32) * (1, 32) -> (8, 32)
        //(8, 32) * (8, 32) => (8, 32)

        
        DenseMatrix64F cableActivitiesT = transpose(cableActivities, null);
        DenseMatrix64F proportionsByCableT = transpose(proportionsByCable, null);
        
        DenseMatrix64F aggDelta0 = matrixVector(proportionsByCable, cableActivitiesT);
     
        DenseMatrix64F aggDelta = new DenseMatrix64F(aggDelta0.getNumRows(), agglomerationEnergy.getNumCols());        
        elementMult(aggDelta0, agglomerationEnergy, aggDelta);
        
        
        scale(AGGLOMERATION_ENERGY_RATE * ENERGY_DECAY_RATE, aggDelta);
        

        sub(agglomerationEnergy, aggDelta, agglomerationEnergy);
        
        /*
        self.agglomeration_energy += (proportions_by_cable * 
                                      coactivities * 
                                      (1. - self.agglomeration_energy) *
                                      self.AGGLOMERATION_ENERGY_RATE)
        */
        DenseMatrix64F oneMinusAgglom = new DenseMatrix64F(agglomerationEnergy.getNumRows(), agglomerationEnergy.getNumCols());
        fill(oneMinusAgglom,1.0);
        
        sub(oneMinusAgglom, agglomerationEnergy, oneMinusAgglom);
                
        DenseMatrix64F aggDelta3 = matrixVector(proportionsByCable, coactivities);

        DenseMatrix64F aggDelta4 = new DenseMatrix64F(aggDelta3.getNumRows(), oneMinusAgglom.getNumCols());
        elementMult(aggDelta3, oneMinusAgglom, aggDelta4);
        
        scale(AGGLOMERATION_ENERGY_RATE, aggDelta4);
        
        add(agglomerationEnergy, aggDelta4, agglomerationEnergy);
        
        
        
        /*
        # For any bundles that are already full, don't change their coactivity                
        # TODO: make this more elegant than enforcing a hard maximum count
        full_bundles = np.zeros((self.max_num_bundles, 1))
        */
        DenseMatrix64F fullBundles = new DenseMatrix64F(1, maxBundles);
        
        /*
        cables_per_bundle = np.sum(self.bundle_map, axis=1)[:,np.newaxis]
        */
        DenseMatrix64F cablesPerBundle = sumRows(bundleMap, null);
        
        /*
        full_bundles[np.where(cables_per_bundle >= 
                              self.max_cables_per_bundle)] = 1.
        */
        double[] cpbD = cablesPerBundle.getData();
        double[] fbD = fullBundles.getData();
        for (int i = 0; i < maxBundles; i++)
            if (cpbD[i] >= maxCablesPerBundle)
                fbD[i] = 1.0;
        
        //self.agglomeration_energy *= 1 - full_bundles
        scale(-1, fullBundles);
        add(fullBundles, 1.0);
                
        //TODO use alternate version of matrixVector that computes in-place without allocating
        //and eliminate need for transpose
        agglomerationEnergy = transpose(matrixVector(transpose(agglomerationEnergy,null), fullBundles), null);       
                
        /*
        new_candidates = np.where(self.agglomeration_energy >= self.JOINING_THRESHOLD)                
        num_candidates =  new_candidates[0].size 
        */
        List<int[]> newCandidates = new ArrayList();
        for (int i = 0; i < agglomerationEnergy.getNumRows(); i++)
            for (int j = 0; j < agglomerationEnergy.getNumCols(); j++)
                if (agglomerationEnergy.get(i, j) >= JOINING_THRESHOLD)
                    newCandidates.add(new int[]{i,j});
        
        int numCandidates = newCandidates.size();
        
        if (numCandidates > 0) {
            //candidate_index = np.random.randint(num_candidates) 
            int candidateIndex = (int)(Math.random() * numCandidates);
            
            //candidate_cable = new_candidates[1][candidate_index]
            //candidate_bundle = new_candidates[0][candidate_index]
            int candidateCable = newCandidates.get(candidateIndex)[1];
            int candidateBundle = newCandidates.get(candidateIndex)[0];
            
            //self.bundle_map[candidate_bundle, candidate_cable] = 1.
            assert(candidateBundle < bundleMap.getNumRows());
            assert(candidateCable < bundleMap.getNumCols());
            bundleMap.set(candidateBundle, candidateCable, 1.0);
            
            //self.nucleation_energy[candidate_cable, 0] = 0.
            nucleationEnergy.set(candidateCable, 0, 0.0);           
            
            //self.agglomeration_energy[:, candidate_cable] = 0.
            for (int i = 0; i < agglomerationEnergy.getNumRows(); i++) {
                agglomerationEnergy.set(i, candidateCable, 0.0);
            }
            
            //print self.name, 'cable', candidate_cable, 'added to bundle', candidate_bundle
            
        }
                    
    }
    
        
        
    /*    

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

    double getCableFractionInBundle(int bundleIndex) {

        int cableCount = 0;
        for (int i = 0; i < bundleMap.getNumCols(); i++) {
            if (bundleMap.get(bundleIndex, i) == 0)
                cableCount++;
        }
        return ((double)cableCount) / ((double)maxCablesPerBundle);
        
    }

    public DenseMatrix64F getBundleActivities() {
        return bundleActivities;
    }

    public DenseMatrix64F getCableActivities() {
        return cableActivities;
    }

    public DenseMatrix64F getBundleMap() {
        return bundleMap;
    }

    public DenseMatrix64F getAgglomerationEnergy() {
        return agglomerationEnergy;
    }

    public DenseMatrix64F getNucleationEnergy() {
        return nucleationEnergy;
    }

    public DenseMatrix64F getNonBundleActivities() {
        return nonBundleActivities;
    }
    
    
    
    
}
