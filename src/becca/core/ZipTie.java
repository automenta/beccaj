package becca.core;

import org.ejml.data.DenseMatrix64F;

import static becca.core.Util.*;
import java.util.ArrayList;
import java.util.List;
import org.ejml.data.RowD1Matrix64F;



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
    private final double AGGLOMERATION_THRESHOLD;
    private final double NUCLEATION_ENERGY_RATE;
    private final double NUCLEATION_THRESHOLD;
    private final double MEAN_EXPONENT;
    private final double ACTIVATION_WEIGHTING_EXPONENT;
    private double ACTIVATED_BUNDLE_MAP_NOISE;
    
    private final DenseMatrix64F bundleMap;
    private DenseMatrix64F bundleActivities;
    private DenseMatrix64F agglomerationEnergy;
    private final DenseMatrix64F nucleationEnergy;
    private DenseMatrix64F cableActivities;
    private DenseMatrix64F nonBundleActivities;
    private boolean bundlesFull;
    private final boolean inBlock;

    private List<Integer> cableIndices = new ArrayList();
    private List<int[]> newCandidates = new ArrayList();
    
    
    private DenseMatrix64F activatedBundleMap;
    private DenseMatrix64F inputInhibitionMap;
    private DenseMatrix64F coactivities;
    private DenseMatrix64F aggMinus;
    private DenseMatrix64F aggPlus;

    public ZipTie(boolean inBlock, int maxCables, int maxBundles, int maxCablesPerBundle) {
        this(inBlock, maxCables, maxBundles, maxCablesPerBundle, BeccaParams.ziptieMeanExponent);
    }

    public ZipTie(boolean inBlock, int maxCables, int maxBundles, int maxCablesPerBundle, double meanExponent) {
        this(inBlock, maxCables, maxBundles, maxCablesPerBundle, meanExponent, BeccaParams.ziptieSpeedUp);
    }
    
    public ZipTie(boolean inBlock, int maxCables, int maxBundles, int maxCablesPerBundle, double meanExponent, double speedup) {
        
        this.maxCables = maxCables;
        this.maxBundles = maxBundles;
        this.inBlock = inBlock;
        
        if (maxCablesPerBundle == 0)
            maxCablesPerBundle = (int)( ((double)maxCables) / ((double)maxBundles) );
        
        this.maxCablesPerBundle = maxCablesPerBundle;
        
        this.numBundles = 0;
        this.ACTIVATED_BUNDLE_MAP_NOISE = BeccaParams.ziptieActivatedBundlemapNoise;
        

        //Identify whether the ziptie is located in a block or in a cog.
        //This changes several aspects of its operation.         
        if (inBlock) {
            //These seem to strike a nice balance between feature quality and learning speed
            this.NUCLEATION_THRESHOLD = BeccaParams.blockziptieNucleationThreshold;
            this.NUCLEATION_ENERGY_RATE = BeccaParams.blockziptieNucleationEnergyRate;
            this.AGGLOMERATION_THRESHOLD = BeccaParams.blockziptieAgglomerationThreshold;
            this.AGGLOMERATION_ENERGY_RATE = BeccaParams.blockziptieAgglomerationEnergyRate;
            this.ACTIVATION_WEIGHTING_EXPONENT = BeccaParams.blockziptieActivationWeightingExponent;
        }
        else {
            //For zipties within cogs
            this.NUCLEATION_THRESHOLD = BeccaParams.cogziptieNucleationThreshold;
            this.NUCLEATION_ENERGY_RATE = BeccaParams.cogziptieNucleationEnergyRate;
            this.AGGLOMERATION_THRESHOLD = BeccaParams.cogziptieAgglomerationThreshold;
            this.AGGLOMERATION_ENERGY_RATE = BeccaParams.cogziptieAgglomerationEnergyRate;
            this.ACTIVATION_WEIGHTING_EXPONENT = BeccaParams.cogziptieActivationWeightingExponent;            
        }
        
        /*
        # Exponent for calculating the generalized mean of signals in 
        # order to find bundle activities
        # real, x != 0
        */
        this.MEAN_EXPONENT = meanExponent;

        //Exponent controlling the strength of inhibition between bundles
        
              
        this.bundlesFull = false;
                
                
        this.bundleActivities = new DenseMatrix64F(this.maxBundles, 1);

        assert(maxCables > 0);
        this.bundleMap = new DenseMatrix64F(maxBundles, maxCables);        

        this.agglomerationEnergy = new DenseMatrix64F(maxBundles, maxCables);
        
        this.nucleationEnergy = new DenseMatrix64F(this.maxCables, 1);
    }
    

    public DenseMatrix64F stepUp(DenseMatrix64F cableActivitiesIn) {
        // Update co-activity estimates and calculate bundle activity """

        this.cableActivities = cableActivitiesIn;

        /*DenseMatrix64F cableActivitiesTranspose = transpose(cableActivities, null);*/
        
        /*
        # Find bundle activities by taking the generalized mean of
        # the signals with a negative exponent.
        # The negative exponent weights the lower signals more heavily.
        # At the extreme value of negative infinity, the 
        # generalized mean becomes the minimum operator.
        # Make a first pass at the bundle activation levels by 
        # multiplying across the bundle map.
         */
        DenseMatrix64F bundleMapT = transpose(bundleMap, null);
        
        
        /*initial_bundle_activities = tools.generalized_mean(
                self.cable_activities, self.bundle_map.T, self.MEAN_EXPONENT)*/
        DenseMatrix64F initialBundleActivities = getGeneralizedMean(cableActivities, bundleMapT, MEAN_EXPONENT); //WAS transposed before
        
                
        //bundle_contribution_map = np.zeros(self.bundle_map.shape)
        //bundle_contribution_map[np.nonzero(self.bundle_map)] = 1
        //TODO do not reallocate this each cycle
        DenseMatrix64F bundleContributionMap = getNonZeroMask(bundleMap);
        
        /*
        # use aggressive lateral 
        # inhibition between bundles so that cables' activity
        # is monopolized by the strongest-activated bundle
        activated_bundle_map = (initial_bundle_activities * 
                                bundle_contribution_map)
        */
        activatedBundleMap = matrixVector(bundleContributionMap, initialBundleActivities, true, activatedBundleMap);
        
        /*
        # Add just a little noise to break ties
        
        activated_bundle_map += .0001 * np.random.random_sample(
                activated_bundle_map.shape)        
        */
        matrixAddNoise(activatedBundleMap, ACTIVATED_BUNDLE_MAP_NOISE);
        
        /*
        # Find the largest bundle activity that each input contributes to
        max_activation = np.max(activated_bundle_map, axis=0) + tools.EPSILON
         */
        //maximum value along row axis
        DenseMatrix64F maxActivation = maxRow(activatedBundleMap);
        assert(maxActivation.getNumCols() == activatedBundleMap.getNumCols());
        add(maxActivation, Util.EPSILON);
                
        /*
        # Divide the energy that each input contributes to each bundle
        input_inhibition_map = np.power(activated_bundle_map / max_activation, 
                                        self.ACTIVATION_WEIGHTING_EXPONENT)
        */                
        inputInhibitionMap = matrixVector(activatedBundleMap, maxActivation, false, null);
        matrixPower(inputInhibitionMap, ACTIVATION_WEIGHTING_EXPONENT);

        //# Find the effective strength of each cable to each bundle after inhibition.
        //inhibited_cable_activities = (input_inhibition_map * self.cable_activities.T)
        
        DenseMatrix64F inhibitedCableActivities = inhibitedCableActivities = matrixVector(inputInhibitionMap, cableActivities/*Transpose*/);               
        /*DenseMatrix64F inhibitedCableActivitiesT = transpose(inhibitedCableActivities, null);*/
        
        /*final_bundle_activities = tools.generalized_mean(inhibited_cable_activities.T, self.bundle_map.T, self.MEAN_EXPONENT)*/        
        //self.bundle_activities = final_bundle_activities
        

        assert(inhibitedCableActivities.getNumRows() == bundleMapT.getNumRows());
        DenseMatrix64F finalBundleActivities = this.bundleActivities = 
                transpose(getGeneralizedMean(inhibitedCableActivities, bundleMapT, MEAN_EXPONENT),null);
                
        
        //# Calculate how much energy each input has left to contribute to the co-activity estimate. 
        //final_activated_bundle_map = (final_bundle_activities * bundle_contribution_map)
        

        //combined_weights = np.sum(final_activated_bundle_map, axis=0)[:,np.newaxis]
        
        if (this.inBlock) {
            /*self.nonbundle_activities 
                = (cable_activities * 
                    2 ** -np.sum(self.bundle_map, 
                    axis=0)[:,np.newaxis])     */
            nonBundleActivities = sumColsT(bundleMap, null, true);

            scale(-1, nonBundleActivities);
            matrixPowerExp(nonBundleActivities, 2);             
            elementMult(nonBundleActivities, cableActivities);
        }
        else {
            /* self.nonbundle_activities = np.maximum(0., cable_activities - combined_weights)    */
            DenseMatrix64F finalActivatedBundleMap = matrixVector(bundleMap, finalBundleActivities);
            DenseMatrix64F combinedWeights = sumColsT(finalActivatedBundleMap, null, true);

            assert(combinedWeights.getNumRows()==cableActivities.getNumRows());
            assert(combinedWeights.getNumCols()==1);
            
            this.nonBundleActivities = cableActivities.copy();
            
            subEquals(this.nonBundleActivities, combinedWeights);
            matrixMaximum(this.nonBundleActivities, 0);            
        }
        this.cableActivities = cableActivitiesIn;

        
        //# As appropriate update the co-activity estimate and create new bundles
        if (!bundlesFull) {
            createNewBundles();
        }
        growBundles();
                 
        
        //return self.bundle_activities[:self.num_bundles,:]
        if (numBundles > 0)
            return extract(cableActivities, 0, numBundles, 0, cableActivities.getNumCols() );
        else
            return new DenseMatrix64F(0, cableActivities.getNumCols());
    }

    
    protected void createNewBundles() {
        
        //""" If the right conditions have been reached, create a new bundle """
        /*
        # Bundle space is a scarce resource
        # Decay the energy        
        self.nucleation_energy -= (self.cable_activities *
                                   self.nucleation_energy * 
                                   self.NUCLEATION_ENERGY_RATE) 
        */
        DenseMatrix64F nucleationEnergyDelta = nucleationEnergy.copy();
        elementMult(nucleationEnergyDelta, cableActivities);        
        scale(NUCLEATION_ENERGY_RATE, nucleationEnergyDelta);
        subEquals(nucleationEnergy, nucleationEnergyDelta);
        
        //--
        /*
        self.nucleation_energy += (self.nonbundle_activities * 
                                   (1. - self.nucleation_energy) *
                                   self.NUCLEATION_ENERGY_RATE)
        */
        DenseMatrix64F oneMinusNE = nucleationEnergy.copy();
        scale(-1, oneMinusNE);
        add(oneMinusNE, 1.0);        
                
        DenseMatrix64F nucleationEnergyDelta2 = oneMinusNE.copy();        
        elementMult(nucleationEnergyDelta2, nonBundleActivities);
        scale(NUCLEATION_ENERGY_RATE, nucleationEnergyDelta2);
        
        addEquals(nucleationEnergy, nucleationEnergyDelta2);
                        
        
        
        /* cable_indices = np.where(self.nucleation_energy > 
                                 self.NUCLEATION_THRESHOLD)  */           
        cableIndices.clear();
        final double[] ned = nucleationEnergy.getData();
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
                this.bundlesFull = true;
            
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
        
        //coactivities = np.dot(self.bundle_activities, self.nonbundle_activities.T)               
        coactivities = multM(bundleActivities, nonBundleActivitiesT, coactivities);                
        
        

        
        /*
        self.agglomeration_energy -= (self.cable_activities.T *
                                      self.agglomeration_energy * 
                                      self.AGGLOMERATION_ENERGY_RATE)
        */
        
        
        DenseMatrix64F cableActivitiesT = transpose(cableActivities, null);
        aggMinus = matrixVector(agglomerationEnergy, cableActivitiesT, true, aggMinus);
        scale(AGGLOMERATION_ENERGY_RATE, aggMinus);     
        subEquals(agglomerationEnergy, aggMinus);
        
        /*
        self.agglomeration_energy += (coactivities * 
                                      (1. - self.agglomeration_energy) *
                                      self.AGGLOMERATION_ENERGY_RATE)  */
        aggPlus = ensureSize(aggPlus, agglomerationEnergy.numRows, agglomerationEnergy.numCols);
        aggPlus.set(agglomerationEnergy);
        scale(-1, aggPlus);
        add(aggPlus, 1);        
        elementMult(aggPlus, coactivities);
        scale(AGGLOMERATION_ENERGY_RATE, aggPlus);
        addEquals(agglomerationEnergy, aggPlus);
        

        
        /*
        # For any bundles that are already full, don't change their coactivity                
        # TODO: make this more elegant than enforcing a hard maximum count
        full_bundles = np.zeros((self.max_num_bundles, 1))
        */
        DenseMatrix64F fullBundles = new DenseMatrix64F(maxBundles, 1);
        
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
        scale(-1, fullBundles); //we can modify full_bundles because it is used nowhere else
        add(fullBundles, 1.0);
                
        //TODO use alternate version of matrixVector that computes in-place without allocating
        //and eliminate need for transpose
        agglomerationEnergy = matrixVector(agglomerationEnergy, fullBundles);       
                
        /*
        new_candidates = np.where(self.agglomeration_energy >= self.JOINING_THRESHOLD)                
        num_candidates =  new_candidates[0].size 
        */
        newCandidates.clear();
        for (int i = 0; i < agglomerationEnergy.getNumRows(); i++)
            for (int j = 0; j < agglomerationEnergy.getNumCols(); j++)
                if (agglomerationEnergy.get(i, j) >= AGGLOMERATION_THRESHOLD)
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
        
    
    
    
    public DenseMatrix64F stepDown(DenseMatrix64F bundleGoals) {
        /*    
        """ 
        Project the bundle goal values to the appropriate cables

        Multiply the bundle goals across the cables that contribute 
        to them, and perform a bounded sum over all bundles to get 
        the estimated activity associated with each cable.
        """
        */
        DenseMatrix64F cableActivityGoals;
        if (bundleGoals.getNumElements() > 0) {
            //bundle_goals = tools.pad(bundle_goals, (self.max_num_bundles, 0))
            bundleGoals = pad(bundleGoals, maxBundles, 0, 0.0);
            
            //cable_activity_goals = tools.bounded_sum(self.bundle_map * bundle_goals, axis=0)
            //self.bundle_map * bundle_goals
            DenseMatrix64F mapgoals = matrixVector(bundleMap, bundleGoals);
            assert(mapgoals.getNumRows() == bundleMap.getNumRows());
            assert(mapgoals.getNumCols() == bundleMap.getNumCols());
                                    
            cableActivityGoals = boundedRowSum(mapgoals);            
            //cableActivityGoals = sumCols(mapgoals, null);
        }
        else {
            cableActivityGoals = new DenseMatrix64F(maxCables, 1);
        }
        return cableActivityGoals;
    }

    public DenseMatrix64F getIndexProjection(int bundleIndex) {
        //""" Project bundle indices down to their cable indices """
        
        //projection = self.bundle_map[bundle_index,:]
        return extract(bundleMap,bundleIndex,bundleIndex+1,0,bundleMap.getNumCols());

        ///OLD CODE ----------------
        
        
        /*DenseMatrix64F bundle = new DenseMatrix64F(maxBundles, 1);
        bundle.set(bundleIndex, 0, 1.0);

        //TODO may not need to invert 3 times
        //projection = np.sign(np.max(self.bundle_map * bundle,axis=0))[np.newaxis, :]        
        DenseMatrix64F bmb = transpose(matrixVector(transpose(bundleMap,null), transpose(bundle,null)),null);
        
        //TODO this may need to be iterated by rows, depending on which axis is meant        
        
        DenseMatrix64F projection = maxRow(bmb);
        matrixSign(projection);
                
        return projection; */
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

    double getCableFractionInBundle(final int bundleIndex) {
        /*
        cable_count = np.nonzero(self.bundle_map[bundle_index,:])[0].size
        cable_fraction = float(cable_count) / float(self.max_cables_per_bundle)
        return cable_fraction
        */
        //bundleindex selects a row of bundleMap
        //counts how many non-zero elements are in that row
        int cableCount = 0;
        
        for (int i = 0; i < bundleMap.getNumCols(); i++) {
            if (bundleMap.get(bundleIndex, i) != 0)
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
