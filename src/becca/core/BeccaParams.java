package becca.core;

public class BeccaParams {
    
    //General
    static int stmSize = 10; //short-term memory size, in cycles (time)
    
    static final boolean RandomGaussian = true;  //gaussian is more computation expensive
    
    static final int NoiseFactorySize = 8192; //pre-allocated random numbers, traversed cyclically.  set to zero to use pure random numbers instead.
    
    static final boolean ExactZero = false;
    
    static final double EpsilonScale = 16.0; //number of times floating point epsilon to use as a minimal threshold.  larger value increases epsilon value threshold, decreasing accuracy / resolution of minimum floating point values. "Precision.EPSILON = Largest double-precision floating-point number such that 1 + EPSILON is numerically equal to 1."        

    
    //Agent ---------------------------------------------------    
    static int agentRecentSurpriseHistorySize = stmSize;
    
    //Block ---------------------------------------------------        
    static boolean BlockGoalBoundedSum = false;
    static int blockMaxCablesPerCog = 16;
    static int blockMaxBundlesPerCog = 8;    
    static double blockFillFractionThreshold = 0;
    static double blockRangeDecayRate =  decayRate(stmSize*512); //original: 0.001
    
    static double blockActivityDecayRate = 1.0; //0.99; //real, 0 < x < 1, higher decays faster
    
    /**
     * lower value makes it easier for the top block to spawn a new top block.
     * when the proportion of bundles nuclated exceeds this threshold, 
     * adds a new top block.
     */
    static double blockInitializationThreshold = 0.5;
    
    
    //Cog --------------------------
    static boolean cogParallel = false;

    
    //DaisyChain ---------------------------------------------------        
    static double daisyCountDecayRate = decayRate(stmSize*2); //real, 0 < x < 1; higher = decays more quickly
    static double daisyChainUpdateRate = decayRate(stmSize*2); //originally: 0.01
    static boolean daisyAllowSelfTransitions = false; //originally: true
    static double daisyAgingTimeConstant = Math.pow(10, 6);
    
    
    //Hub -----------------------------------
    static int hubTraceLength = stmSize;
    static double hubInitialReward = 0.0;
    static double hubUpdateRate = decayRate(stmSize*2);
    
    static double hubRewardDecayRate = decayRate(stmSize); //Recalculate in terms of how many cycles before reward at that history is insignificant (< epsilon).  this cycle time should be approximately equal to hubTraceLength, or can be exactly that
    
    static double hubForgettingRate = decayRate(stmSize*10);
    static double hubExploration = .1;
    
    //ZipTie ----------------------------------------
    static double ziptieSpeedUp = 1.0;
    static int ziptieMeanExponent = -2;
    static double ziptieActivatedBundlemapNoise = 0.05;
    
    //      ZipTie, in Block --------------
    static double blockziptieNucleationThreshold = 0.1;
    static double blockziptieNucleationEnergyRate = 1E-5;
    static double blockziptieAgglomerationThreshold = 0.1;
    static double blockziptieAgglomerationEnergyRate = 1E-3;
    static int blockziptieActivationWeightingExponent = 8;
    
    //      ZipTie, in Cog -----------------
    static double cogziptieNucleationThreshold = 0.1;
    static double cogziptieNucleationEnergyRate = 1E-5;
    static double cogziptieAgglomerationThreshold = 0.1;
    static double cogziptieAgglomerationEnergyRate = 1E-4;
    static int cogziptieActivationWeightingExponent = 8;
    
    
    
   final static double DecayEpsilon = 0.0001; // value for which the result would be insignificant
   
    public static double decayRate(int cycles) {
        
        return 1.0 - Math.pow(DecayEpsilon, (1.0/((double)cycles)));
    }
    
    
}
