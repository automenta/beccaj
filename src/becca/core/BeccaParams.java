package becca.core;

public class BeccaParams {
    
    static int stmSize = 16; //short-term memory size, in cycles (time)

    //gaussian is more computation expensive
    static boolean RandomGaussian = false;
    
    //Agent ---------------------------------------------------    
    static int agentRecentSurpriseHistorySize = stmSize;
    
    //Block ---------------------------------------------------        
    
    static int blockMaxCablesPerCog = 8;
    static int blockMaxBundlesPerCog = 4;    
    static double blockFillFractionThreshold = 0.1;
    static double blockRangeDecayRate =  Math.pow(10, -3);
    static double blockActivityDecayRate = 0.99; //real, 0 < x < 1, higher decays faster
    static double blockInitializationThreshold = 0.5;
    
    
    //Cog --------------------------
    static boolean cogParallel = false;

    
    //DaisyChain ---------------------------------------------------        
    static double daisyCountDecayRate = decayRate(stmSize*2); //real, 0 < x < 1; higher = decays more quickly
    static double daisyChainUpdateRate = decayRate(stmSize*2); //originally: 0.01
    static boolean daisyAllowSelfTransitions = true; //originally: true
    
    
    //Hub -----------------------------------
    static int hubTraceLength = stmSize;
    static double hubInitialReward = 0.0;
    static double hubUpdateRate = decayRate(stmSize*2);
    
    static double hubRewardDecayRate = decayRate(stmSize); //Recalculate in terms of how many cycles before reward at that history is insignificant (< epsilon).  this cycle time should be approximately equal to hubTraceLength, or can be exactly that
    
    static double hubForgettingRate = decayRate(stmSize*2);
    static double hubExploration = .1;
    
    //ZipTie ----------------------------------------
    static double ziptieSpeedUp = 1.0;
    static double ziptieMeanExponent = -2;
    static double ziptieActivatedBundlemapNoise = 0.01;
    
    //      ZipTie, in Block --------------
    static double blockziptieNucleationThreshold = 0.1;
    static double blockziptieNucleationEnergyRate = 1E-5;
    static double blockziptieAgglomerationThreshold = 0.1;
    static double blockziptieAgglomerationEnergyRate = 1E-3;
    static double blockziptieActivationWeightingExponent = 512;
    
    //      ZipTie, in Cog -----------------
    static double cogziptieNucleationThreshold = 0.1;
    static double cogziptieNucleationEnergyRate = 1E-5;
    static double cogziptieAgglomerationThreshold = 0.1;
    static double cogziptieAgglomerationEnergyRate = 1E-4;
    static double cogziptieActivationWeightingExponent = 512;
    
   
    public static double decayRate(int cycles) {
        double epsilon = 0.0001; // value for which the result would be insignificant
        return 1.0 - Math.pow(epsilon, (1.0/((double)cycles)));
    }
    
    
}
