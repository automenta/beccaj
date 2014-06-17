package becca.core;

public class BeccaParams {
    
    //Agent ---------------------------------------------------    
    static int agentRecentSurpriseHistorySize = 100;
    
    //Block ---------------------------------------------------        
    
    static int blockMaxCablesPerCog = 16;
    static int blockMaxBundlesPerCog = 8;    
    static double blockFillFractionThreshold = 0.5;
    static double blockRangeDecayRate =  Math.pow(10, -3);
    static double blockActivityDecayRate = 0.5; //real, 0 < x < 1, higher decays faster
    
    //DaisyChain ---------------------------------------------------        
    static double daisyCountDecayRate = 0.1; //real, 0 < x < 1; higher = decays more quickly
    static double daisyChainUpdateRate = 0.1; //originally: 0.01
    static boolean daisyAllowSelfTransitions = false; //originally: true
    
    
    //Hub -----------------------------------
    static int hubTraceLength = 24;
    static double hubInitialReward = 0.0;
    static double hubUpdateRate = Math.pow(10, -3);
    
    static double hubRewardDecayRate = .6; //Recalculate in terms of how many cycles before reward at that history is insignificant (< epsilon).  this cycle time should be approximately equal to hubTraceLength, or can be exactly that
    
    static double hubForgettingRate = Math.pow(10, -3);
    static double hubExploration = .05;
    
    //ZipTie ----------------------------------------
    static double ziptieJoiningThreshold = 0.05;
    static double ziptieSpeedUp = 1.0;
    static double ziptieMeanExponent = -4;
    static double ziptieActivatedBundlemapNoise = 0.01;
    
    //      ZipTie, in Block --------------
    static double blockziptieNucleationThreshold = 0.1;
    static double blockziptieNucleationEnergyRate = 1E-5;
    static double blockziptieAgglomerationThreshold = 0.1;
    static double blockziptieAgglomerationEnergyRate = 1E-3;
    static double blockziptieActivationWeightingExponent = 5;
    
    //      ZipTie, in Cog -----------------
    static double cogziptieNucleationThreshold = 0.1;
    static double cogziptieNucleationEnergyRate = 1E-5;
    static double cogziptieAgglomerationThreshold = 0.1;
    static double cogziptieAgglomerationEnergyRate = 1E-4;
    static double cogziptieActivationWeightingExponent = 5;
    
    
    
}
