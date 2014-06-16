package becca.core;

public class BeccaParams {
    
    //Agent ---------------------------------------------------    
    static int agentRecentSurpriseHistorySize = 100;
    
    //Block ---------------------------------------------------        
    static double blockActivityDecayRate = 0.1; //real, 0 < x < 1, higher decays faster
    static int blockMaxCablesPerCog = 16;
    static int blockMaxBundlesPerCog = 4;    
    static double blockFillFractionThreshold = 0.7;
    static double blockRangeDecayRate =  Math.pow(10, -3);
    
    //DaisyChain ---------------------------------------------------        
    static double daisyCountDecayRate = 0.1; //real, 0 < x < 1; higher = decays more quickly
    static double daisyChainUpdateRate = 0.05; //originally: 0.01
    static boolean daisyAllowSelfTransitions = true; //originally: true
    
    
    //Hub -----------------------------------
    static double hubInitialReward = 1.0;
    static double hubUpdateRate = Math.pow(10, -2);
    static double hubRewardDecayRate = .7;
    static double hubForgettingRate = Math.pow(10, -5);
    static int hubTraceLength = 16;
    static double hubExploration = .1;
    
    //ZipTie ----------------------------------------
    static double ziptieJoiningThreshold = 0.05;
    static double ziptieSpeedUp = 1.0;
    static double ziptieMeanExponent = -4;
    static double ziptieActivatedBundlemapNoise = 0.001;
    
    //      ZipTie, in Block --------------
    static double blockziptieNucleationThreshold = 0.01;
    static double blockziptieNucleationEnergyRate = 1E-5;
    static double blockziptieAgglomerationThreshold = 0.01;
    static double blockziptieAgglomerationEnergyRate = 1E-3;
    static double blockziptieActivationWeightingExponent = 6;
    
    //      ZipTie, in Cog -----------------
    static double cogziptieNucleationThreshold = 0.1;
    static double cogziptieNucleationEnergyRate = 1E-5;
    static double cogziptieAgglomerationThreshold = 0.1;
    static double cogziptieAgglomerationEnergyRate = 1E-4;
    static double cogziptieActivationWeightingExponent = 6;
    
    
    
}
