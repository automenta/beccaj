/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package conceptor;

import static becca.core.Util.broadcastRows;
import static becca.core.Util.matrixRandomGaussian;
import static becca.core.Util.matrixRandomUniform;
import static becca.core.Util.matrixVector;
import becca.gui.MatrixPanel;
import static conceptor.Util.tanh;
import conceptor.chaos.LorenzSystem;
import java.util.LinkedList;
import java.util.List;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.ops.CommonOps;
import static org.ejml.ops.CommonOps.addEquals;
import static org.ejml.ops.CommonOps.identity;
import static org.ejml.ops.CommonOps.insert;
import static org.ejml.ops.CommonOps.invert;
import static org.ejml.ops.CommonOps.mult;
import static org.ejml.ops.CommonOps.scale;
import static org.ejml.ops.CommonOps.subEquals;
import static org.ejml.ops.CommonOps.transpose;


/**
 * Translated from: fig3fig18fig19A_DemoApertureChaos_Main.py
 */
public class ConceptorChaos {

    //#%%% Experiment control
    double randState = 1;
    double newNets = 1;
    boolean newSystemScalings = true;
    double newChaosData = 1;
    
    //#%%% Setting system params
    int Netsize = 20;
    
    //#% network size
    double NetSR = 1.5;
    
    //#% spectral radius
    double NetinpScaling = 1.5;
    //#% scaling of pattern feeding weights
    double BiasScaling = 0.2;
    //#% size of bias
    //#%%% Weight learning
    double TychonovAlphaEqui = .0001;
    //#% regularizer for equi weight training
    int learnLength = 1000;
    
    double alpha = 4;
    
    double noiseMix = 0.05;
    
    //#% for learning W and output weights
    double TychonovAlphaReadout = 0.01;
    
    
    private double Netconnectivity;
    private final DenseMatrix64F WinRaw;
    private final DenseMatrix64F WstarRaw;
    private final DenseMatrix64F WbiasRaw;
    private DenseMatrix64F Wstar;
    private DenseMatrix64F Win;
    private DenseMatrix64F Wbias;
    private boolean collectReservoir = true;
    
    int inputSize = 2;
    private DenseMatrix64F W;
    
    List<DenseMatrix64F> xHistory = new LinkedList();
    List<DenseMatrix64F> xOldHistory = new LinkedList();
    List<DenseMatrix64F> inputHistory = new LinkedList();
    private DenseMatrix64F Wout;
    
    //#% Initializations for random numbers
    //plt.randn('state', randstate)
    //np.random.rand('twister', randstate)
    
    public ConceptorChaos() {
        //#%%%% Demo of aperture sweep
        
        JPanel p = new JPanel();
        p.setLayout(new BoxLayout(p, BoxLayout.PAGE_AXIS));
        
        
        //#% Create raw weights
        if (Netsize <= 20) 
            Netconnectivity = 1.0;
        else
            Netconnectivity = 10.0/((double)Netsize);

        WinRaw = new DenseMatrix64F(Netsize, 2);
        becca.core.Util.matrixRandomGaussian(WinRaw, 1.0, 0);
        
        WstarRaw = Util.newInternalWeights(Netsize, Netconnectivity);

        WbiasRaw = new DenseMatrix64F(Netsize, 1);
        becca.core.Util.matrixRandomGaussian(WbiasRaw, 1.0, 0);
         

        //#% Scale raw weights and initialize weights
        if (newSystemScalings) {
            Wstar = WstarRaw.copy(); scale(NetSR, Wstar);
            Win = WinRaw.copy(); scale(NetinpScaling, Win);
            Wbias = WbiasRaw.copy(); scale(BiasScaling, Wbias);
        }
        
        MatrixPanel mp = new MatrixPanel();
        mp.setLayout(new BoxLayout(mp, BoxLayout.PAGE_AXIS));
        mp.addMatrix("Wstar", Wstar);
        mp.addMatrix("Win", Win);
        mp.addMatrix("Wbias", Wbias);
        
        p.add(mp);
        
        
        //#% % learn equi weights
        //#% harvest data from network externally driven by patterns
        //allTrainArgs = np.zeros(Netsize, np.dot(Np, learnLength))
        //allTrainOldArgs = np.zeros(Netsize, np.dot(Np, learnLength))
        //allTrainOuts = np.zeros(1., np.dot(Np, learnLength))

        
        //patternCollectors = cell(1., Np)
        //xCollectorsCentered = cell(1., Np)
        //xCollectors = cell(1., Np)
        //patternRs = cell(1., Np)
        //startXs = np.zeros(Netsize, Np)
        //#% collect data from driving native reservoir with different drivers
        //for p in np.arange(1., (Np)+1):
        //    patt = patts.cell[int(p)-1]
        //    #% current pattern generator
        DenseMatrix64F x = new DenseMatrix64F(Netsize, 1);
        
        DynamicPlot plot1 = new DynamicPlot("in0", "in1");
        DynamicPlot plot2 = new DynamicPlot("out0", "out1");
        //DynamicPlot plot3 = new DynamicPlot("NRMSE");
        p.add(plot1);
        p.add(plot2);
        //p.add(plot3);
        
        for (int t = 0; t < learnLength; t++) {            

            //u = the current input pattern
            DenseMatrix64F u = DenseMatrix64F.wrap(inputSize, 1, nextInput(t));
            DenseMatrix64F uNoise = new DenseMatrix64F(inputSize, 1);
            matrixRandomGaussian(uNoise, 0.5, 0);
            
            scale(1.0 - noiseMix, u);
            scale(noiseMix, uNoise);
            addEquals(u, uNoise);
            
            plot1.update(t, u.get(0, 0), u.get(1, 0));
            
            
            
            DenseMatrix64F xOld = x;
            
            //x = tanh(Wstar * x + Win * u + Wbias);            
            //x = np.tanh((np.dot(Wstar, x)+np.dot(Win, u)+Wbias))
            
            DenseMatrix64F WstarX = new DenseMatrix64F(Netsize, 1);
            mult(Wstar, x, WstarX);
            
            DenseMatrix64F WinU = new DenseMatrix64F(Netsize, 1);
            mult(Win, u, WinU);
            
            x = WstarX;
            addEquals(x, WinU);
            addEquals(x, Wbias);            
            Util.tanh(x);
            
           
            
            if (collectReservoir) {
                xHistory.add(x);
                xOldHistory.add(xOld);
                inputHistory.add(u);
                
                //xCollector[:,int((n-washoutLength))-1] = x
                //xOldCollector[:,int((n-washoutLength))-1] = xOld
                //if p == 2. or p == 3.:
                //#% the Lorenz and MG observers are rescaled to [0 1]
                //#% the other two are already in that range
                //  pCollector[0,int((n-washoutLength))-1] = np.dot(0.5, u[0]+1.)
                //else:
                //  pCollector[0,int((n-washoutLength))-1] = u[0]
                
            }
        }
        
        
        DenseMatrix64F C = calculateConceptor();

        
        DenseMatrix64F generated = generate(C, 0, 1000);
        
        for (int t = 0; t < generated.numCols; t++) {
            plot2.update(t, generated.get(0, t), generated.get(1, t));
        }
        

        
        
        //FINISHED
        
        JFrame f = new JFrame();
        f.getContentPane().add(new JScrollPane(p));
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setSize(800, 1000);
        f.setVisible(true);
    
    }
    
    protected DenseMatrix64F getHistoryMatrix(List<DenseMatrix64F> history) {
        int historyLength = Math.min(history.size(), learnLength);
        DenseMatrix64F r = new DenseMatrix64F(history.get(0).numRows, historyLength);
        
        int n = 0;
        
        for (DenseMatrix64F h : history.subList(Math.max(0, history.size()-historyLength), history.size()-1)) {

            CommonOps.insert(h, r, 0, n++);
        }
        return r;
    }
    
    public DenseMatrix64F colMean(DenseMatrix64F m) {
        DenseMatrix64F result = new DenseMatrix64F(m.numRows, 1);
        for (int i = 0; i < m.numRows; i++) {
            double mean = 0;
            for (int c = 0; c < m.numCols; c++) {
                mean += m.get(i, c);
            }
            mean /= ((double)m.numCols);
            result.set(i, 0, mean);
        }
        return result;
    }
    
    public DenseMatrix64F generate(DenseMatrix64F conceptor, int washoutTime, int generateTime) {
        DenseMatrix64F output = new DenseMatrix64F(inputSize, generateTime);

        DenseMatrix64F x = new DenseMatrix64F(Netsize, 1);       
        matrixRandomUniform(x, 1.0, 0);
        
        
        
        int outputNum = 0;
        
        
        for (int t = 0; t < washoutTime+generateTime; t++) {
            //x = Cmix * tanh(W * x + Wbias);
            DenseMatrix64F Wx = new DenseMatrix64F(W.numRows, x.numCols);
            mult(W, x, Wx);
            addEquals(Wx, Wbias);
            tanh(Wx);
            mult(conceptor, Wx, x);
            
            if (t > washoutTime) {
                DenseMatrix64F o = new DenseMatrix64F(Wout.numRows, x.numCols);
                mult(Wout, x, o);
                
                insert(o, output, 0, outputNum++);
            }
        }
        
        
        return output;
        
    }
    
    /**
     * calculates conceptor matrix for the current reservoir
     * @return 
     */
    public DenseMatrix64F calculateConceptor() {
        DenseMatrix64F xCollector = getHistoryMatrix(xHistory);
        DenseMatrix64F xOldCollector = getHistoryMatrix(xOldHistory);
        
        //mean(A,2) is a column vector containing the mean of each row
        {
        /*
            DenseMatrix64F xCollectorMean =  colMean(xCollector);        
            xCollectorMean = broadcastRows(xCollectorMean, learnLength);

            DenseMatrix64F xCollectorCentered = xCollector.copy();
            CommonOps.subEquals(xCollectorCentered, xCollectorMean);
        */
        }
        
        DenseMatrix64F pCollector = getHistoryMatrix(inputHistory);
        
        DenseMatrix64F allTrainArgs = xCollector;
        DenseMatrix64F allTrainOldArgs = xOldCollector;
        DenseMatrix64F allTrainOuts = pCollector;
        
        
        DenseMatrix64F R, R2;
        R = new DenseMatrix64F(xCollector.numRows, xCollector.numRows);
        mult(xCollector, transpose(xCollector,null), R);
        scale(1.0/((double)learnLength), R);
        
    
        
        //[Ux, Sx, Vs] = svd(R)
        //[U,S,V] = svd(X) produces a diagonal matrix S of the same dimension as X, with nonnegative diagonal elements in decreasing order, and unitary matrices U and V so that X = U*S*V'.



        /*
        In mathematics, the conjugate transpose, Hermitian transpose, Hermitian conjugate, bedaggered matrix, or adjoint matrix of an m-by-n matrix A with complex entries is the n-by-m matrix A* obtained from A by taking the transpose and then taking the complex conjugate of each entry (i.e., negating their imaginary parts but not their real parts).            
        */                        
        //Wout = np.dot(
        //          np.dot(
        //              linalg.inv( (np.dot(allTrainArgs, allTrainArgs.conj().T)+
        //          np.dot(TychonovAlphaReadout, np.eye(Netsize)))), allTrainArgs),

        //          allTrainOuts.conj().T
        //       ).conj().T        
        /*
          Wout = (
                   inv(allTrainArgs * allTrainArgs' + TychonovAlphaReadout * eye(Netsize))
                    * allTrainArgs * allTrainOuts'
                 )';
        */
        DenseMatrix64F allTrainArgsSq = new DenseMatrix64F(Netsize, Netsize);
        mult(allTrainArgs, transpose(allTrainArgs, null), allTrainArgsSq);


        DenseMatrix64F tychonovAlphaEye = identity(Netsize);
        scale(TychonovAlphaReadout, tychonovAlphaEye);

        DenseMatrix64F r1 = allTrainArgsSq;
        addEquals(r1, tychonovAlphaEye);            
        invert(r1);

        DenseMatrix64F allTrainArgOut = new DenseMatrix64F(r1.numRows, allTrainArgs.numCols);
        mult(r1, allTrainArgs, allTrainArgOut);

        Wout = new DenseMatrix64F(Netsize, allTrainOuts.numRows);

        mult(allTrainArgOut, transpose(allTrainOuts,null), Wout);

        Wout = transpose(Wout, null);
            
            
        //#%%% compute W
        //Wtargets = atanh(allTrainArgs)-matcompat.repmat(Wbias, 1., np.dot(Np, learnLength))
        //Wtargets = (atanh(allTrainArgs) - repmat(Wbias,1,Np*learnLength));
        //Wtargets = (atanh(x) - Wbias)
        DenseMatrix64F Wtargets = allTrainArgs.copy();
        Util.atanh(Wtargets);
        subEquals(Wtargets, broadcastRows(Wbias,Wtargets.numCols));


        //#% training error
        /*
        The normalized root-mean-square deviation or error (NRMSD or NRMSE) is the RMSD divided by the range of observed values of a variable being predicted
        */
        //NRMSE_readout = nrmse(np.dot(Wout, allTrainArgs), allTrainOuts)
        //NRMSE_readout = nrmse(Wout*allTrainArgs, allTrainOuts);
        //NRMSE_readout = nrmse(Wout*x, input);
        //np.disp(sprintf('NRMSE readout: %g', NRMSE_readout))
        DenseMatrix64F Woutx = matrixVector(Wout, allTrainArgs);
        double mse = 0;
        double minInput = CommonOps.elementMin(allTrainOuts);
        double maxInput = CommonOps.elementMax(allTrainOuts);
        for (int i = 0; i < Woutx.numRows; i++) {
            for (int c = 0; c < Woutx.numCols; c++) {                
                double delta = Woutx.get(i, c) - allTrainOuts.get(i, 0);
                mse += delta*delta;
            }
        }
        mse/=((double)Woutx.elements);
        double rmse = Math.sqrt(mse);
        double nrmse = (minInput!=maxInput) ? rmse / (maxInput-minInput) : rmse;

        //W = np.dot(np.dot(linalg.inv((np.dot(allTrainOldArgs, allTrainOldArgs.conj().T)+np.dot(TychonovAlphaEqui, np.eye(Netsize)))), allTrainOldArgs), Wtargets.conj().T).conj().T
        //W = (inv(allTrainOldArgs * allTrainOldArgs' + TychonovAlphaEqui * eye(Netsize)) * allTrainOldArgs * Wtargets')';
        DenseMatrix64F allTrainOldArgsSq = new DenseMatrix64F(Netsize, Netsize);
        mult(allTrainOldArgs, transpose(allTrainOldArgs,null), allTrainOldArgsSq);

        DenseMatrix64F tychonovAlphaEquiEye = identity(Netsize);
        scale(TychonovAlphaEqui, tychonovAlphaEquiEye);

        DenseMatrix64F r2 = allTrainOldArgsSq;
        addEquals(r2, tychonovAlphaEquiEye);

        invert(r2);

        DenseMatrix64F r2factor = new DenseMatrix64F(r2.numRows, allTrainOldArgs.numCols);
        mult(r2, allTrainOldArgs, r2factor);

        W = new DenseMatrix64F(Netsize, Netsize);
        mult(r2factor, transpose(Wtargets,null), W);
        W = transpose(W, null);

     
            
        //% % compute projectors
        SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(R.numRows, R.numCols, true, true, false);
        svd.decompose(R);
        //used later when calculating Conceptors
        
        

        
        
        DenseMatrix64F U = svd.getU(null, false);
        DenseMatrix64F S = svd.getW(null);
        DenseMatrix64F V = svd.getV(null, false);
        //Snew = (S * inv(S + alpha^(-2) * eye(Netsize)));
        //C = U * Snew * U';
        
        DenseMatrix64F Sfactor = S.copy();
        DenseMatrix64F Seye = CommonOps.identity(Netsize);
        //scale(Math.pow(alpha, -2), Seye);         //fig3 does not scale Seye
        addEquals(Sfactor, Seye);
        invert(Sfactor);
        
        DenseMatrix64F Snew = new DenseMatrix64F(S.numRows, Sfactor.numCols);
        mult(S, Sfactor, Snew);
        
        DenseMatrix64F Cfactor = new DenseMatrix64F(U.numRows, Snew.numCols);
        mult(U, Snew, Cfactor);
        
        DenseMatrix64F C = new DenseMatrix64F(Cfactor.numRows, U.numRows);
        mult(Cfactor, transpose(U,null), C);
        
        return C;        
    }
    
    
    public static double sqr(double x) { return x*x; }
    
    public double[] nextInput(double t) {
             /*return new double[] { 
                 0.5 + Math.sin(t/5.0)*0.5, 
                 0.5 + Math.cos(sqr(t/2.0)/100.0)*0.5 };       */
             return new double[] { 
                 0.5 + Math.sin(t/15.0)*0.5, 
                 0.5}; //0.5 + Math.cos(t/2.5)*0.5 };       

    }

    public void initChaos() {
        {
            int L = 1000;
            int points = 200;
        
            /*
            ls = [10.036677794959058; 9.98674414052542; 
    29.024692318601613] + 0.01 * randn(3,1);

delta = 1 / incrementsperUnit; % length of discrete approximation update interval
                    */
        
            double sigma = 10.0, b = 4.0/3, r = 28.0; 
            double[] x0 = new double[] { 
               10.036677794959058, 
               9.98674414052542, 
               29.024692318601613 };
            
            DenseMatrix64F LorenzSequence =  new DenseMatrix64F(
                    LorenzSystem.generate(points, 0.05, sigma, b, r, x0));

            Util.normalizeColumn(LorenzSequence, 0);
            Util.normalizeColumn(LorenzSequence, 1);
            Util.normalizeColumn(LorenzSequence, 2);
            
            //mp.add(new Plot2D(LorenzSequence, 0, 2, 200, 200));
            //mp.add(new Plot2D(LorenzSequence, 0, 1, 200, 200));
            
        }
        
        {
            // OTHER CHAOS FUNCTIONS
            //#% Set pattern handles
            //if newChaosData:
            //    patts = cell(1., 4.)
            //    L = washoutLength+learnLength
            //    LorenzSeq = generateLorenzSequence2D(200., 15., L, 5000.)
            //    patts.cell[1] = lambda n: 2.*LorenzSeq[:,int(n)-1]-1.
            //    RoesslerSeq = generateRoesslerSequence2D(200., 150., L, 5000.)
            //    patts.cell[0] = lambda n: RoesslerSeq[:,int(n)-1]
            //    MGSeq = generateMGSequence2D(17., 10., 3., L, 5000.)
            //    patts.cell[2] = lambda n: 2.*MGSeq[:,int(n)-1]-1.
            //    HenonSeq = generateHenonSequence2D(L, 1000.)
            //    patts.cell[3] = lambda n: HenonSeq[:,int(n)-1]
            //    Np = 4.
            //        

        }
        
    }
    

//#%%
//plt.figure(10.)
//plt.clf
//for p in np.arange(1., 5.0):
//    plt.subplot(2., 2., p)
//    [U, S, V] = plt.svd(patternRs.cell[int(p)-1])
//    plt.plot(np.log10(np.diag(S)))
//    
//#%% Plotting sweeps through aperture
//#%bestAlphas = [2000 80 350 300];
//#%bestAlphas = [10^3 10^2.4 10^2.7 10^2.8];
//bestAlphas = np.array(np.hstack((10.**3., 10.**2.6, 10.**3.1, 10.**2.8)))
//factorsShort = np.array(np.hstack((7., 6.012, 7., 7.)))
//halfPlotNumberShort = 2.
//exponentsShort = np.arange(-halfPlotNumberShort, (halfPlotNumberShort)+1)
//NalphasShort = 2.*halfPlotNumberShort+1.
//allAlphasShort = np.zeros(4., NalphasShort)
//for i in np.arange(1., (NalphasShort)+1):
//    allAlphasShort[:,int(i)-1] = bestAlphas.conj().T*factorsShort.conj().T**exponentsShort[int(i)-1]
//    
//sigmaPL = np.zeros(NalphasShort, Netsize, 4.)
//#%%
//testLengthes = 900.*np.array(np.hstack((1., 1., 1., 1.)))
//plotLengthesDelayEmbed = 500.*np.array(np.hstack((1., 1., 1., 1.)))
//delays = np.array(np.hstack((2., 2., 3., 1.)))
//#%%
//plt.figure(7.)
//plt.clf
//set(plt.gcf, 'WindowStyle', 'normal')
//set(plt.gcf, 'Position', np.array(np.hstack((600., 500., 1200., 200.))))
//plotInd = 0.
//for p in np.array(np.hstack((1., 3., 4.))):
//    plotInd = plotInd+1.
//    C = Cs.cell[0,int(p)-1]
//    testLength = testLengthes[int(p)-1]
//    plotLengthDelayEmbed = plotLengthesDelayEmbed[int(p)-1]
//    delay = delays[int(p)-1]
//    alpha = allAlphasShort[int(p)-1,2]
//    apSweepPL = np.zeros(1., testLength)
//    Calpha = PHI(C, alpha)
//    [U, S, V] = plt.svd(Calpha)
//    sigmaPL[int(i)-1,:,int(p)-1] = np.diag(S).conj().T
//    x = startXs[:,int(p)-1]
//    for n in np.arange(1., (testLength)+1):
//        z = np.tanh((np.dot(W, x)+Wbias))
//        x = np.dot(Calpha, z)
//        apSweepPL[0,int(n)-1] = np.dot(Wout, x)
//        
//    plt.subplot('Position', np.array(np.hstack(((plotInd-1.)*1./3., 0., 1./6., 1.))))
//    if p == 4.:
//        plt.plot(apSweepPL[0,0:plotLengthDelayEmbed], apSweepPL[0,int(1.+delay)-1:plotLengthDelayEmbed+delay], 'b.')
//    else:
//        plt.plot(apSweepPL[0,0:plotLengthDelayEmbed], apSweepPL[0,int(1.+delay)-1:plotLengthDelayEmbed+delay], 'b', 'MarkerSize', 1.)
//        
//    
//    rectangle('Position', np.array(np.hstack((0.01, 0.85, 0.45, 0.14))), 'FaceColor', 'w', 'EdgeColor', 'none')
//    plt.text(0.05, 0.93, num2str(alpha, 2.), 'FontSize', 18., 'FontWeight', 'bold')
//    set(plt.gca, 'YLim', np.array(np.hstack((0., 1.))), 'XLim', np.array(np.hstack((0., 1.))), 'xtick', np.array([]), 'ytick', np.array([]))
//    plt.subplot('Position', np.array(np.hstack(((plotInd-1.)*1./3.+1./6., 0., 1./6., 1.))))
//    if p == 4.:
//        plt.plot(patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed)-1:0-delay](), patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed+delay)-1:](), 'g.')
//    else:
//        plt.plot(patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed)-1:0-delay](), patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed+delay)-1:](), 'g', 'MarkerSize', 1.)
//        
//    
//    set(plt.gca, 'YLim', np.array(np.hstack((0., 1.))), 'XLim', np.array(np.hstack((0., 1.))), 'xtick', np.array([]), 'ytick', np.array([]))
//    
//#%%
//for p in np.arange(1., 5.0):
//    C = Cs.cell[0,int(p)-1]
//    testLength = testLengthes[int(p)-1]
//    plotLengthDelayEmbed = plotLengthesDelayEmbed[int(p)-1]
//    delay = delays[int(p)-1]
//    alphas = allAlphasShort[int(p)-1,:]
//    apSweepPL = np.zeros(NalphasShort, testLength)
//    quotaPL = np.zeros(1., NalphasShort)
//    for i in np.arange(1., (NalphasShort)+1):
//        alpha = alphas[int(i)-1]
//        Calpha = PHI(C, alpha)
//        [U, S, V] = plt.svd(Calpha)
//        sigmaPL[int(i)-1,:,int(p)-1] = np.diag(S).conj().T
//        quotaPL[0,int(i)-1] = matdiv(np.trace(Calpha), Netsize)
//        x = startXs[:,int(p)-1]
//        for n in np.arange(1., (testLength)+1):
//            z = np.tanh((np.dot(W, x)+Wbias))
//            x = np.dot(Calpha, z)
//            apSweepPL[int(i)-1,int(n)-1] = np.dot(Wout, x)
//            
//        
//    plt.figure(p)
//    plt.clf
//    set(plt.gcf, 'WindowStyle', 'normal')
//    set(plt.gcf, 'Position', np.array(np.hstack((900.+(p-1.)*50., 200.+(p-1.)*100., 600., 400.))))
//    for i in np.arange(1., (NalphasShort)+1):
//        plt.subplot('Position', np.array(np.hstack((np.mod((i-1.), 3.)/3., (-np.ceil((i/3.))+2.)/2., 1./3., 1./2.))))
//        if p == 4.:
//            plt.plot(apSweepPL[int(i)-1,0:plotLengthDelayEmbed], apSweepPL[int(i)-1,int(1.+delay)-1:plotLengthDelayEmbed+delay], '.')
//        else:
//            plt.plot(apSweepPL[int(i)-1,0:plotLengthDelayEmbed], apSweepPL[int(i)-1,int(1.+delay)-1:plotLengthDelayEmbed+delay], 'MarkerSize', 1.)
//            
//        
//        rectangle('Position', np.array(np.hstack((0.01, 0.9, 0.48, 0.08))), 'FaceColor', 'w', 'EdgeColor', 'none')
//        rectangle('Position', np.array(np.hstack((0.01, 0.8, 0.28, 0.1))), 'FaceColor', 'w', 'EdgeColor', 'none')
//        plt.text(0.05, 0.95, num2str(alphas[int(i)-1], 2.), 'FontSize', 18., 'FontWeight', 'bold')
//        plt.text(0.05, 0.85, num2str(quotaPL[0,int(i)-1], 2.), 'FontSize', 18., 'FontWeight', 'bold')
//        set(plt.gca, 'YLim', np.array(np.hstack((0., 1.))), 'XLim', np.array(np.hstack((0., 1.))), 'xtick', np.array([]), 'ytick', np.array([]))
//        
//    plt.subplot('Position', np.array(np.hstack((np.mod((5.0.), 3.)/3., (-np.ceil((6./3.))+2.)/2., 1./3., 1./2.))))
//    if p == 4.:
//        plt.plot(patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed)-1:0-delay](), patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed+delay)-1:](), 'g.')
//    else:
//        plt.plot(patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed)-1:0-delay](), patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed+delay)-1:](), 'g', 'MarkerSize', 1.)
//        
//    
//    set(plt.gca, 'YLim', np.array(np.hstack((0., 1.))), 'XLim', np.array(np.hstack((0., 1.))), 'xtick', np.array([]), 'ytick', np.array([]))
//    
//#%%
//for p in 2.:
//    C = Cs.cell[0,int(p)-1]
//    testLength = testLengthes[int(p)-1]
//    plotLengthDelayEmbed = plotLengthesDelayEmbed[int(p)-1]
//    delay = delays[int(p)-1]
//    alphas = allAlphasShort[int(p)-1,:]
//    apSweepPL = np.zeros(NalphasShort, testLength)
//    quotaPL = np.zeros(1., NalphasShort)
//    for i in np.arange(1., (NalphasShort)+1):
//        alpha = alphas[int(i)-1]
//        Calpha = PHI(C, alpha)
//        [U, S, V] = plt.svd(Calpha)
//        sigmaPL[int(i)-1,:,int(p)-1] = np.diag(S).conj().T
//        quotaPL[0,int(i)-1] = matdiv(np.trace(Calpha), Netsize)
//        x = startXs[:,int(p)-1]
//        for n in np.arange(1., (testLength)+1):
//            z = np.tanh((np.dot(W, x)+Wbias))
//            x = np.dot(Calpha, z)
//            apSweepPL[int(i)-1,int(n)-1] = np.dot(Wout, x)
//            
//        
//    plt.figure(5.)
//    plt.clf
//    set(plt.gcf, 'WindowStyle', 'normal')
//    set(plt.gcf, 'Position', np.array(np.hstack((600., 300., 1200., 200.))))
//    for i in np.arange(1., 6.0):
//        plt.subplot('Position', np.array(np.hstack((np.mod((i-1.), 6.)/6., 0., 1./6., 1.))))
//        if p == 4.:
//            plt.plot(apSweepPL[int(i)-1,0:plotLengthDelayEmbed], apSweepPL[int(i)-1,int(1.+delay)-1:plotLengthDelayEmbed+delay], 'k.')
//        else:
//            plt.plot(apSweepPL[int(i)-1,0:plotLengthDelayEmbed], apSweepPL[int(i)-1,int(1.+delay)-1:plotLengthDelayEmbed+delay], 'b', 'MarkerSize', 1.)
//            
//        
//        plt.text(0.05, 0.93, num2str(alphas[int(i)-1], 2.), 'FontSize', 18., 'FontWeight', 'bold')
//        set(plt.gca, 'YLim', np.array(np.hstack((0., 1.))), 'XLim', np.array(np.hstack((0., 1.))), 'xtick', np.array([]), 'ytick', np.array([]))
//        
//    plt.subplot('Position', np.array(np.hstack((5./6., 0., 1./6., 1.))))
//    if p == 4.:
//        plt.plot(patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed)-1:0-delay](), patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed+delay)-1:](), 'k.')
//    else:
//        plt.plot(patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed)-1:0-delay](), patternCollectors.cell[0,int(p)-1,0,int(0-plotLengthDelayEmbed+delay)-1:](), 'g', 'MarkerSize', 1.)
//        
//    
//    set(plt.gca, 'YLim', np.array(np.hstack((0., 1.))), 'XLim', np.array(np.hstack((0., 1.))), 'xtick', np.array([]), 'ytick', np.array([]))
//    
//#%%
//#%% Plotting attenuations
//testLengthesAtt = 200.*np.array(np.hstack((1., 1., 1., 1.)))
//factors = np.dot(np.dot(10.**(2./8.), 0.8), np.array(np.hstack((1., 1., 1., 1.))))
//halfPlotNumber = 10.
//exponents = np.arange(-halfPlotNumber, (halfPlotNumber)+1)
//Nalphas = 2.*halfPlotNumber+1.
//allAlphas = np.zeros(4., Nalphas)
//allQuotas = np.zeros(4., Nalphas)
//allNorms = np.zeros(4., Nalphas)
//attenuationPL = np.zeros(4., Nalphas)
//diffPL = np.zeros(4., Nalphas)
//zsPL = np.zeros(4., Nalphas)
//for i in np.arange(1., (Nalphas)+1):
//    allAlphas[:,int(i)-1] = bestAlphas.conj().T*factors.conj().T**exponents[int(i)-1]
//    
//for p in np.arange(1., 5.0):
//    C = Cs.cell[0,int(p)-1]
//    testLength = testLengthesAtt[int(p)-1]
//    delay = delays[int(p)-1]
//    alphas = allAlphas[int(p)-1,:]
//    Nalphas = matcompat.size(alphas, 2.)
//    sigmaPL = np.zeros(Nalphas, Netsize)
//    for i in np.arange(1., (Nalphas)+1):
//        alpha = alphas[int(i)-1]
//        Calpha = PHI(C, alpha)
//        allQuotas[int(p)-1,int(i)-1] = matdiv(np.trace(Calpha), Netsize)
//        allNorms[int(p)-1,int(i)-1] = linalg.norm(Calpha, 'fro')**2.
//        [U, S, V] = plt.svd(Calpha)
//        sigmaPL[int(i)-1,:] = np.diag(S).conj().T
//        x = startXs[:,int(p)-1]
//        att = 0.
//        diff = 0.
//        zs = 0.
//        for n in np.arange(1., (testLength)+1):
//            z = np.tanh((np.dot(W, x)+Wbias))
//            x = np.dot(Calpha, z)
//            att = att+matdiv(linalg.norm((x-z))**2., linalg.norm(z)**2.)
//            zs = zs+linalg.norm(z)**2.
//            diff = diff+linalg.norm((x-z))**2.
//            
//        attenuationPL[int(p)-1,int(i)-1] = matdiv(att, testLength)
//        diffPL[int(p)-1,int(i)-1] = matdiv(diff, testLength)
//        zsPL[int(p)-1,int(i)-1] = matdiv(zs, testLength)
//        
//    
//#%%
//plt.figure(6.)
//plt.clf
//set(plt.gcf, 'WindowStyle', 'normal')
//set(plt.gcf, 'Position', np.array(np.hstack((700., 200., 1000., 180.))))
//fs = 18.
//yValuesForDots = np.dot(-4., np.array(np.hstack((1., 1., 1., 1.))))
//plotInd = 0.
//for p in np.array(np.hstack((2., 1., 3., 4.))):
//    plotInd = plotInd+1.
//    plt.subplot(1., 4., plotInd)
//    plt.hold(on)
//    plt.plot(np.log10(allAlphas[int(p)-1,:]), np.log10((diffPL[int(p)-1,:]/zsPL[int(p)-1,:])), 'k', 'LineWidth', 2.)
//    if p == 2.:
//        plt.plot(np.log10(allAlphasShort[int(p)-1,:]), np.dot(yValuesForDots[int(p)-1], np.ones(1., 5.)), 'b.', 'MarkerSize', 35.)
//    else:
//        plt.plot(np.log10(allAlphasShort[int(p)-1,2]), yValuesForDots[int(p)-1], 'b.', 'MarkerSize', 35.)
//        
//    
//    plt.hold(off)
//    set(plt.gca, 'Box', 'on', 'FontSize', fs, 'XLim', np.array(np.hstack((1., 5.))), 'XTick', np.arange(1., 6.0), 'YLim', np.array(np.hstack((-8., 0.))))
//    if p == 1.:
//        plt.title('Roessler')
//    elif p == 2.:
//        plt.title('Lorenz')
//        plt.xlabel('log10 aperture')
//        
//    elif p == 3.:
//        plt.title('Mackey-Glass')
//        #%ylabel('log10 attenuation');
//        
//    else:
//        plt.title('H\E9non')
//        
//    
//    
//#%%
//plt.figure(16.)
//plt.clf
//set(plt.gcf, 'WindowStyle', 'normal')
//set(plt.gcf, 'Position', np.array(np.hstack((700., 200., 1000., 600.))))
//fs = 22.
//yValuesForDots = np.dot(-4., np.array(np.hstack((1., 1., 1., 1.))))
//plotInd = 0.
//for p in np.array(np.hstack((2., 1., 3., 4.))):
//    plotInd = plotInd+1.
//    plt.subplot(2., 2., plotInd)
//    plt.hold(on)
//    plt.plot(np.log10(allAlphas[int(p)-1,:]), np.log10((diffPL[int(p)-1,:]/zsPL[int(p)-1,:])), 'k', 'LineWidth', 2.)
//    plt.plot(np.log10(allAlphasShort[int(p)-1,:]), np.dot(yValuesForDots[int(p)-1], np.ones(1., 5.)), 'b.', 'MarkerSize', 35.)
//    plt.hold(off)
//    set(plt.gca, 'Box', 'on', 'FontSize', fs, 'XLim', np.array(np.hstack((1., 5.))), 'XTick', np.arange(1., 6.0), 'YLim', np.array(np.hstack((-8., 0.))))
//    if p == 1.:
//        plt.title('Roessler')
//    elif p == 2.:
//        plt.title('Lorenz')
//        
//    elif p == 3.:
//        plt.title('Mackey-Glass')
//        plt.xlabel('log10 aperture')
//        plt.ylabel('log10 attenuation')
//        
//    else:
//        plt.title('H\E9non')
//        
//    
    
    public static void main(String[] args) {
        new ConceptorChaos();
    }
   
    
}
