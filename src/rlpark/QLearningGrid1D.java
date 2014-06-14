/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rlpark;

/**
 *
 * @author me
 */

import becca.test.World;
import becca.world.Grid1DSimple;
import java.util.Arrays;
import java.util.Random;

import rlpark.plugin.rltoys.agents.rl.LearnerAgentFA;
import rlpark.plugin.rltoys.algorithms.control.ControlLearner;
import rlpark.plugin.rltoys.algorithms.control.acting.EpsilonGreedy;
import rlpark.plugin.rltoys.algorithms.control.qlearning.QLearning;
import rlpark.plugin.rltoys.algorithms.control.qlearning.QLearningControl;
import rlpark.plugin.rltoys.algorithms.functions.stateactions.TabularAction;
import rlpark.plugin.rltoys.algorithms.functions.states.Projector;
import rlpark.plugin.rltoys.algorithms.traces.RTraces;
import rlpark.plugin.rltoys.envio.actions.Action;
import rlpark.plugin.rltoys.envio.observations.Legend;
import rlpark.plugin.rltoys.envio.policy.Policy;
import rlpark.plugin.rltoys.envio.rl.TRStep;
import rlpark.plugin.rltoys.experiments.runners.Runner;
import rlpark.plugin.rltoys.experiments.runners.AbstractRunner.RunnerEvent;
import rlpark.plugin.rltoys.math.vector.RealVector;
import rlpark.plugin.rltoys.math.vector.implementations.PVector;
import rlpark.plugin.rltoys.math.vector.implementations.SVector;
import rlpark.plugin.rltoys.problems.ProblemDiscreteAction;
import zephyr.plugin.core.api.Zephyr;
import zephyr.plugin.core.api.monitoring.annotations.Monitor;
import zephyr.plugin.core.api.signals.Listener;
import zephyr.plugin.core.api.synchronization.Clock;



@Monitor
public class QLearningGrid1D implements Runnable {

    
    public class WorldProblem implements ProblemDiscreteAction {
        private final World world;
        private final Action[] actions;
        private TRStep step;
        double[] sensor;
        private double reward, rewardTotal;

        private Projector getMarkovProjector() {
            final int size = sensor.length;
            
            final PVector p = new PVector(sensor);
            final SVector projection = new SVector(size);
            return new Projector() {
              @Override
              public int vectorSize() {
                return size;
              }

              @Override
              public double vectorNorm() {
                return 1;
              }

              @Override
              public RealVector project(double[] obs) {
                  projection.clear();
                  projection.add(p);
                  return projection;
              }
            };
            
        }

        private double getReward() {
            double r = rewardTotal;
            rewardTotal = 0;
            return r;
        }

        public class WorldAction implements Action {
            private final int ii;
            public WorldAction(int i) {
                this.ii = i;
            }
            public String toString() {return "action" +  ii; }
            public int i() { return ii; }            
        }
        
        
        private WorldProblem(World w) {
            this.world = w;
            
            actions = new Action[w.getNumActions()];
            for (int i = 0; i < w.getNumActions(); i++) {
                final int ii = i;
                actions[i] = new WorldAction(i);                
            }
            sensor = new double[w.getNumSensors()];
        }

        @Override
        public Action[] actions() {
            return actions;
        }

        @Override
        public TRStep initialize() {
            step = new TRStep(new double[world.getNumActions()], 0);
            return step;
        }

        double[] d = null;
        
        @Override
        public TRStep step(Action action) {
            if (d==null) {
                d = new double[world.getNumActions()];
            }
            else {
                Arrays.fill(d, 0, d.length, 0.0);
            }
            d[((WorldAction)action).i()] = 1.0;
            
            reward = world.step(d, sensor);
            rewardTotal+=reward;
            
            
            step = new TRStep(step, action, sensor, reward);
            return step;                   
        }

        @Override
        public TRStep forceEndEpisode() {
            step = step.createEndingStep();            
            return step;
        }

        @Override
        public TRStep lastStep() {
            return step;
        }

        @Override
        public Legend legend() {
            return null;
        }
            
        
    }
    
    private final WorldProblem problem;
    private final ControlLearner control;
    private final Clock clock = new Clock("QLearningMaze");
    private final Projector projector;
    private final PVector occupancy;
    private final LearnerAgentFA agent;

    public QLearningGrid1D() {
        problem = new WorldProblem(new Grid1DSimple(4, 1000, 0.0, 0.0002));
        
        projector = problem.getMarkovProjector();
        
        occupancy = new PVector(projector.vectorSize());
        
        TabularAction toStateAction = new TabularAction(problem.actions(), projector.vectorNorm(), projector.vectorSize());
        
        double alpha = .15 / projector.vectorNorm();
        double gamma = 1.0;
        double lambda = 0.6;
        QLearning qlearning = new QLearning(problem.actions(), alpha, gamma, lambda, toStateAction, new RTraces());
        double epsilon = 0.3;
        Policy acting = new EpsilonGreedy(new Random(0), problem.actions(), toStateAction, qlearning, epsilon);
        control = new QLearningControl(acting, qlearning);
        agent = new LearnerAgentFA(control, projector);
        //mazeValueFunction = new MazeValueFunction(problem, qlearning, toStateAction, qlearning.greedy());
        Zephyr.advertise(clock, this);
        
        new DynamicChart() {

            @Override
            public double getReward() {
                return problem.getReward();
            }
            
        };
    }

    @Override
    public void run() {
        Runner runner = new Runner(problem, agent);
        runner.onEpisodeEnd.connect(new Listener<RunnerEvent>() {
            @Override
            public void listen(RunnerEvent eventInfo) {
                System.out.println(String.format("Episode %d: %d steps", eventInfo.nbEpisodeDone, eventInfo.step.time));
            }
        });
        while (clock.tick()) {
            runner.step();
            occupancy.addToSelf(agent.lastState());
        }
    }

    public static void main(String[] args) {
        new QLearningGrid1D().run();
    }
}
