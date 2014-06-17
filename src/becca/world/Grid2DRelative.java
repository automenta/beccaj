package becca.world;

import becca.core.BeccaAgent;
import becca.gui.AgentPanel;
import becca.test.Agent;
import becca.test.Simulation;
import becca.test.World;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import javax.swing.JPanel;

public class Grid2DRelative implements World {

    private final int width, height;

    private final double REWARD_MAGNITUDE;
    private final double JUMP_FRACTION;
    private final double ENERGY_COST_FACTOR;
    private final double MATCH_REWARD_FACTOR;

    private double[] action;

    private final int totalTime;
    private int time;

    private final double noise;
    private int focusPositionX;
    private int focusPositionY;

    private double positionX;
    private double positionY;
    private final Image2DPanel image;
    private final double POSITION_VELOCITY;
    private double[] sensor;

    public class Image2DPanel extends JPanel {
        private final BufferedImage bi;
        private final Graphics2D g2;

        public Image2DPanel() {
            super(new BorderLayout());

            bi
                    = new BufferedImage(width, height,
                            BufferedImage.TYPE_INT_ARGB);

            g2 = bi.createGraphics();
        }

        public void updateImage() {
            g2.clearRect(0, 0, width, height);
            
            int px = 0, py = 0;
            for (int i = 0; i < sensor.length; i++) {
                g2.setPaint(Color.getHSBColor(0.5f, 0.5f, (float)sensor[i]));
                g2.fillRect(px, py, 1, 1);
                px++;
                if (px == width) {
                    px = 0;
                    py++;
                }
            }
            //g2.setPaint(Color.RED);
            //g2.fillRect(focusPositionX-1, focusPositionY-1, 2, 2);

            g2.setPaint(Color.BLUE);
            g2.fillRect((int)Math.round(positionX)-1, (int)Math.round(positionY)-1, 2, 2);
        }
        @Override
        public void paint(Graphics g) {
            g.drawImage(bi, 0, 0, getWidth(), getHeight(), null);
            /*
                        dx, dy, dx+cw, dy+ch,
                        sx, sy, sx+cw, sy+ch,
                        null);*/
            
            validate();
            repaint();
        }        

        

    }

    public Grid2DRelative(int width, int height, int totalTime, double noise, double focusVelocity) {
        this.time = 1;
        this.width = width;
        this.height = height;
        this.ENERGY_COST_FACTOR = 0.01;
        this.MATCH_REWARD_FACTOR = 1.0;
        this.REWARD_MAGNITUDE = 1;
        this.JUMP_FRACTION = 0.002;
        this.POSITION_VELOCITY = 0.1;
        this.noise = noise;

        image = new Image2DPanel();
        AgentPanel.window(image, true);
        
        this.totalTime = totalTime;
        randomFocus();
    }

    protected void randomFocus() {
        this.focusPositionX = (int) (Math.random() * width);
        this.focusPositionY = (int) (Math.random() * height);
    }

    @Override
    public String getName() {
        return getClass().toString();
    }

    @Override
    public int getNumSensors() {
        return height * width;
    }

    @Override
    public int getNumActions() {
        return 4;
    }

    @Override
    public boolean isActive() {
        return time < totalTime;
    }

    double[] action2 = null;

    @Override
    public double step(double[] action, double[] sensor) {

        
        time++;

        this.sensor = sensor;
        this.action = action;

        //# At random intervals, jump to a random position in the world
        if (Math.random() < JUMP_FRACTION) {
            randomFocus();
        }

        /*        
         # Assign basic_feature_input elements as binary. 
         # Represent the presence or absence of the current position in the bin.
         */
        //blur the action
        /*if (action2 == null) action2 = new double[action.length];
         for (int i = 0; i < action2.length; i++) {
         action2[i] = action[i];
         if (i > 0) action2[i] += 0.5 * action[i-1];
         if (i < action2.length-1) action2[i] += 0.5 * action[i+1];
         } */
        /*if (action[0] > 0.5) {
         //nothing
         }*/
        if ((action[0] > 0.5) && !(action[1] > 0.5)) {
            positionX-=POSITION_VELOCITY * action[0];
        }
        if ((action[1] > 0.5) && !(action[0] > 0.5)) {
            positionX+=POSITION_VELOCITY * action[1];
        }
        if (positionX < 0) {
            positionX = 0;
        }
        if (positionX >= width) {
            positionX = width - 1;
        }
        if ((action[2] > 0.5) && !(action[3] > 0.5)) {
            positionY-=POSITION_VELOCITY * action[2];
        }
        if ((action[3] > 0.5) && !(action[2] > 0.5)) {
            positionY+=POSITION_VELOCITY * action[3];
        }
        if (positionY < 0) {
            positionY = 0;
        }
        if (positionY >= height) {
            positionY = height - 1;
        }

        double match = 0;
        double energyCost = 0;

        double dx = (positionX - focusPositionX);
        double dy = (positionY - focusPositionY);
        double dist = Math.sqrt(dx * dx + dy * dy);
        match = 1.0 / (1.0 + dist);

        double reward = REWARD_MAGNITUDE * ((MATCH_REWARD_FACTOR * match) - (energyCost * ENERGY_COST_FACTOR)) - 0.5;

        int ix = (int)Math.round(positionX);
        int iy = (int)Math.round(positionY);
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                final double exp = 1.0; //sharpen

                double fdx = x - focusPositionX;
                double fdy = y - focusPositionY;
                double fdist = Math.sqrt(fdx * fdx + fdy * fdy);

                int i = x + y * width;
                double v = Math.pow(1.0 / (1.0 + fdist/2.0), exp) + (Math.random() * noise);
                if (v < 0.0) {
                    v = 0;
                }

                
                /*if ((ix==x) && (iy==y))
                    v += 0.5;*/
                
                if (v > 1.0) v = 1.0;
                sensor[i] = v;

                
            }

        }

        image.updateImage();
        
        return reward;
    }

    public static void main(String[] args) throws Exception {
        Class<? extends Agent> a = BeccaAgent.class;
        //Class<? extends Agent> a = QLAgent.class;

        new Simulation(a, new Grid2DRelative(12, 10, 11990000, 0.02, 0.01));

    }
}
