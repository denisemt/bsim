package bsim.particle;

import java.awt.*;
import java.util.Vector;

import javax.vecmath.Vector3d;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimUtils;
import bsim.ode.BSimOdeSolver;
import bsim.ode.BSimOdeSystem;

import org.opensourcephysics.display.Dataset;
import org.opensourcephysics.display.DrawingFrame;
import org.opensourcephysics.display.PlottingPanel;
import org.opensourcephysics.numerics.ODE;
import org.opensourcephysics.numerics.ode_interpolation.StateHistory;
import org.opensourcephysics.numerics.ode_solvers.InterpolatorEventSolver;
import org.opensourcephysics.numerics.ode_solvers.rk.*;

import com.kabouterlabs.jodeint.codepack.CodepackLibrary;
import org.bridj.Pointer;

import com.mathworks.toolbox.javabuilder.*;
import stiffODE.*;

/**
 *
 * Class representing a bacterium whose run-tumble motion
 * is affected in a simple way by a single goal chemical.
 */
public class BSimBacterium extends BSimParticle {

	/*
	 * MOVEMENT including CHEMOTAXIS
	 */
	public static enum MotionState {
		/* Quotes from 'Motile behavior of bacteria', Berg */
		/**
		 * "When the motors turn counterclockwise, the filaments rotate in parallel in a bundle that
		 * pushes the cell body steadily forward, and the cell is said to 'run'"
		 */ RUNNING,
		/**
		 * "When the motors turn clockwise, the flagellar filaments work independently, and the cell body
		 * moves erratically with little net displacement; the cell is then said to 'tumble'"
		 */ TUMBLING }
	protected MotionState motionState;

	/**
	 * Magnitude of the flagellar force produced by the cell whilst RUNNING.
	 * Calculated from Stokes law F = 6*PI*radius*viscosity*speed with a radius of 1 micron,
	 * a viscosity of 2.7e-3 Pa s and a speed of 20 microns/s (conditions of
	 * 'Chemotaxis in Escherichia Coli', Berg et al.)
	 */
	protected double forceMagnitude = 1; // pN
	/**
	 * Direction that the cell exerts its flagellar force.
	 */
	protected Vector3d direction;

	/** Bacteria tend to swim towards higher concentrations of this chemical field. */
	protected BSimChemicalField goal;
	/** Memory of previous concentrations of the goal field. */
	protected double[] memory; // molecules/(micron)^3
	/*
	 * 'Temporal comparisons in bacterial chemotaxis', Segall, Berg et al.:
	 * Cells continuously compare the stimulus experienced during past second with
	 * that experienced during the previous 3 seconds and respond to the difference
	 */
	protected double shortTermMemoryDuration = 1; // seconds
	protected double longTermMemoryDuration = 3; // seconds
	/** sim.timesteps(shortTermMemoryDuration) */
	protected double shortTermMemoryLength;
	/** sim.timesteps(longTermMemoryDuration) */
	protected double longTermMemoryLength;
	/** Sensitivity to differences in sequential averages (molecules/(micron)^3). */
	protected double sensitivity = 1;

	/*
	 *
	 * Bottom p86, 'Random Walks in Biology', Berg: (no chemical fields)
	 * "The distribution of run (or tumble) intervals is exponential..
	 * the probability per unit time that a run (or tumble) will end is constant."
	 *
	 * 'Chemotaxis in Escherichia Coli', Berg et al.:
	 * "When a bacterium moves up the gradient the probability per unit time of the
	 * termination of a run decreases; when it moves down the gradient the probability
	 * reverts to the value appropriate to an isotropic concentration of similar concentration."
	 *
	 * The probabilities per unit time for ending a run/tumble could plausibly depend on the value
	 * of the concentration at the particle's location and the derivatives with respect to time and
	 * space. Here we allow only a boolean test for whether the particle is moving
	 * up a chemical gradient in time. This corresponds to the model of Schnizter, Berg et al.,
	 * 'Strategies  for Chemotaxis', p23 but instead of reducing the run termination probability in
	 * the case of an increasing chemical gradient by an amount proportional to the difference
	 * between the sequential averages, we set it to a constant pEndRunUp (distinct from pEndRunElse).
	 *
	 * Values from 'Chemotaxis in Escherichia Coli', Berg et al.
	 */
	/** Probability per per unit time of ending a run when moving up a chemical gradient */
	protected double pEndRunUp = 1/1.07; // 1/seconds
	// mean time to end a run when moving up a gradient = 1/pEndRunUp = 1.07 seconds
	/** Probability per per unit time of ending a run otherwise */
	protected double pEndRunElse = 1/0.86; // 1/seconds
	// mean time to end a run otherwise = 1/pEndRunElse = 0.86 seconds
	/** Probability per per unit time of ending a tumble */
	protected double pEndTumble = 1/0.14; // 1/seconds
	// mean time to end a tumble = 1/pEndTumble = 0.14 seconds
	/** Probability per per unit time of ending a run */
	public double pEndRun() {
		if(goal != null && movingUpGradient()) return pEndRunUp;
		else return pEndRunElse;
	}
	/** Probability per per unit time of ending a tumble. */
	public double pEndTumble() { return pEndTumble; }
	/* Setters */
	public void pEndRunUp(double d) { pEndRunUp = d; }
	public void pEndRunElse(double d) { pEndRunElse = d; }
	public void pEndTumble(double d) { pEndTumble = d; }

	public void setMotionState(MotionState s) { motionState = s; }
	public void setForceMagnitude(double d) { forceMagnitude = d; }
	/**
	 * Set the direction of the cell to the direction of the vector v.
	 */
	public void setDirection(Vector3d v) {
		Vector3d x = new Vector3d(v);
		x.normalize();
		this.direction = x;
	}
	/**
	 * Set this chemical field as the goal field.
	 */
	public void setGoal(BSimChemicalField goal) {
		this.goal = goal;
		setMemoryDuration(shortTermMemoryDuration, longTermMemoryDuration);
		memory = new double[sim.timesteps(getMemoryDuration())];
		for(int i=0;i<memory.length;i++) memory[i] = goal.getConc(position);
	}
	public void setMemoryDuration(double shortTermMemoryDuration, double longTermMemoryDuration) {
		this.shortTermMemoryDuration = shortTermMemoryDuration;
		this.shortTermMemoryLength = sim.timesteps(shortTermMemoryDuration);
		this.longTermMemoryDuration = longTermMemoryDuration;
		this.longTermMemoryLength = sim.timesteps(longTermMemoryDuration);
	}

	public Vector3d getDirection() { return direction; }
	public MotionState getMotionState() { return motionState; }
	public double getMemoryDuration() { return shortTermMemoryDuration + longTermMemoryDuration; }

	/**
	 * Applies the flagellar force.
	 */
	public void flagellarForce() {
		Vector3d f = new Vector3d();
		f.scale(forceMagnitude, direction);
		addForce(f);
	}

	/**
	 * Causes the cell to rotate such that Var(theta(dt)) = 4*D*dt.
	 */
	public void rotationalDiffusion() {
		double dTheta = rng.nextGaussian()*Math.sqrt(4*BSim.BOLTZMANN*sim.getTemperature()*sim.getDt()/rotationalStokesCoefficient())*Math.pow(10,9);
		BSimUtils.rotatePerp(direction, dTheta);
	}

	public double rotationalStokesCoefficient() {
		return 8.0*Math.PI*sim.getVisc()*Math.pow(radius,3); // Pa sec microns^3
	}

	/**
	 * Return a tumble angle in radians distributed according to Fig. 3, 'Chemotaxis
	 * in Escherichia Coli', Berg et al. (claim from 'AgentCell: a digital single-cell
	 * assay for bacterial chemotaxis', Emonet et al.).
	 */
	public double tumbleAngle() {
		double tumbleShape = 4;
		double tumbleScale = 18.32;
		double tumbleLocation = -4.60;

		double tumbleAngle;
		do {
			tumbleAngle = BSimUtils.sampleGamma(tumbleShape, tumbleScale) + tumbleLocation;
		} while (tumbleAngle > 180);

		return Math.toRadians(tumbleAngle);
	}

	/**
	 * p23, 'Strategies for Chemotaxis', Schnizter, Berg et al.
	 * Compare two sequential averages.
	 */
	public boolean movingUpGradient() {
		double shortTermCounter = 0, shortTermMean = 0;
		double longTermCounter = 0, longTermMean = 0;

		System.arraycopy(memory, 0, memory, 1, memory.length - 1);
		memory[0] = goal.getConc(position);

		for(int i=0; i<memory.length; i++) {
			if(i < shortTermMemoryLength) {
				shortTermCounter = shortTermCounter + memory[i];
			}
			else {
				assert (i < shortTermMemoryLength + longTermMemoryLength);
				longTermCounter = longTermCounter + memory[i];
			}
		}
		shortTermMean = shortTermCounter/shortTermMemoryLength;
		longTermMean = longTermCounter/longTermMemoryLength;

		return shortTermMean - longTermMean > sensitivity;
	}

	/*
	 * GROWTH and REPLICATION
	 *
	 * See 'Cell Shape Dynamics in Escherichia Coli', Reshes et al
	 *
	 */

    /******************************************************************************************************************/

    protected LambdaOdeSystem lamSys;
    protected double[] y, yNew;
    public double lambda;

    class LambdaOdeSystem implements BSimOdeSystem {

        @Override
        public double[] derivativeSystem(double x, double[] y) {
            double[] dy = new double[22]; //lambda could be in the last entry (23) (we calculate lambda while solving the ODEs anyway)

            //PARAMETERS
            // 0: thetar= 426.8693338968694;  % ribosome transcription threshold [molecs/cell]
            // 1: k_cm= 0.005990373118888;    % chloramphenicol-binding rate [1/(min microM)]
            // 2: s0= 1.0e4;                  % external nutrient [molecs]
            // 3: gmax= 1260.0;               % max. transl. elongation rate [aa/(min molecs)]
            // 4: cl= 0;                      % ribosome-bound mRNA??????????????
            // 5: thetax= 4.379733394834643;  % non-ribosomal transcription threshold [molecs/cell]
            // 6: Kt= 1.0e3;                  % nutrient import threshold [molecs]
            // 7: M= 1.0e8;                   % total cell mass [aa]
            // 8: we= 4.139172187824451;      % max. enzyme transcription rate [molecs/(min cell)]
            // 9: Km= 1.0e3;                  % enzymatic threshold [molecs/cell]
            //10: vm= 5800.0;                 % max. enzymatic rate [1/min]
            //11: nx= 300.0;                  % length of non-ribosomal protein [aa/molecs]
            //12: Kq= 1.522190403737490e+05;  % q-autoinhibition threshold [molecs/cell]
            //13: Kp= 180.1378030928276;      % K_p = k_1*k_2/(k_(-1) + k_2) (translation)
            //14: vt= 726.0;                  % max. nutrient import rate [1/min]
            //15: wr= 929.9678874564831;      % max. ribosome transcription rate [molecs/(min cell)]
            //16: wq= 948.9349882947897;      % max. q-transcription rate [molecs/(min cell)]
            //17: wp= 0.0;                    % ?????? for another gene
            //18: nq= 4;                      % q-autoinhibition Hill coeff. [none]
            //19: nr= 7549.0;                 % ribosome length [aa/molecs]
            //20: ns= 0.5;                    % nutrient efficiency [none]
            double[] para = new double[] {426.8693338968694, 0.005990373118888, 10000, 1260.0, 0.0, 4.379733394834643, 1000, 100000000, 4.139172187824451, 1000, 5800.0, 300.0, 152219.0403737490, 180.1378030928276, 726.0, 929.9678874564831, 948.9349882947897, 0.0, 4, 7549.0, 0.5};

            //RATES
            //define rate constants
            //0: b= 0;       % ?????
            //1: dm= 0.1;    % mRNA-degradation rate
            //2: kb= 1;      % mRNA-ribosome binding rate
            //3: ku= 1.0;    % mRNA-ribosome unbinding rate
            //4: f= cl*k_cm; % ????? for chloramphenicol
            double[] rates = new double[] {0.0, 0.1, 1.0, 1.0, para[4]*para[1]};

            //ODEs
            // translating ODE system from matlab-code
            double Kgamma, gamma, R_t, ttrate, nucat, lam;

            Kgamma= para[3]/para[13];               // transl. elongation treshold (half-max. elongation)
            gamma= para[3]*y[21]/(Kgamma + y[21]);  // transl. elongation rate
            R_t = y[3] + y[0] + y[2] + y[4] + y[6]; // sum of translating ribosomes (sum c_x)
            ttrate= R_t*gamma;                      // translation-transcription rate?
            lam = ttrate/para[7];                   // growth-rate (M is total protein mass)
            //setLambda(lam);
            setSurfaceAreaGrowthRate(lam);
            //System.out.printf("Lambda in 'derivative system' %.4f ", lam);
            //fr= nr*(r + rmr + rmp + rmt + rmm + rmq + zmr + zmp + zmt + zmm + zmq) / ...    % r+rmr+rmp+rmt+rmm+rmq is total amount of ribosomes
            //        ( nr*(r + rmr + rmp + rmt + rmm + rmq + zmr + zmp + zmt + zmm + zmq) ...
            //+ nx * (p + q + et + em));
            nucat= y[1]*para[10]*y[16]/(para[9] + y[16]);

            dy[0]= +rates[2]*y[20]*y[19]+rates[0]*y[8]-rates[3]*y[0]-gamma/para[19]*y[0]-rates[4]*y[0]-lam*y[0];
            dy[1]= +gamma/para[11]*y[6]-lam*y[1];
            dy[2]= +rates[2]*y[20]*y[18]+rates[0]*y[9]-rates[3]*y[2]-gamma/para[11]*y[2]-rates[4]*y[2]-lam*y[2];
            dy[3]= +rates[2]*y[20]*y[17]+rates[0]*y[10]-rates[3]*y[3]-gamma/para[11]*y[3]-rates[4]*y[3]-lam*y[3];
            dy[4]= +rates[2]*y[20]*y[12]+rates[0]*y[11]-rates[3]*y[4]-gamma/para[11]*y[4]-rates[4]*y[4]-lam*y[4];
            dy[5]= +gamma/para[11]*y[4]-lam*y[16];
            dy[6]= +rates[2]*y[20]*y[13]+rates[0]*y[7]-rates[3]*y[6]-gamma/para[11]*y[6]-rates[4]*y[6]-lam*y[6];
            dy[7]= +rates[4]*y[6]-rates[0]*y[7]-lam*y[7];
            dy[8]= +rates[4]*y[0]-rates[0]*y[8]-lam*y[8];
            dy[9]= +rates[4]*y[2]-rates[0]*y[9]-lam*y[9];
            dy[10]= +rates[4]*y[3]-rates[0]*y[10]-lam*y[10];
            dy[11]= +rates[4]*y[4]-rates[0]*y[11]-lam*y[11];
            dy[12]= +(para[8]*y[21]/(para[5] + y[21]))+rates[3]*y[4]+gamma/para[11]*y[4]-rates[2]*y[20]*y[12]-rates[1]*y[12]-lam*y[12];
            dy[13]= +(para[8]*y[21]/(para[5] + y[21]))+rates[3]*y[6]+gamma/para[11]*y[6]-rates[2]*y[20]*y[13]-rates[1]*y[13]-lam*y[13];
            dy[14]= +gamma/para[11]*y[3]-lam*y[14];
            dy[15]= +gamma/para[11]*y[2]-lam*y[15];
            dy[16]= +(y[5]*para[14]*para[2]/(para[6] + para[2]))-nucat-lam*y[16];
            dy[17]= +(para[16]*y[21]/(para[5] + y[21])/(1 + Math.pow(y[14]/para[12], para[18])))+rates[3]*y[3]+gamma/para[11]*y[3]-rates[2]*y[20]*y[17]-rates[1]*y[17]-lam*y[17];
            dy[18]= +(para[17]*y[21]/(para[5] + y[21]))+rates[3]*y[2]+gamma/para[11]*y[2]-rates[2]*y[20]*y[18]-rates[1]*y[18]-lam*y[18];
            dy[19]= +(para[15]*y[21]/(para[0] + y[21]))+rates[3]*y[0]+gamma/para[19]*y[0]-rates[2]*y[20]*y[19]-rates[1]*y[19]-lam*y[19];
            dy[20]= +rates[3]*y[0]+rates[3]*y[4]+rates[3]*y[6]+rates[3]*y[2]+rates[3]*y[3]+gamma/para[19]*y[0]+gamma/para[19]*y[0]+gamma/para[11]*y[4]+gamma/para[11]*y[6]+gamma/para[11]*y[2]+gamma/para[11]*y[3]-rates[2]*y[20]*y[19]-rates[2]*y[20]*y[12]-rates[2]*y[20]*y[13]-rates[2]*y[20]*y[18]-rates[2]*y[20]*y[17]-lam*y[20];
            dy[21]= +para[20]*nucat-ttrate-lam*y[21];
            //dy[22] = lambda;
            return dy;
        }

        public double[] getICs() { // overloading for children!
            double[] ics = new double[22];
            // 0: rmr
            // 1: em
            // 2: rmp
            // 3: rmq
            // 4: rmt
            // 5: et
            // 6: rmm
            // 7: zmm
            // 8: zmr
            // 9: zmp
            //10: zmq
            //11: zmt
            //12: mt
            //13: mm
            //14: q
            //15: p
            //16: si
            //17: mq
            //18: mp
            //19: mr
            //20: r
            //21: a
            //22: lambda


            ics[0] = 355.9405;
            ics[1] = 10369;   // metabolic enzyme
            ics[2] = 0.0;   // ribosome binding r + m_p -> c_p
            ics[3] = 1673.3;   // ribosome binding r + m_q -> c_q
            ics[4] = 67.3151;   // ribosome binding r + m_t -> c_t
            ics[5] = 10369;    // transporter enzyme
            ics[6] = 67.3151;   // ribosome binding r + m_m
            ics[7] = 0.0;
            ics[8] = 0.0;
            ics[9] = 0.0;
            ics[10] = 0.0;
            ics[11] = 0.0;
            ics[12] = 9.1942;    // free mRNA of transporter enzymes
            ics[13] = 9.1942;    // free mRNA of metabolic enzymes
            ics[14] = 257760;     // house-keeping proteins (growth-independent)
            ics[15] = 0.0;     //??????????
            ics[16] = 128.4047;    // internal nutrient
            ics[17] = 228.5501;    // free mRNA of house-keeping proteins q
            ics[18] = 0.0;    // free mRNA of ???????
            ics[19] = 24.7372;    // free mRNA of ribosomes r
            ics[20] = 15.0902;  // ribosomes
            ics[21] = 2.3441; // energy

			/*
			ics[0] = 0.0;
			ics[1] = 0.0;   // metabolic enzyme
			ics[2] = 0.0;   // ribosome binding r + m_p -> c_p
			ics[3] = 0.0;   // ribosome binding r + m_q -> c_q
			ics[4] = 0.0;   // ribosome binding r + m_t -> c_t
			ics[5] = 0.0;    // transporter enzyme
			ics[6] = 0.0;   // ribosome binding r + m_m
			ics[7] = 0.0;
			ics[8] = 0.0;
			ics[9] = 0.0;
			ics[10] = 0.0;
			ics[11] = 0.0;
			ics[12] = 0.0;    // free mRNA of transporter enzymes
			ics[13] = 0.0;    // free mRNA of metabolic enzymes
			ics[14] = 0.0;     // house-keeping proteins (growth-independent)
			ics[15] = 0.0;     //??????????
			ics[16] = 0.0;    // internal nutrient
			ics[17] = 0.0;    // free mRNA of house-keeping proteins q
			ics[18] = 0.0;    // free mRNA of ???????
			ics[19] = 0.0;    // free mRNA of ribosomes r
			ics[20] = 10.0;  // ribosomes
			ics[21] = 1000.0; // energy
			*/
            //ics[22] = 0.5;   // LAMBDA
            return ics;
        }

        public int getNumEq(){
            return 22;
        }
    }

	protected double surfaceAreaGrowthRate = 0; // microns^2/s
	public void setSurfaceAreaGrowthRate() { surfaceAreaGrowthRate = 0.025; } // 5 typical vesicle surface areas/second
	public void setSurfaceAreaGrowthRate(double d) { surfaceAreaGrowthRate = d; }

	/*
	 * Generation time T
	 * 	= S(T) - S(0) / surfaceAreaGrowthRate
	 * 	= S(T)/(2*surfaceAreaGrowthRate)
	 * 	= (2*pi*replicationRadius^2)/surfaceAreaGrowthRate
	 */
	protected double replicationRadius = Math.sqrt(2); // microns, so birth radius = 1 micron
	protected void setReplicationRadius(double r) { replicationRadius = r; }
	/** The external list of children. Required when bacteria reach the replicationRadius */
	@SuppressWarnings("rawtypes")
	protected Vector childList;
	public void setChildList(@SuppressWarnings("rawtypes") Vector v) { childList = v; }

	/** Sets the radius so that the surface area of the bacterium is randomly distributed between surfaceArea(replicationRadius)/2 and surfaceArea(replicationRadius) */
	public void setRadius() {
		setRadiusFromSurfaceArea(surfaceArea(replicationRadius)/2 + Math.random()*surfaceArea(replicationRadius)/2);
	}

	public void grow() {
        System.out.printf(" Time %.2f ", sim.getTime());
        System.out.printf(" Growth rate in grow %.4f ", surfaceAreaGrowthRate);

		double dS = surfaceAreaGrowthRate*sim.getDt();
		setRadiusFromSurfaceArea(getSurfaceArea() + dS);

		if(pVesicle > 0 && Math.random() < pVesicle*(dS/typicalVesicleSurfaceArea))
			vesiculate();

		if (radius > replicationRadius)
			replicate();
	}

	@SuppressWarnings("unchecked")
	public void replicate() {
		setRadiusFromSurfaceArea(surfaceArea(replicationRadius)/2);
		BSimBacterium child = new BSimBacterium(sim, new Vector3d(position));
		child.setRadius(radius);
		child.setSurfaceAreaGrowthRate(surfaceAreaGrowthRate);
		child.setChildList(childList);
		/* Overwrite to allow inheritance of other properties */
		childList.add(child);
	}


	/*
	 * VESICULATION
	 */
	protected double vesicleRadius = 0.02; // microns
	public double vesicleRadius() { return vesicleRadius; }
	public void vesicleRadius(double d) { vesicleRadius = d; }
	protected double typicalVesicleSurfaceArea = 0.005; // microns^2
	/*
	* 'Some Characteristics of the Outer Membrane Material Released by Growing
	* Enterotoxigenic Escherichia Coli', Gankema et al.:
	* 'The medium vesicles.. accounted for 3 to 5% of the total cellular outer membrane'
	* Mean growth before producing a vesicle = 1/pVesicle
	* For pVesicle = 5%, Mean growth before production = 20 vesicle surface areas ~ 4 seconds
	*/
	/** Probability per typical vesicle surface area growth of producing a vesicle */
	protected double pVesicle = 0; // 1/(typical vesicle surface areas)
	public void pVesicle(double d) { pVesicle = d; }
	/** The external list of vesicles. Required when bacteria vesiculate */
	@SuppressWarnings("rawtypes")
	protected Vector vesicleList;
	public void setVesicleList(@SuppressWarnings("rawtypes") Vector v) { vesicleList = v; }

	@SuppressWarnings("unchecked")
	public void vesiculate() {
		double r = vesicleRadius();
		vesicleList.add(new BSimVesicle(sim, new Vector3d(position), r));
		setRadiusFromSurfaceArea(getSurfaceArea()-surfaceArea(r));
	}



	/**
	 * Creates a RUNNING bacterium at the specified position, facing in a
	 * random direction
	 */
	public BSimBacterium(BSim sim, Vector3d position) {
		super(sim, position, 1); // default radius 1 micron
		setMotionState(MotionState.RUNNING);
		setDirection(new Vector3d(0.5-Math.random(),0.5-Math.random(),0.5-Math.random()));
        lamSys = new LambdaOdeSystem();
        y = lamSys.getICs();
		/*
        MWNumericArray yy = null;
		MWNumericArray t0 = null;
		MWNumericArray tf = null;
		Object[] ynew = null;
		Class1 theMagic = null;*/
		//theMagic = new Class1();

        //lambda = getLambda();

        //radau5 solver
		//bsim.object.BSimBacterium b = new bsim.object.BSimBacterium();
        // Construct ODEs for solving all contact constraints
        //ODE odes = new RepressilatorODESystem();
		/*
        // Solve the contact constraint system
        final double initTime = sim.getTime(), maxTime = sim.getTime()+sim.getDt();

        // Solver and its parameters
        InterpolatorEventSolver solver = new InterpolatorEventSolver(new Radau5(), odes);

        double stepSize = 0.01; // The initial step size (used by fixed step methods)
        double plotStepSize = 0.1; // The step size for plotting the solution
        final double absTol = 1.0e-6, relTol = 1.0e-3; // The tolerance for adaptive methods

        // Initialize and customize the solver
		solver.initialize(stepSize);       // This step size affects the solver internal step size
		solver.setStepSize(plotStepSize);  // This step size is the reading step size
		solver.setTolerances(absTol, relTol);
		solver.setHistoryLength(Double.POSITIVE_INFINITY); // Don't recall all past values
		solver.setEnableExceptions(false); // Do not throw exceptions when an error occurs

		double[] state = odes.getState();

        // Solve for the whole [initTime,maxTime] interval at once
		while (solver.getCurrentTime() < maxTime) {
                //            System.out.println("Hello - stepping solver");
                //System.out.println("Advancing the solution from " + solver.getCurrentTime());

			solver.step();
			if (solver.getErrorCode() != InterpolatorEventSolver.ERROR.NO_ERROR) {
				System.err.println("Error when advancing the solution from " + solver.getCurrentTime());
				return;
            }
        }

        // Compute max error at each plot point
        StateHistory history = solver.getStateHistory();
		*/
        /*{  // plot the graphs
            // Prepare the graphics
            Dataset stripChart = new Dataset(Color.BLUE, Color.BLUE, true);
            stripChart.setMarkerShape(Dataset.NO_MARKER);

            double time = initTime;
            double[] interpolated = new double[odes.getState().length];
            while (time <= maxTime) {
                history.interpolate(time, interpolated);
                stripChart.append(time, interpolated[0]);
                time += plotStepSize;
            }

            PlottingPanel plottingPanel = new PlottingPanel("time", "state", "ODE Test"); //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
            DrawingFrame plottingFrame = new DrawingFrame("ODE Test", plottingPanel); //$NON-NLS-1$
            plottingFrame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
            plottingPanel.addDrawable(stripChart);

            plottingPanel.render();
            plottingFrame.setLocation(0, 0);
            plottingFrame.setSize(700, 700);
            plottingFrame.setVisible(true);
        }*/
    }
	/*
    class RepressilatorODESystem implements ODE {

        // Everything gets upset without mState
        private double[] mState;

        // Let's try to redefine mState in the constructor
        // (useful in case we want some complicated ICs...)
        public RepressilatorODESystem(){
            this.mState = lamSys.getICs();
        }

        /**
         * Initial conditions
         */
        /*
        @Override
        public double[] getState() {
            return mState;
        }
		*/
        /**
        * Classic Garcia-Ojalvo repressilator ODEs
        */
        /*
        @Override
        public void getRate(double[] y, double[] dy) {

            double[] para = new double[] {426.8693338968694, 0.005990373118888, 10000, 1260.0, 0.0, 4.379733394834643, 1000, 100000000, 4.139172187824451, 1000, 5800.0, 300.0, 152219.0403737490, 180.1378030928276, 726.0, 929.9678874564831, 948.9349882947897, 0.0, 4, 7549.0, 0.1};

            //RATES
            //define rate constants
            //0: b= 0;       % ?????
            //1: dm= 0.1;    % mRNA-degradation rate
            //2: kb= 1;      % mRNA-ribosome binding rate
            //3: ku= 1.0;    % mRNA-ribosome unbinding rate
            //4: f= cl*k_cm; % ????? for chloramphenicol
            double[] rates = new double[] {0.0, 0.1, 1.0, 1.0, para[4]*para[1]};

            //ODEs
            // translating ODE system from matlab-code
            double Kgamma, gamma, R_t, ttrate, nucat, lam;

            Kgamma= para[3]/para[13];               // transl. elongation treshold (half-max. elongation)
            gamma= para[3]*y[21]/(Kgamma + y[21]);  // transl. elongation rate
            R_t = y[3] + y[0] + y[2] + y[4] + y[6]; // sum of translating ribosomes (sum c_x)
            ttrate= R_t*gamma;                      // translation-transcription rate?
            lam = ttrate/para[7];                   // growth-rate (M is total protein mass)
            setSurfaceAreaGrowthRate(lam);
            System.out.printf("Lambda in 'derivative system' %.4f at t %.2f ", lam, sim.getTime());
            nucat= y[1]*para[10]*y[16]/(para[9] + y[16]);

            dy[0]= rates[2]*y[20]*y[19]+rates[0]*y[8]-rates[3]*y[0]-gamma/para[19]*y[0]-rates[4]*y[0]-lam*y[0];
            dy[1]= gamma/para[11]*y[6]-lam*y[1];
            dy[2]= rates[2]*y[20]*y[18]+rates[0]*y[9]-rates[3]*y[2]-gamma/para[11]*y[2]-rates[4]*y[2]-lam*y[2];
            dy[3]= rates[2]*y[20]*y[17]+rates[0]*y[10]-rates[3]*y[3]-gamma/para[11]*y[3]-rates[4]*y[3]-lam*y[3];
            dy[4]= rates[2]*y[20]*y[12]+rates[0]*y[11]-rates[3]*y[4]-gamma/para[11]*y[4]-rates[4]*y[4]-lam*y[4];
            dy[5]= gamma/para[11]*y[4]-lam*y[16];
            dy[6]= rates[2]*y[20]*y[13]+rates[0]*y[7]-rates[3]*y[6]-gamma/para[11]*y[6]-rates[4]*y[6]-lam*y[6];
            dy[7]= rates[4]*y[6]-rates[0]*y[7]-lam*y[7];
            dy[8]= rates[4]*y[0]-rates[0]*y[8]-lam*y[8];
            dy[9]= rates[4]*y[2]-rates[0]*y[9]-lam*y[9];
            dy[10]= rates[4]*y[3]-rates[0]*y[10]-lam*y[10];
            dy[11]= rates[4]*y[4]-rates[0]*y[11]-lam*y[11];
            dy[12]= (para[8]*y[21]/(para[5] + y[21]))+rates[3]*y[4]+gamma/para[11]*y[4]-rates[2]*y[20]*y[12]-rates[1]*y[12]-lam*y[12];
            dy[13]= (para[8]*y[21]/(para[5] + y[21]))+rates[3]*y[6]+gamma/para[11]*y[6]-rates[2]*y[20]*y[13]-rates[1]*y[13]-lam*y[13];
            dy[14]= gamma/para[11]*y[3]-lam*y[14];
            dy[15]= gamma/para[11]*y[2]-lam*y[15];
            dy[16]= (y[5]*para[14]*para[2]/(para[6] + para[2]))-nucat-lam*y[16];
            dy[17]= (para[16]*y[21]/(para[5] + y[21])/(1 + Math.pow((y[14]/para[12]), para[18])))+rates[3]*y[3]+gamma/para[11]*y[3]-rates[2]*y[20]*y[17]-rates[1]*y[17]-lam*y[17];
            dy[18]= (para[17]*y[21]/(para[5] + y[21]))+rates[3]*y[2]+gamma/para[11]*y[2]-rates[2]*y[20]*y[18]-rates[1]*y[18]-lam*y[18];
            dy[19]= (para[15]*y[21]/(para[0] + y[21]))+rates[3]*y[0]+gamma/para[19]*y[0]-rates[2]*y[20]*y[19]-rates[1]*y[19]-lam*y[19];
            dy[20]= rates[3]*y[0]+rates[3]*y[4]+rates[3]*y[6]+rates[3]*y[2]+rates[3]*y[3]+gamma/para[19]*y[0]+gamma/para[19]*y[0]+gamma/para[11]*y[4]+gamma/para[11]*y[6]+gamma/para[11]*y[2]+gamma/para[11]*y[3]-rates[2]*y[20]*y[19]-rates[2]*y[20]*y[12]-rates[2]*y[20]*y[13]-rates[2]*y[20]*y[18]-rates[2]*y[20]*y[17]-lam*y[20];
            dy[21]= para[20]*nucat-ttrate-lam*y[21];
            //this.mState = dy;
        }
    }*/

    public double getLambda(){
        return surfaceAreaGrowthRate;
    }


	@Override
	public void action() {
		super.action();

		switch(motionState) {
		case RUNNING:
			if(Math.random() < pEndRun()*sim.getDt())
				motionState = MotionState.TUMBLING;
			break;
		case TUMBLING:
			if(Math.random() < pEndTumble()*sim.getDt()) {
				/* Change the direction at the end of a tumble phase */
				BSimUtils.rotatePerp(direction, tumbleAngle());
				motionState = MotionState.RUNNING;
			}
			break;
		default:
			assert false : motionState;
		}

		if(motionState == MotionState.RUNNING) {
			rotationalDiffusion();
			flagellarForce();
		}
		MWNumericArray yy = null;
		MWNumericArray t0 = null;
		MWNumericArray tf = null;
		Object[] ynew = null;
		Class1 theMagic = null;

		yy = new MWNumericArray(y, MWClassID.DOUBLE);
		t0 = new MWNumericArray(Double.valueOf(sim.getTime()), MWClassID.DOUBLE);
		tf = new MWNumericArray(Double.valueOf(sim.getDt()+sim.getTime()), MWClassID.DOUBLE);

        //yNew = theMagic.stiffODE(22, yy, t0, tf);
        y = yNew;

		/*
		// Construct ODEs for solving all contact constraints
		ODE odes = new RepressilatorODESystem();

		// Solve the contact constraint system
		final double initTime = sim.getTime(), maxTime = 10000;

		// Solver and its parameters
		InterpolatorEventSolver solver = new InterpolatorEventSolver(new Radau5(), odes);

		double stepSize = 0.0001; // The initial step size (used by fixed step methods)
		double plotStepSize = 0.1; // The step size for plotting the solution
		final double absTol = 1.0e-6, relTol = 1.0e-3; // The tolerance for adaptive methods

		// Initialize and customize the solver
		solver.initialize(stepSize);       // This step size affects the solver internal step size
		solver.setStepSize(plotStepSize);  // This step size is the reading step size
		solver.setTolerances(absTol, relTol);
		solver.setHistoryLength(Double.POSITIVE_INFINITY); // Don't recall all past values
		solver.setEnableExceptions(true); // Do not throw exceptions when an error occurs

		//double[] state = odes.getState();

		// Solve for the whole [initTime,maxTime] interval at once
		while (solver.getCurrentTime() < maxTime) {
			//System.out.println("Hello - stepping solver");
			//System.out.println("Advancing the solution from " + solver.getCurrentTime());

			solver.step();
			if (solver.getErrorCode() != InterpolatorEventSolver.ERROR.NO_ERROR) {
				System.err.println("Error when advancing the solution from " + solver.getCurrentTime());
				return;
			}
		}*/

		grow();
	}
	

}