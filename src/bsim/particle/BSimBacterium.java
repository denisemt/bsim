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

//import com.kabouterlabs.jodeint.codepack.CodepackLibrary;
//import org.bridj.Pointer;

import com.mathworks.toolbox.javabuilder.*;
import IMPLsolver.*;

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
	/******************************************************/
	public double[] y;
	public double getY(int i){return this.y[i];}

	public double[] getICs(){
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

		//half of steady state

		ics[0] = 355.9405 / 2.0;
		ics[1] = 10369 / 2.0;   // metabolic enzyme
		ics[2] = 0.0;   // ribosome binding r + m_p -> c_p
		ics[3] = 1673.3 / 2.0;   // ribosome binding r + m_q -> c_q
		ics[4] = 67.3151 / 2.0;   // ribosome binding r + m_t -> c_t
		ics[5] = 10369 / 2.0;    // transporter enzyme
		ics[6] = 67.3151 / 2.0;   // ribosome binding r + m_m
		ics[7] = 0.0;
		ics[8] = 0.0;
		ics[9] = 0.0;
		ics[10] = 0.0;
		ics[11] = 0.0;
		ics[12] = 9.1942 / 2.0;    // free mRNA of transporter enzymes
		ics[13] = 9.1942 / 2.0;    // free mRNA of metabolic enzymes
		ics[14] = 257760 / 2.0;     // house-keeping proteins (growth-independent)
		ics[15] = 0.0;     //??????????
		ics[16] = 128.4047 / 2.0;    // internal nutrient
		ics[17] = 228.5501 / 2.0;    // free mRNA of house-keeping proteins q
		ics[18] = 0.0;    // free mRNA of ???????
		ics[19] = 24.7372 / 2.0;    // free mRNA of ribosomes r
		ics[20] = 15.0902 / 2.0;  // ribosomes
		ics[21] = 2.3441 / 2.0; // energy

			/*
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
			ics[1] = 0.0;    // metabolic enzyme
			ics[2] = 0.0;    // ribosome binding r + m_p -> c_p
			ics[3] = 0.0;    // ribosome binding r + m_q -> c_q
			ics[4] = 0.0;    // ribosome binding r + m_t -> c_t
			ics[5] = 0.0;    // transporter enzyme
			ics[6] = 0.0;    // ribosome binding r + m_m
			ics[7] = 0.0;
			ics[8] = 0.0;
			ics[9] = 0.0;
			ics[10] = 0.0;
			ics[11] = 0.0;
			ics[12] = 0.0;    // free mRNA of transporter enzymes
			ics[13] = 0.0;    // free mRNA of metabolic enzymes
			ics[14] = 0.0;    // house-keeping proteins (growth-independent)
			ics[15] = 0.0;    //??????????
			ics[16] = 0.0;    // internal nutrient
			ics[17] = 0.0;    // free mRNA of house-keeping proteins q
			ics[18] = 0.0;    // free mRNA of ???????
			ics[19] = 0.0;    // free mRNA of ribosomes r
			ics[20] = 10.0;   // ribosomes
			ics[21] = 1000.0; // energy
			*/
		return ics;
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
        //System.out.printf(" Time %.2f ", sim.getTime());
        //System.out.printf(" Growth rate in grow %.10f ", surfaceAreaGrowthRate);

		double dS = surfaceAreaGrowthRate*sim.getDt();
		setRadiusFromSurfaceArea(getSurfaceArea() + dS);
		System.out.print("current radius (from surfaceArea): ");
		System.out.print(radius);

		if(pVesicle > 0 && Math.random() < pVesicle*(dS/typicalVesicleSurfaceArea))
			vesiculate();

		// if division: take half of the nutrients
		if (radius > replicationRadius) {
			for (int k = 0; k < 22; k++) {
				this.y[k] = this.y[k]/2.0;
			}
            replicate();
        }
	}

	@SuppressWarnings("unchecked")
	public void replicate() {
		setRadiusFromSurfaceArea(surfaceArea(replicationRadius)/2);
		// if I change BSimBacterium constructor to decide if its the very first mothercell (ICs half of steady state)
        // or every other cell afterwards (ICs half of y of mother), distribution of nutrients between mother daughter
        // should be possible
		BSimBacterium child = new BSimBacterium(sim, new Vector3d(position));
		child.setRadius(radius);
		// maybe set here half of vector y of mother for child (as ICS)
		// and set surfaceAreaGrowthRate calculated from half of y as growthrate for children
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
        y = getICs();

    }

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

		//**************************************************************************+
		// use matlab control to solve the ODE system inside each bacterium
		double t0 = sim.getTime()/60.0;
		double tf = (sim.getDt() + sim.getTime())/60.0;

		double[] YNEW = new double[22];

		try{
			Class1 theMagic = new IMPLsolver.Class1();
			Object[] yNew = theMagic.IMPLsolver(22, y, t0, tf);
			for (int k = 0; k < 22; k++) {
				YNEW[k] = Double.valueOf(yNew[k].toString());
			}
			MWArray.disposeArray(yNew);
			theMagic.dispose();
		}catch (Exception e){
			System.out.println("Exception: " + e.toString());
		}

		y = YNEW;

		// 0: rmr, 1: em, 2: rmp, 3: rmq, 4: rmt, 5: et, 6: rmm, 7: zmm, 8: zmr, 9: zmp, 10: zmq, 11: zmt, 12: mt,
		// 13: mm, 14: q, 15: p, 16: si, 17: mq, 18: mp, 19: mr, 20: r, 21: a

		double Mref = 1e8;
		double R_t = y[3] + y[0] + y[2] + y[4] + y[6];           // sum of translating ribosomes (sum c_x)
		//aatot= nx*(q+et+em)+nr*(rmq+rmr+rmp+rmt+rmm);
		double M = 300.0 * (y[14] + y[5] + y[1]) + 7459.0 * R_t; // total cell mass
		double Kgamma= 1260.0 / 180.1378030928276*M/Mref;        // transl. elongation treshold (half-max. elongation)
		double gamma = 1260.0 * y[21] / (Kgamma + y[21]);        // transl. elongation rate

		double ttrate = R_t * gamma;                             // translation-transcription rate?
		double lam = ttrate / M;                                 // growth-rate (M is total protein mass)
		setSurfaceAreaGrowthRate(lam);
		System.out.print(" Lambda: ");
		System.out.print(lam);

		System.out.print(" total mass: ");
		System.out.print(M);

		double Volume = 1e-8 * M;
		double r = Math.pow(3.0 / 4.0 * Volume / Math.PI, (1.0 / 3));
		System.out.print(" radius (from M): ");
		System.out.print(r);
		surfaceArea(r);

		grow();
	}
	

}