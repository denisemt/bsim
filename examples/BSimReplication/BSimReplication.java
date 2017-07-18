package BSimReplication;

import java.awt.Color;
import java.util.Vector;

import javax.vecmath.Vector3d;

import org.opensourcephysics.numerics.ODE;
import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimTicker;
import bsim.draw.BSimP3DDrawer;
import bsim.particle.BSimBacterium;

import com.mathworks.toolbox.javabuilder.*;
import makesqr.*;
import stiffODE.*;
import stiffODE.Class1;

/**
 * Growth and replication test.</br>
 * 
 * Simple example of growing and replicating population of bacteria swimming in an environment.
 * Visualisation and testing of built in {@link BSimBacterium} growth/replication methods.
 */
public class BSimReplication {

	public static void main(String[] args) {

		/*********************************************************
		 * Set simulation properties
		 */
		/*
		BSim sim = new BSim();
		sim.setDt(10000);
		sim.setSimulationTime(10000);
		//BSimReplication b = new BSimReplication();
		//ODE odes = b.new RepressilatorODESystem();
						
		/*********************************************************
		 * Set up the bacteria
		 */
		/*
		final Vector<BSimBacterium> bacteria = new Vector<BSimBacterium>();
		final Vector<BSimBacterium> children = new Vector<BSimBacterium>();
		while(bacteria.size() < 1) {		
			BSimBacterium b = new BSimBacterium(sim, new Vector3d(Math.random()*sim.getBound().x, Math.random()*sim.getBound().y, Math.random()*sim.getBound().z));
			b.setRadius();
			b.setSurfaceAreaGrowthRate(2.0);
			b.setChildList(children);
			bacteria.add(b);		
		}

		/*********************************************************
		 * Set up the ticker
		 */
		/*
		sim.setTicker(new BSimTicker() {
			@Override
			public void tick() {
				for(BSimBacterium b : bacteria) {
					b.action();		
					b.updatePosition();					
				}
				bacteria.addAll(children);
				children.clear();
				System.out.printf("number %d ", bacteria.size());
				System.out.printf("\n");
			}		
		});

		/*********************************************************
		 * Set up the drawer
		 */
		/*
		BSimP3DDrawer drawer = new BSimP3DDrawer(sim, 800,600) {
			@Override
			public void scene(PGraphics3D p3d) {						
				for(BSimBacterium b : bacteria) {
					draw(b, Color.GREEN);
				}			
			}
		};
		sim.setDrawer(drawer);				

		// Run the simulation
		sim.preview();
		*/
		//TEST MATLAB MAGIC
		/*
		String[] arg = {String.valueOf(5)};

		MWNumericArray n = null;
		Object[] result = null;
		Class1 theMagic = null;

		if (args.length == 0)
		{
			System.out.println("Error: must input a positive integer");
			//return;
		}

		try
		{
			n = new MWNumericArray(Double.valueOf(arg[0]), MWClassID.DOUBLE);
			//int m = 1;
			theMagic = new Class1();

			result = theMagic.makesqr(1, n);
			System.out.println(result[0]);
		}
		catch (Exception e)
		{
			System.out.println("Exception: " + e.toString());
		}
		finally
		{
			MWArray.disposeArray(n);
			MWArray.disposeArray(result);
			//System.out.printf("Result %.4f", result[0]);
			theMagic.dispose();
		}
		*/

		MWNumericArray yy = null;
		MWNumericArray t0 = null;
		MWNumericArray tf = null;
		Object[] yNew = null;
		Class1 theMagic = null;

		System.out.printf("Marke 1 ");

		try {
			double[] y = new double[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10.0, 1000.0};
			System.out.print("y: ");
			System.out.println(y[3]);
			System.out.print("\n");

			double time = 0.0;
			double endTime = 1e7;

			System.out.printf("Marke 2 ");

			System.out.printf("yy: %.4f");
			System.out.print("\n");

			int[] j = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21};
			System.out.print("j: ");
			System.out.println(j[5]);
			System.out.print("\n");

			yy = new MWNumericArray(j, MWClassID.INT8);//, MWClassID.DOUBLE);


			/*
			for (int i=0; i<22; i++){
				yy.set(i,y[i]);
			}*/
			//double[] yyA = toDoubleArray(yy);
			System.out.print(yy[21]);
			System.out.print("\n");

			t0 = new MWNumericArray(time, MWClassID.DOUBLE);
			tf = new MWNumericArray(endTime, MWClassID.DOUBLE);

			System.out.println(t0);
			System.out.println(tf);
			System.out.print("\n");

			System.out.printf("Marke 3 ");


			theMagic = new Class1();

			System.out.printf("Marke 4 ");


			yNew = theMagic.stiffODE(1, yy, t0, tf);
			System.out.printf("Marke 5 ");


			System.out.println(yNew[0]);
			System.out.print("\n");
		}
		catch (Exception e)
		{
			System.out.println("Exception: " + e.toString());
		}
	}
	/*
	class RepressilatorODESystem implements ODE {

		// Everything gets upset without mState
		private double[] mState;

		// Let's try to redefine mState in the constructor
		// (useful in case we want some complicated ICs...)
		public RepressilatorODESystem(){
			//this.mState = lamSys.getICs();
		}

		/**
		 * Initial conditions

		@Override
		public double[] getState() {
			return mState;
		}


		  Classic Garcia-Ojalvo repressilator ODEs

		@Override
		public void getRate(double[] y, double[] dy) {

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
			//setSurfaceAreaGrowthRate(lam);
			//System.out.printf("Lambda in 'derivative system' %.4f at t %.2f ", lam, sim.getTime());
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
			//this.mState = dy;
		} */
}
