package BSimReplication;

import java.awt.Color;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.util.Vector;

import javax.vecmath.Vector3d;

import bsim.export.BSimLogger;
//import org.opensourcephysics.numerics.ODE;
//import processing.core.PGraphics3D;
import bsim.BSim;
import bsim.BSimTicker;
//import bsim.draw.BSimP3DDrawer;
import bsim.particle.BSimBacterium;
import bsim.*;


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
		BSim sim = new BSim();
		sim.setDt(60);
		sim.setSimulationTime(36000);
						
		/*********************************************************
		 * Set up the bacteria
		 */

		final Vector<BSimBacterium> bacteria = new Vector<BSimBacterium>();
		final Vector<BSimBacterium> children = new Vector<BSimBacterium>();
		while(bacteria.size() < 1) {		
			BSimBacterium b = new BSimBacterium(sim, new Vector3d(Math.random()*sim.getBound().x, Math.random()*sim.getBound().y, Math.random()*sim.getBound().z));
			b.setRadius();
			b.setSurfaceAreaGrowthRate();
			b.setChildList(children);
			bacteria.add(b);		
		}

		/*********************************************************
		 * Set up the ticker
		 */

		sim.setTicker(new BSimTicker() {
			@Override
			public void tick() {
				for(BSimBacterium b : bacteria) {
					b.action();		
					b.updatePosition();
					//System.out.printf("first entry %.2f ", b.getY(0));
				}
				bacteria.addAll(children);
				children.clear();
				System.out.printf("Time: %s ", sim.getTime());
				System.out.printf("number %d ", bacteria.size());
				System.out.printf("\n");
			}		
		});

		/*********************************************************
		 * Set up the drawer


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
		sim.preview();*/

		BSimLogger CellNumber = new BSimLogger(sim, "cellnumbers.csv") {

			@Override
			public void before() {
				super.before();
				// Write a header containing the names of variables we will be exporting
				write("time, ,cellnumber");
			}

			@Override
			public void during() {
				// Write the time and number into file
				//BSimBacterium first = bacteria.firstElement();
				//double[] Y = first.y;
				//System.out.printf("Y in Replication: %d ", Y[0]);
				write(sim.getFormattedTime()+","+bacteria.size());
			}
		};
		CellNumber.setDt(60);
		sim.addExporter(CellNumber);

        BSimLogger tars = new BSimLogger(sim, "cells_trajectories.csv") {

        	DecimalFormat formatter = new DecimalFormat("##########.##################", DecimalFormatSymbols.getInstance( Locale.ENGLISH ));
            @Override
            public void before() {
                super.before();

                String buffer = "time, ";

                if(bacteria.size() > 0) {
                    for (int j = 0; j < 22; j++) {
                        buffer += ",tra(" + j + ")";
                    }
                }

                write(buffer);
            }

            @Override
            public void during() {
                String buffer = "" + sim.getFormattedTime();

                for(BSimBacterium b : bacteria) {
                    for(int j = 0; j < 22; j++) {
                        buffer += "," + formatter.format(b.getY(j));
                    }
                }

                write(buffer);
            }
        };
        tars.setDt(60);
        sim.addExporter(tars);

		BSimLogger rad = new BSimLogger(sim, "cells_radii.csv") {

			DecimalFormat formatter = new DecimalFormat("###.##################", DecimalFormatSymbols.getInstance( Locale.ENGLISH ));
			@Override
			public void before() {
				super.before();

				String buffer = "time, , SurfaceRadius, Current(Model)Radius";

				write(buffer);
			}

			@Override
			public void during() {
				String buffer = "" + sim.getFormattedTime();

				for(BSimBacterium b : bacteria) {
					buffer += "," + formatter.format(b.getRadius()) + "," + formatter.format(b.getCurrentRadius());
				}

				write(buffer);
			}
		};
		rad.setDt(60);
		sim.addExporter(rad);

        BSimLogger lambda = new BSimLogger(sim, "cells_lam.csv") {

            DecimalFormat formatter = new DecimalFormat("###.####################", DecimalFormatSymbols.getInstance( Locale.ENGLISH ));
            @Override
            public void before() {
                super.before();
                String buffer = "time, , lambda";
                write(buffer);
            }

            @Override
            public void during() {
                String buffer = "" + sim.getFormattedTime();
                for(BSimBacterium b : bacteria) {
                    buffer += "," + formatter.format(b.getLambda());
                }
                write(buffer);
            }
        };
        lambda.setDt(60);
        sim.addExporter(lambda);

		sim.export();
	}

}
