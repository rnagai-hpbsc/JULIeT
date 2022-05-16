package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;
import java.util.Scanner;

/** Make the differential cross section table of Pair Creation Interactions
  for a given DynamicParticle, that allows to input a different mass and lifetime 
  The generated data will be used for PairCreationFit object.*/

public class MakeDynamicPairCreationTable {

    public static void main(String[] args) throws IOException {

        String fileName = null;
        int flavor = 1;
        int producedFlavor = 0;     // Recoiled particle flavor
        int doublet =1;        // Charged Lepton
        int material = 0;      // Ice-> 0, Rock->1
        double energy = 1.0e9; // 1EeV= 10^9 GeV
        double epsilon = 0.01;
        double mass = 1.0;      // initial value [GeV]
        double ln10 = Math.log(10.0);
        double[ ] logEArray = new double[140];
        double[ ] logDyArray = new double[35];
        double[ ][ ] yDsigmaArray = new double[140][35];



        if(args.length!=5){
            System.out.println("Usage: MakeDynamicPairCreationTable file-name flavor producedFlavor material mass");
            System.exit(0);
        }else{
            fileName = args[0];
            flavor = Integer.valueOf(args[1]).intValue();
            producedFlavor = Integer.valueOf(args[2]).intValue();
            material = Integer.valueOf(args[3]).intValue();
            mass  = Double.valueOf(args[4]).doubleValue();
        }


        // Generate the particle class.
        DynamicParticle lepton = 
            new DynamicParticle(flavor, doublet, energy); // Charged Leptons
        lepton.setMass(mass);
        System.err.format("mass is set to %e [GeV]\n",lepton.getMass());


        System.err.println("The Particle Name is " + 
                lepton.particleName(lepton.getFlavor(), lepton.getDoublet()));

        // Generate the ParticlePoint class.
        ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,material);

        for(int i=0;i<s.NumberOfSpecies[material];i++){
            System.out.println("Charge " + s.getCharge(i));
            System.out.println("Atomic Number " + s.getAtomicNumber(i));
        }

        //Generate object of the PairCreation interaction.
        PairCreation pairC = new PairCreation(lepton, s, producedFlavor);
        pairC.setRoundOffError(0);
        System.err.println(pairC.interactionName( ));
        System.err.format("the incoming particle's mass in pair creation %e [GeV]\n",pairC.mass);

        numRecipes.Integration.setRelativeAccuracy(epsilon);
        // integration accuracy for save CPUtime

        pairC.setIncidentParticleEnergy(0);
        //double Ymax = pairC.getYmaxCharge(0)*0.99999*1.e-4*pairC.getProducedMass()/0.511e-3;
	double Ymax = pairC.getYmaxCharge(0)*0.99999;
	double Ymin = pairC.getYmin();
	double LogMidY = (Math.log10(Ymax)+Math.log10(Ymin))/2.;
	double dULogY = (Math.log10(Ymax) - LogMidY)/5.;
	double dLLogY = (LogMidY - Math.log10(Ymin))/30.;
	System.err.println(Ymax + " " + pairC.getProducedMass()/0.511e-3 );

        int iLogE;
        for(iLogE=0;iLogE<lepton.getDimensionOfLogEnergyMatrix();iLogE+=5){
            pairC.setIncidentParticleEnergy(iLogE);
            energy = pairC.getIncidentParticleEnergy();
            System.out.println("The Incident energy " + iLogE + ": " + 
                    energy + " GeV" );
            logEArray[iLogE/5] = 
                Math.log(energy)/ln10;

            // differential cross sections 
            int jLogY; 
            int y_min_index = 0; 
            for(jLogY=0;jLogY<35;jLogY++){
                //double Y = Math.pow(10.0,-0.2*(double )jLogY)*Ymax;
		double Y = jLogY < 5 ? Math.pow(10.,LogMidY+dULogY*(5.-jLogY)) : Math.pow(10.,LogMidY-dLLogY*(jLogY-5.));
                //double Y = Math.pow(10.0,0.2*(35.-(double )jLogY))*Ymin;
                logDyArray[jLogY] = Math.log(Y)/ln10;
                double dSigma;
                if(y_min_index == 0) {
                    dSigma = pairC.getDSigmaDyHighMass(Y);
                    System.err.println("Calculated dSigma: " + dSigma);
                }
                else dSigma = 1.0e-50;
                if(dSigma<=0.0){
                    dSigma = 1.0e-50;
                    y_min_index = 1;
                }
                yDsigmaArray[iLogE/5][jLogY]= 
                    Math.log(dSigma*Y/(energy*1.0e-32))/ln10;
                System.err.println("  Transfer Matrix " + jLogY + " Y " + Y + " " +
                        yDsigmaArray[iLogE/5][jLogY]);
                //System.err.println("  Transfer Matrix "  + Y + " " + dSigma + " " + yDsigmaArray[iLogE/5][jLogY]); 

            }

        }

        DataOutputStream out =  
            new DataOutputStream(new FileOutputStream(fileName));

        int jLogY;
        for(jLogY=0;jLogY<35;jLogY++) out.writeDouble(logDyArray[jLogY]);
        for(iLogE=0;iLogE<lepton.getDimensionOfLogEnergyMatrix();iLogE+=5){
            out.writeDouble(logEArray[iLogE/5]);
            for(jLogY=0;jLogY<35;jLogY++){
                out.writeDouble(yDsigmaArray[iLogE/5][jLogY]);
                //System.out.println(yDsigmaArray[iLogE/5][jLogY]);
            }
        }

        out.close( );

    }
}

