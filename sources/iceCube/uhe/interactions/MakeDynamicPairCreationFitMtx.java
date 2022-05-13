package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Make the energy transfer matrix of PairCreationFit Distribution Interactions */
public class MakeDynamicPairCreationFitMtx {

    public static void main(String[] args) throws IOException {

        String fileName = null;
        int flavor = 1;          // Incident particle's flavor
        int producedFlavor = 0;  // Produced particle's flavor
        int doublet =1;        // Charged Lepton
        int material = 0;      // Ice->0, Rock->1
        double mass = 1.0;     // mass [GeV]
        double energy = 1.0e9; // 1EeV= 10^9 GeV
        double epsilon = 1.0e-3;
        double energyCut = 1.0e4; // 10000 GeV
        double ln10 = Math.log(10.0);


        if(args.length!=5){
            System.out.println("Usage: MakePairCreationFitMtx file-name flavor producedFlavor material mass");
            System.exit(0);
        }else{
            fileName = args[0];
            flavor = Integer.valueOf(args[1]).intValue();
            producedFlavor = Integer.valueOf(args[2]).intValue();
            material = Integer.valueOf(args[3]).intValue();
            mass = Double.valueOf(args[4]).doubleValue();
        }


        // Generate the particle class.
        DynamicParticle lepton = 
            new DynamicParticle(flavor, doublet, energy); // Charged Leptons
        lepton.setMass(mass);
        lepton.setIsNP(1);
        System.err.format("mass is set to %e [GeV]\n",lepton.getMass());
        if (lepton.getIsNP()>0){
            System.out.println("New Physics setting!");
        }


        System.err.println("The Particle Name is " + 
                lepton.particleName(lepton.getFlavor(), lepton.getDoublet()));

        // Generate the ParticlePoint class.
        ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,material);

        for(int i=0;i<s.NumberOfSpecies[material];i++){
            System.out.println("Charge " + s.getCharge(i));
            System.out.println("Atomic Number " + s.getAtomicNumber(i));
        }

        //Generate object of the PairCreationFit interaction.
        PairCreationFit pairC = new PairCreationFit(lepton, s, producedFlavor);
        System.err.println(pairC.interactionName( ));
        pairC.setRoundOffError(0);

        //Generate object of the Interaction Matrix
        InteractionsMatrix pairCMtx = new InteractionsMatrix(pairC);


        numRecipes.Integration.setRelativeAccuracy(epsilon);
        // integration accuracy for save CPUtime

        int iLogE;
        double raw = 0.0;
        for(iLogE=0;iLogE<lepton.getDimensionOfLogEnergyMatrix();iLogE++){
            pairC.setIncidentParticleEnergy(iLogE);
            energy = pairC.getIncidentParticleEnergy();
            if(iLogE%50==0) System.err.println("The Incident energy " + energy + " GeV");

            // Total Cross Section
            pairCMtx.setSigmaMatrix(iLogE);
            //System.err.format("  Total cross section done %e cm2\n",pairCMtx.getSigmaMatrix(iLogE));


            // Transfer Matrix
            int jLogE;
            for(jLogE=0;jLogE<lepton.getDimensionOfLogEnergyMatrix();jLogE++){
                double ret = pairCMtx.setTransferMatrixHighMass(iLogE,jLogE);
                if (jLogE==iLogE) raw = ret;
                //if(jLogE%200==0) System.err.format("  Transfer Matrix (%d) %e\n",
                //				   jLogE, pairCMtx.getTransferMatrix(iLogE,jLogE));
            }

            // check sigma matrix
            pairCMtx.verifySigmaMatrix(iLogE);

	    double allsum = pairCMtx.getSigmaMatrix(iLogE);
	    double parsum = pairCMtx.getLeptonTransferMatrix(iLogE,iLogE);
	    if (allsum<parsum){
	    	System.err.println("Warning! Total Xsec less than survival at " + iLogE + ": "  + allsum + " " + parsum);
	    }else{
            if (raw > allsum) {
	    	    System.err.println("Total Xsec, survival" + iLogE + ": "  + allsum + " " + parsum + " *** raw: " + raw + " ("+String.format("%.2f",(raw/allsum*100.))+"%)");
            }else {
	    	    System.err.println("Total Xsec, survival" + iLogE + ": "  + allsum + " " + parsum + " dif: " + (allsum-parsum) + " ("+ String.format("%.2f",((allsum-parsum)/allsum*100.))  +"%)");
            }
	    }
        }

        FileOutputStream out = new FileOutputStream(fileName);
        InteractionsMatrixOutput.outputInteractionsMatrix(pairCMtx, out);
        out.close( );

    }
}
