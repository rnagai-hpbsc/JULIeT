package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Make the energy transfer matrix of PhotoNuclear Distribution Interactions */
public class MakeDynamicPhotoNuclearMtx {

    public static void main(String[] args) throws IOException {

        String fileName = null;
        int flavor = 1;
        int doublet =1;        // Charged Lepton
        int material = 0;      // Rock (1) Ice (0)
        double mass = 1.0;
        double energy = 1.0e9; // 1EeV= 10^9 GeV
        double epsilon = 0.01;
        double energyCut = 1.0e4; // 10000 GeV


        if(args.length!=4){
            System.out.println("Usage: MakeNeutrinoChargeMtx file-name flavor material mass");
            System.exit(0);
        }else{
            fileName = args[0];
            flavor = Integer.valueOf(args[1]).intValue();
            material = Integer.valueOf(args[2]).intValue();
            mass  = Double.valueOf(args[3]).doubleValue();
        }


        // Generate the particle class.
        DynamicParticle lepton = 
            new DynamicParticle(flavor, doublet, energy); // Charged Leptons
        lepton.setMass(mass);
        System.err.format("mass is set to %e [GeV]\n",lepton.getMass());
        lepton.setIsNP(1);
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

        //Generate object of the Photo-Nuclear interaction.
        PhotoNuclear photoNucl = new PhotoNuclear(lepton, s);
        System.err.println(photoNucl.interactionName( ));

        energyCut = 1.0;//Math.pow(10.0,lepton.getLogEnergyMinimum( )-3.0); //GeV
        photoNucl.setEnergyCut(energyCut);
        photoNucl.setRoundOffError(0);

        //Generate object of the Interaction Matrix
        InteractionsMatrix photoNuclMtx = new InteractionsMatrix(photoNucl);


        numRecipes.Integration.setRelativeAccuracy(epsilon);
        // integration accuracy for save CPUtime

        int iLogE;
        double raw = 0.0;
        for(iLogE=0;iLogE<lepton.getDimensionOfLogEnergyMatrix();iLogE++){
            photoNucl.setIncidentParticleEnergy(iLogE);
            if(iLogE%50==0) System.err.println("The Incident energy " + 
                    photoNucl.getIncidentParticleEnergy() + " GeV");

            // Total Cross Section
            photoNuclMtx.setSigmaMatrix(iLogE);
            //System.err.println("  Total cross section done");


            // Transfer Matrix
            int jLogE;
            for(jLogE=0;jLogE<lepton.getDimensionOfLogEnergyMatrix();jLogE++){
                double ret = photoNuclMtx.setTransferMatrixHighMass(iLogE,jLogE);
                if (jLogE==iLogE) raw = ret; 
                //if(jLogE%100==0) System.err.println("  Transfer Matrix " + jLogE);
            }

            // check sigma matrix
            photoNuclMtx.verifySigmaMatrix(iLogE);

	    double allsum = photoNuclMtx.getSigmaMatrix(iLogE);
	    double parsum = photoNuclMtx.getLeptonTransferMatrix(iLogE,iLogE);
	    if (allsum<parsum){
	    	System.err.println("Warning! Total Xsec less than survival at " + iLogE + ": "  + allsum + " " + parsum);
	    }else{
            if (raw > allsum) {
	    	    System.err.println("Total Xsec, survival" + iLogE + ": "  + allsum + " " + parsum + " *** raw: " + raw + " ("+String.format("%.2f",(raw/allsum*100.))+"%)");
            }else {
	    	    System.err.println("Total Xsec, survival" + iLogE + ": "  + allsum + " " + parsum + " dif: " + (allsum-parsum)+" ("+ String.format("%.2f",((allsum-parsum)/allsum*100.))  +"%)");
            }
	    }

        }

        FileOutputStream out = new FileOutputStream(fileName);
        InteractionsMatrixOutput.outputInteractionsMatrix(photoNuclMtx, out);
        out.close( );

    }
}
