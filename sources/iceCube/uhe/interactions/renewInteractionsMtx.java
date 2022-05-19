package iceCube.uhe.interactions;

import iceCube.uhe.interactions.*;
import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import java.io.*;

public class renewInteractionsMtx {

    public static void main(String[] args) throws IOException {
        String intName = null;
        int flavor = 1;
        int producedFlavor = 0;
        int doublet = 1;
        int material = 0;
        double mass = 1.0;
        double energy = 1.0e9;
        double epsilon = 1.0e-3;
        double energyCut = 1.0e4;
        
        if (args.length!=5){
            System.out.println("Usage: renewInteractionMtx interaction-name flavor producedFlavor material mass");
            System.exit(0);
        } else {
            intName = args[0];
            flavor = Integer.valueOf(args[1]).intValue();
            producedFlavor = Integer.valueOf(args[2]).intValue();
            material = Integer.valueOf(args[3]).intValue();
            mass = Double.valueOf(args[4]).doubleValue();
        }

        String[] storeDir = {"iceCube/uhe/interactions/ice/new/","iceCube/uhe/interactions/rock/new/"};

        DynamicParticle lepton = new DynamicParticle(flavor, doublet, energy);
        lepton.setMass(mass);
        lepton.setIsNP(1);
        String renewFilename =  "stau" + getFileName(intName, producedFlavor, mass);
        System.out.println("Renew file: " + renewFilename);

        ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0, material);

        InteractionsMatrix intMtx = selectIntMatrix(intName, lepton, s, producedFlavor);
        fillMatrixElements(intMtx, storeDir[material] + renewFilename + "_sigma.txt", 0);
        for (int i=0; i<700; i++) { 
            double allsum = intMtx.getSigmaMatrix(i);
            System.out.print(allsum + ((i+1)==700?"\n":","));
        }
        fillMatrixElements(intMtx, storeDir[material] + renewFilename + "_inela.txt", 1);
        for (int i=0; i<700; i++) {
            double inela = intMtx.getInelasticityMatrix(i);
            System.out.print(inela + ((i+1)==700?"\n":","));
        }
        fillMatrixElements(intMtx, storeDir[material] + renewFilename + "_trans.txt", 2);
        fillMatrixElements(intMtx, storeDir[material] + renewFilename + "_surviv.txt", 3);

        FileOutputStream out = new FileOutputStream(storeDir[material] + renewFilename);
        InteractionsMatrixOutput.outputInteractionsMatrix(intMtx, out);
        out.close( );

    }

    public static void fillMatrixElements(InteractionsMatrix intMtx, String filename, int func) throws IOException {
        File file = new File(filename);
        System.out.println("Read " + filename + "..."); 
        BufferedReader reader = null;

        try {
            reader = new BufferedReader(new FileReader(file));

            String text;
            int j = 0;
            while ((text = reader.readLine()) != null) {
                String[] nums = text.split(",");
                for (int i=0; i<nums.length; i++) {
                    double num = Double.parseDouble(nums[i]);
                    selectFillMatrix(intMtx, func, i, j, num);
                }
                j++; 
            }
        } catch (FileNotFoundException e) {
            System.out.println("Not found");
        } catch (IOException e) {
            System.out.println("Falied");
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }

    public static void selectFillMatrix(InteractionsMatrix intMtx, int func, int i, int j, double value) {
        switch (func) {
            case 0: intMtx.fillSigmaMatrix(i, value); break;
            case 1: intMtx.fillInelasticityMatrix(i, value); break;
            case 2: intMtx.fillTransferMatrix(i, j, value); break; 
            case 3: intMtx.fillLeptonTransferMatrix(i, j, value); break;
        }
    }

    public static InteractionsMatrix selectIntMatrix(String mtxName, DynamicParticle lepton, ParticlePoint s, int producedFlavor) throws IOException {
        switch (mtxName) {
            case "pairC" : 
                PairCreationFit pairC = new PairCreationFit(lepton, s, producedFlavor);
                pairC.setRoundOffError(0);
                InteractionsMatrix pairCMtx = new InteractionsMatrix(pairC);
                return pairCMtx;
            case "brems" : 
                Bremsstrahlung brems = new Bremsstrahlung(lepton, s);
                InteractionsMatrix bremsMtx = new InteractionsMatrix(brems);
                return bremsMtx;
            case "photonuc" : 
                PhotoNuclear photonuc = new PhotoNuclear(lepton, s);
                InteractionsMatrix photonucMtx = new InteractionsMatrix(photonuc);
                return photonucMtx;
        }
        return null;
    }

    public static String getFileName(String intName, int producedFlavor, double mass) {
        String filename = "";
        switch (intName) {
            case "pairC" : 
                switch (producedFlavor){
                    case 0: filename += "ToePairCreation"; break;
                    case 1: filename += "TomuPairCreation"; break;
                    case 2: filename += "TotauPairCreation"; break;
                }
                break;
            case "brems": filename += "Bremsstrahlung"; break;
            case "photonuc": filename += "PhotoNuclear"; break;
        }
        filename += "Mtx" + (int)mass + "GeV";
        return filename;
    }

}
