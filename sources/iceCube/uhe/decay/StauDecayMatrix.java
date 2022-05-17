package iceCube.uhe.decay;

import java.io.*;

import iceCube.uhe.interactions.*;
import iceCube.uhe.particles.*;


/** 
  <pre>
  Matrix of the Energy Transfer by the STAU decays
  The matrix elements are calculated by the methods supplied
  by the Decay class.

  /------------------------------------------------------------\
  logEmin | 0   0  ...................................   EmindN/dEmin  |
  logE1   | 0   0  ............................E1dN/dE1  EmindN/dEmin  |
  .     |                                                            |x binWidth x ln(10)
  .     |                                                            |
  logEmax |EmaxdN/dEmax................................  EmindN/dEmin  |
  \------------------------------------------------------------/

  The bin width and logEmin are defined in the Particle class.

  Actually each matrix element dN/dLogEdecay is calculated by the integral,
  \int dN/dY dY from logY - 0.5xbinWidth to logY + 0.5*binWidth where
  logY = logEdecay - logEtau. This is more accurate way than simply
  calculating dN/dLogEdecay.

  The transfer matrix to the tau
  is acquired by the method getStauToTauDecayMatrix( ),
  while the decay into gravitino
  is obtained by getTauToGravitinoDecayMatrix( ).
  </pre>
  */
public class StauDecayMatrix {

    /** Array for the energy tansfer probability from sTau to Tau.*/
    double[][] tauMtx;
    /** Array for the energy tansfer probability from sTau to gravitino */
    double[][] gravitinoMtx;
    /** Array for the lifetime of tauon considering the Lorentz duration. */
    double[] lifetimeMtx;
    /** Array Dimension */
    int dimension = Particle.getDimensionOfLogEnergyMatrix();
    /** The bin width of the matrix element.*/
    double delta = Particle.getDeltaLogEnergy();
    /** The Particle class. Should be tau. This is checked by the constructor.*/
    Particle p;

    /** mass of tau */
    final static double Mtau   = Particle.particleMasses[2][1]; // tau mass [GeV]
    /** mass ratio square of Mtau/Mstau */
    private double rTau = 1.0e-5;  // initial value. Will be set in the constructor
    /** mass of gravitino */
    double Mgrav = 1.0;  // 10 GeV. You can change it by calling setGravitinoMass()
    /** mass ratio square of Mgrav/Mstau */
    private double rGrav = 4.0e-4;  // initial value. Will be set in the constructor or when change Mgrav


    /** Constructor: Generate the matrix array */
    public StauDecayMatrix(Particle p){
        if(p.getFlavor()==2){ // Requires tau flavor
            this.p = p;
            tauMtx = new double[dimension][dimension];
            gravitinoMtx = new double[dimension][dimension];
            lifetimeMtx = new double[dimension];
            rTau = setMassSquaredRatio(Mtau);
            rGrav = setMassSquaredRatio(Mgrav);
        }else{
            System.err.println("This particle " + 
                    p.particleName(p.getFlavor(), p.getDoublet()) + 
                    " is not TAUs!!");
            System.exit(0);
        }
    }


    private double setMassSquaredRatio(double mass){
        return mass*mass/(p.getMass()*p.getMass());
    }


    /** Set the gravitino mass. Also reculculate the squared ratio.*/
    public void setGravitinoMass(double mass){
        Mgrav = mass;
        rGrav = setMassSquaredRatio(Mgrav);
    }


    /** Return the gravitino mass. */
    public double getGravitinoMass(){
        return Mgrav;
    }


    /** Calculate the decay matrix  from stau to neutrinos, and gravitino*/
    public void setStauDecayMatrix(int iLogE, int jLogE){

        double logEnergy = p.getLogEnergyMinimum()+
            p.getDeltaLogEnergy()*(double )iLogE;
        double energy = Math.pow(10.0,logEnergy);
        p.putEnergy(energy);
        p.putLogEnergy(logEnergy);

        if(isValidIndex(iLogE) && isValidIndex(jLogE)){
            if(jLogE>iLogE){
                tauMtx[iLogE][jLogE]=0.0;
                gravitinoMtx[iLogE][jLogE]=0.0;
            }else{
                double logY = delta*(double )(jLogE-iLogE);
                double logYLow = logY-0.5*delta;
                double logYUp = logY+0.5*delta;
                //System.err.println("logYLow " + logYLow + " logYUp " + logYUp);


                double yLow = Math.pow(10.0,logYLow);

                double yUp = Math.pow(10.0,logYUp);


                // sTau -> tau decay
                double decayElement = 0.0;

                double yUpRange = yUp;
                double yLowRange = yLow;
                if(yUp>=Decay.getYmax(rGrav)) 
                    yUpRange = Decay.getYmax(rGrav);
                if(yLow <= (1.0-Decay.getYmax(rTau))) // lower limit. new for sTau
                    yLowRange = 1.0-Decay.getYmax(rTau);
                if(yLowRange<yUpRange){
                    decayElement = (yUpRange-yLowRange)*
                        Decay.getTauHadronDecayProbToW(yLowRange,rTau+rGrav);
                    tauMtx[iLogE][jLogE] += decayElement;
                    // to tau 
                }
                yUpRange = yUp;
                yLowRange = yLow;
                if(yUp>=Decay.getYmax(rTau))  // uper limit. new for sTau
                    yUpRange = Decay.getYmax(rTau);
                if(yLow<=(1.0-Decay.getYmax(rGrav))) 
                    yLowRange = 1.0-Decay.getYmax(rGrav);
                if(yLowRange<yUpRange){
                    decayElement = (yUpRange-yLowRange)*
                        Decay.getTauHadronDecayProbFromW(yLowRange,rTau+rGrav);
                    gravitinoMtx[iLogE][jLogE] = decayElement;
                }else{
                    gravitinoMtx[iLogE][jLogE] = 0.0;
                }

            }
        }
    }

    /** Calculate the life time matrix considering the Lorentz duration */
    public void setLifeTimeMatrix(int iLogE){
        double logEnergy = p.getLogEnergyMinimum()+
            p.getDeltaLogEnergy()*(double )iLogE;
        double energy = Math.pow(10.0,logEnergy);
        p.putEnergy(energy);
        p.putLogEnergy(logEnergy);
        if(isValidIndex(iLogE)){
            //System.out.println("LifeTime:" + p.getLifeTime());
            lifetimeMtx[iLogE] = p.getLifeTime()*p.getEnergy( )/p.getMass( );
        }
    }

    /** get the element of the decay matrix  of stau-> tau */
    public double getStauToTauDecayMatrix(int iLogE, int jLogE){
        if(isValidIndex(iLogE) && isValidIndex(jLogE)){
            return tauMtx[iLogE][jLogE];
        }
        else{
            return 0.0;
        }
    }

    /** get the element of the decay matrix  of stau-> gravitino */
    public double getStauToGravitinoDecayMatrix(int iLogE, int jLogE){
        if(isValidIndex(iLogE) && isValidIndex(jLogE)){
            return gravitinoMtx[iLogE][jLogE];
        }
        else{
            return 0.0;
        }
    }

    /** get the element of the LifeTime matrix */
    public double getLifeTimeMatrix(int iLogE){
        if(isValidIndex(iLogE)){
            return lifetimeMtx[iLogE];
        }
        else{
            return 0.0;
        }
    }


    /** Checking the energy index */
    public boolean isValidIndex(int iLogE){
        if(0<=iLogE && iLogE<dimension) return true;
        else return false;
    }

}
