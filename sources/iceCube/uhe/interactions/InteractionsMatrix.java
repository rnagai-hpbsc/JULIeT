package iceCube.uhe.interactions;

import java.io.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;


/** 
<pre>
    Matrix of the Energy Transfer by particle Interactions.
    The matrix elements are calculated by the methods supplied
    by the Interaction class.

        /------------------------------------------------------------\
logEmin | 0   0  ...................................   EmindS/dEmin  |
logE1   | 0   0  ............................E1dS/dE1  EmindS/dEmin  |
  .     |                                                            |x binWidth x ln(10)
  .     |                                                            |
logEmax |EmaxdS/dEmax................................  EmindS/dEmin  |
        \------------------------------------------------------------/

The bin width and logEmin are defined in the Particle class.

Actually each matrix element dSigma/dLogE is calculated by the integral,
\int dSigma/dY dY from logY - 0.5xbinWidth to logY + 0.5*binWidth where
logY = logEproduced - logEincoming. This is more accurate way than simply
calculating dSigma/dLogE.

The transfer matrix of the recolied lepton energy (1-y)*Eincoming 
is acquired by the method getLeptonTransferMatrix( )
while the energy of the transfered matrix for the counter current 
y*Eincoming is obtained by getTransferMatrix( ).
</pre>
*/

public class InteractionsMatrix implements Serializable {

    Interactions interactions;
    private static final long serialVersionUID = -3532858511737660878L;

    double[][] transferMtx;
    double[][] transferAMtx;
    double[] sigmaMtx;
    double[] inelasticityMtx;
    double[][] rawTAMtx; // added 

    /** Array Dimension */
    int dimension = Particle.getDimensionOfLogEnergyMatrix();
    /** Bin width of the matrix element. */
    double delta = Particle.getDeltaLogEnergy();

    /** Constructor: Generate the matrix array */
    public InteractionsMatrix(Interactions interactions){
	this.interactions = interactions;
	transferMtx = new double[dimension][dimension];
	transferAMtx = new double[dimension][dimension];
	sigmaMtx = new double[dimension];
	inelasticityMtx = new double[dimension];
    rawTAMtx = new double[dimension][dimension];
    }




    /** Calculate the transfer matrix */
    public void setTransferMatrix(int iLogE, int jLogE){

	interactions.setIncidentParticleEnergy(iLogE);/** Setup the primary E*/

	if(isValidIndex(iLogE) && isValidIndex(jLogE)){
	    if(jLogE>iLogE){
		transferMtx[iLogE][jLogE]=0.0;
		transferAMtx[iLogE][jLogE]=0.0;
	    }
	    else{
		double logY = delta*(double )(jLogE-iLogE);
		double logYLow = logY-0.5*delta;
		double logYUp = logY+0.5*delta;
		//System.err.println("logYLow " + logYLow + " logYUp " + logYUp);

		double yLow = Math.pow(10.0,logYLow);
		double yLowRange = yLow;
		if(yLow<=interactions.getYmin()+2.0*interactions.roundOffError)  
		    yLowRange = interactions.getYmin( ) + 2.0*interactions.roundOffError;

		double yUp = Math.pow(10.0,logYUp);
		double yUpRange = yUp;
		if(yUp>=interactions.getYmax()-2.0*interactions.roundOffError)   
		    yUpRange = interactions.getYmax( ) - 2.0*interactions.roundOffError;


		if(yLowRange<yUpRange){
		    transferMtx[iLogE][jLogE]=
			interactions.integralDSigmaDy(yLowRange,yUpRange);
		}else{
		    transferMtx[iLogE][jLogE]=0.0;
		}

		yLowRange = yLow; yUpRange = yUp;
		if(yLow<=(1.0-interactions.getYmax()+2.0*interactions.roundOffError)) 
		    yLowRange = 1.0-interactions.getYmax( )+2.0*interactions.roundOffError;
		if(yUp>=(1.0-interactions.getYmin()-2.0*interactions.roundOffError))  
		    yUpRange = 1.0-interactions.getYmin( )-2.0*interactions.roundOffError;

		if(yLowRange<yUpRange){
		    transferAMtx[iLogE][jLogE]=
			interactions.integralDSigmaDz(yLowRange,yUpRange);
		}else{
		    transferAMtx[iLogE][jLogE]=0.0;
		}

	    }
	}
    }

    /** Calculate the transfer matrix for High mass particles **/
    public double setTransferMatrixHighMass(int iLogE, int jLogE){

        interactions.setIncidentParticleEnergy(iLogE);/** Setup the primary E*/

        if(isValidIndex(iLogE) && isValidIndex(jLogE)){
            if(jLogE>iLogE){
                transferMtx[iLogE][jLogE] = 0.0;
                transferAMtx[iLogE][jLogE] = 0.0;
            }
            else {
                double logY = delta*(double)(jLogE-iLogE);
                double logYLow = logY - 0.5*delta;
                double logYUp  = logY + 0.5*delta;

                
                double yLow = Math.pow(10.0,logYLow);
                double yLowRange = yLow;
                if(yLow<=interactions.getYmin())
                    yLowRange = interactions.getYmin( );

                double yUp = Math.pow(10.0,logYUp);
                double yUpRange = yUp;
                if(yUp>=interactions.getYmax()-2.0*interactions.roundOffError)
                    yUpRange = interactions.getYmax( ) - 2.0*interactions.roundOffError;

                if(yLowRange<yUpRange){
                    transferMtx[iLogE][jLogE]=
                        interactions.integralDSigmaDy(yLowRange,yUpRange);
                }else{
                    transferMtx[iLogE][jLogE] = 0.0;
                }

                /*yLowRange = yLow; yUpRange = yUp;
                if(yLow<=(1.0-interactions.getYmax()+2.0*interactions.roundOffError))
                    yLowRange = 1.0 - interactions.getYmax( ) + 2.0*interactions.roundOffError;
                if(yUp>=(1.0-interactions.getYmin()))
                    yUpRange = 1.0 - interactions.getYmin();

                if(yLowRange<yUpRange){
                    transferAMtx[iLogE][jLogE] = 
                        interactions.integralDSigmaDz(yLowRange,yUpRange);*/
                double zLowRange = yLow; 
                double zUpRange = yUp;
                if(zLowRange<=interactions.getZmin())
                    zLowRange = interactions.getZmin();
                if(zUpRange>=interactions.getZmax())
                    zUpRange = interactions.getZmax();
                if(1.-zUpRange<1.-zLowRange){
                    transferAMtx[iLogE][jLogE] = 
                        interactions.integralDSigmaDy(1.-zUpRange,1.-zLowRange);
                }else{
                    transferAMtx[iLogE][jLogE]=0.0;
                }

                rawTAMtx[iLogE][jLogE] = transferAMtx[iLogE][jLogE];
                double rawTransferAMtx = transferAMtx[iLogE][jLogE];

                //if (transferAMtx[iLogE][jLogE] > getSigmaMatrix(iLogE))
                //    transferAMtx[iLogE][jLogE] = getSigmaMatrix(iLogE)*0.99999;

                return rawTransferAMtx;
            }
            return 0.0;
        }
        return 0.0;
    }

    /** Get the element of the transfter matrix  jlogE ~ log10(y*E) */
    public double getTransferMatrix(int iLogE, int jLogE){
	if(isValidIndex(iLogE) && isValidIndex(jLogE)){
	    return transferMtx[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /** Get the element of the transfter matrix  jlogE ~ log10((1-y)*E) */
    public double getLeptonTransferMatrix(int iLogE, int jLogE){
	if(isValidIndex(iLogE) && isValidIndex(jLogE)){
	    return transferAMtx[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }



    /** Calculate the total cross section matrix */
    public void setSigmaMatrix(int iLogE){
	    
	interactions.setIncidentParticleEnergy(iLogE);/** Setup the primary E*/
	if(isValidIndex(iLogE)){
	    sigmaMtx[iLogE] = interactions.getSigma( );
	    inelasticityMtx[iLogE] = 
	    interactions.getYDSigmaDy(interactions.getYmin()+interactions.roundOffError,
				      interactions.getYmax()-interactions.roundOffError);
	}
    }

    public void verifySigmaMatrix(int iLogE){
        int flag = 0;
        for (int jLogE=0;jLogE<dimension;jLogE++){
            if(transferAMtx[iLogE][jLogE]>sigmaMtx[iLogE]) flag += 1;
        }
        if (flag>0){
            sigmaMtx[iLogE] = 0.0;
            for (int jLogE=0;jLogE<dimension;jLogE++){
                sigmaMtx[iLogE] += transferAMtx[iLogE][jLogE];
            }
        }
    }

    public void fillSigmaMatrix(int iLogE, double value){
        sigmaMtx[iLogE] = value;
    }

    public void fillInelasticityMatrix(int iLogE, double value){
        inelasticityMtx[iLogE] = value;
    }

    public void fillTransferMatrix(int iLogE, int jLogE, double value){
        transferMtx[iLogE][jLogE] = value;
    }

    public void fillLeptonTransferMatrix(int iLogE, int jLogE, double value){
        transferAMtx[iLogE][jLogE] = value;
    }


    /** Get the element of the total cross section matrix */
    public double getSigmaMatrix(int iLogE){
	if(isValidIndex(iLogE)){
	    return sigmaMtx[iLogE];
	}
	else{
	    return 0.0;
	}
    }

    /** Get the element of the inelastisity matrix */
    public double getInelasticityMatrix(int iLogE){
	if(isValidIndex(iLogE)){
	    return inelasticityMtx[iLogE];
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


    /** Get the flavor of the particle propagating */
    public int getFlavor(){
	int flavor = interactions.p.getFlavor();
	return flavor;
    }

    /** Get the doublet of the particle propagating */
    public int getDoublet(){
	int doublet = interactions.p.getDoublet();
	return doublet;
    }

    /** Get the flavor of the produced particle */
    public int getProducedFlavor(){
	int producedFlavor = interactions.producedFlavor;
	return producedFlavor;
    }

}


