package iceCube.uhe.particles;

import iceCube.uhe.particles.*;
import iceCube.uhe.geometry.*;
import geometry.*;

import java.io.*;
import java.util.*;

/**
  This particle class inherited from Particle.java describes
  a particle of dyamical mass and lifetime.
  <P>
  Written originaly for the IceCube sTau analysis
  by S.Yoshida 2020/03/12

*/

public class DynamicParticle extends Particle implements Serializable {

    private static final long serialVersionUID = 20210313L;
    double mass = 1.0; // [GeV], initial value. hiding the mass variable in Particle class
    double lifetime = Double.POSITIVE_INFINITY; // lifetime. highding the same variable in Particle class
    int isNP = 1;

    /** Constructor.
      <pre>
      flavor ... flavor valuable
      doublet ... doublet valuable
      energy ... initial Energy [GeV]
      </pre>
      */
    public DynamicParticle(int flavor, int doublet, double energy){
        super(flavor,doublet,energy);
    }

    /** Constructor.
      <pre>
      flabor ... flavor valuable
      doblet ... doublet valuable
      </pre>
      */
    public DynamicParticle(int flavor, int doublet){
        super(flavor,doublet);
    }

    /**
      setting mass [GeV]
      */

    public void setMass(double mass){
        this.mass = mass;
    }

    /** Return mass. Overriding the method of Particle class */
    public double getMass(){
        return mass;
    }

    /**
      setting lifetime [sec]
      */
    public void setLifeTime(double lifetime){
        System.out.println("set lifetime to: " + lifetime);
        this.lifetime = lifetime;
        System.out.println("set lifetime to: " + this.lifetime);
    }

    public double getLifeTime(){
        return lifetime;
    }

}
