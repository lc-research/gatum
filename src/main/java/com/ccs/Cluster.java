/*
 *  Copyright (c) 2018, LC-Research. (http://www.lc-research.com)
 *
 *  LC-Research licenses this file to you under the Apache License V 2.0.
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
 *  Unless required by applicable law or agreed to in writing, software distributed under the
 *  License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 *  CONDITIONS OF ANY KIND, either express or implied.  See the License for the
 *  specific language governing permissions and limitations under the License.
 */

package com.ccs;

import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.geometry.GeometryUtil;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import javax.vecmath.Point3d;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;

import static com.ccs.Constants.*;

/**
 * <h1>Cluster!</h1>
 * This class models the properties and associated functions of a molecule
 *
 *
 * @author  M. Dias
 * @version 1.0
 */
public class Cluster implements Cloneable{

    private double[] rolj,eolj,rhs;
    private IAtomContainer Molecule;
    private IAtomContainer originalMolecule;
    private double pot,dpotx,dpoty,dpotz,dmax;

    private double NTHETA=3.6297469670513871E-308;
    private double NPHI=6.9533558064128858E-310;
    private double NGAMMA=6.9533558064128858E-310;
    private final double BOND=1.0976E-10;
    private double EOGAS=0.06900;
    private double DIPOL=(0.204956E-30/(2*4*PI*XEO))*XE*XE;
    private final int IBST_MAX=500;
    private final double ROGAS =3.6600;
    private final double CONVE=4.2*0.01036427;
    private final double CONVR=0.890898718;
    private final double PC=-0.4825;
    private final double PC_CENTER=-(PC);
    private final double DIPOL_D=1.710E-30;

    public IAtomContainer getMolecule()
    {
        return Molecule;
    }
    public IAtomContainer getOriginalMolecule()
    {
        return originalMolecule;
    }
    public double getDmax() {
        return this.dmax = dmax;
    }

    public double getPot() {
        return pot;
    }

    public double getDpotx() {
        return dpotx;
    }

    public double getDpoty() {
        return dpoty;
    }

    public double getDpotz() {
        return dpotz;
    }

    public Cluster(IAtomContainer Mol)

    {


        try {
            Molecule=(IAtomContainer)Mol.clone();
        }
        catch (CloneNotSupportedException e)
        {
            //logger.debug(e);
        }
        int NumberOfAtoms=Molecule.getAtomCount();
        rolj=new double[NumberOfAtoms];
        eolj=new double[NumberOfAtoms];
        rhs=new double[NumberOfAtoms];
    }
    /**
     * This method is used to calculate the LJ and ion-dipole potential
     * by the molecule on the gas atom: the gas atom is kept at the origin of a 3D cartesian plane.
     * @param x position of the gas atom ralative to the x axis
     * @param y position of the gas atom ralative to the y axis
     * @param  z position of the gas atom ralative to the z axis
     */
    public void potential(double x, double y, double z, char bufferGas)
    {

        double rx=0, ry=0, rz=0, e00=0, de00x=0, de00y=0, de00z=0;
        double sum1=0, sum2=0, sum3=0, sum4=0, sum5=0, sum6=0;
        double xx, xx2, yy, yy2, zz, zz2, rxyz2, rxyz;
        double rxyz3, rxyz8, rxyz14;

        dmax=2*ROMAX;
        int i=0;

        double de00;
        double rxyz3i, rxyz5i;


        IAtom a=null;

        if(bufferGas=='H') {
            for (Iterator i$ = Molecule.atoms().iterator(); i$.hasNext(); ) {


                a = (IAtom) i$.next();

                xx = x - a.getPoint3d().x;
                xx2 = xx * xx;
                yy = y - a.getPoint3d().y;
                yy2 = yy * yy;
                zz = z - a.getPoint3d().z;
                zz2 = zz * zz;
                rxyz = Math.sqrt(xx2 + yy2 + zz2);

                if (rxyz < dmax)
                    dmax = rxyz;




                // LJ potential
                e00 += ((eolj[i] * 4) * ((Math.pow(rolj[i], 12) / Math.pow(rxyz, 12)) - (Math.pow(rolj[i], 6) / Math.pow(rxyz, 6))));
                // LJ derivative
                de00 = (eolj[i] * 4) * (((6 * Math.pow(rolj[i], 6)) / Math.pow(rxyz, 8)) - ((12 * Math.pow(rolj[i], 12)) / Math.pow(rxyz, 14)));

                de00x += (de00 * xx);
                de00y += (de00 * yy);
                de00z += (de00 * zz);

                // ion-induced dipole potential
                a.setCharge(1.0 / (Molecule.getAtomCount()));

                if (a.getCharge() != 0) {
                    rxyz3i = a.getCharge() / Math.pow(rxyz, 3);
                    rxyz5i = -3 * a.getCharge() / Math.pow(rxyz, 5);
                    rx += (xx * rxyz3i);
                    ry += (yy * rxyz3i);
                    rz += (zz * rxyz3i);

                    // ion-induced dipole derivative
                    sum1 += (rxyz3i + (xx2 * rxyz5i));
                    sum2 += (xx * yy * rxyz5i);
                    sum3 += (xx * zz * rxyz5i);
                    sum4 += (rxyz3i + (yy2 * rxyz5i));
                    sum5 += (yy * zz * rxyz5i);
                    sum6 += (rxyz3i + (zz2 * rxyz5i));
                }


                i++;
            }

            pot = e00 - (DIPOL * ((rx * rx) + (ry * ry) + (rz * rz)));
            dpotx = de00x - (DIPOL * ((2 * rx * sum1) + (2 * ry * sum2) + (2 * rz * sum3)));
            dpoty = de00y - (DIPOL * ((2 * rx * sum2) + (2 * ry * sum4) + (2 * rz * sum5)));
            dpotz = de00z - (DIPOL * ((2 * rx * sum3) + (2 * ry * sum5) + (2 * rz * sum6)));

            return;

        }
        else
        {

            double[][] pottry=new double[3][3];
            double[] pot_mol=new double[4];
            double[][] dpotxtry=new double[3][3];
            double[]dpotx_mol=new double[3];
            double[][] dpotytry=new double[3][3];
            double[]dpoty_mol=new double[3];
            double[][] dpotztry=new double[3][3];
            double[] dpotz_mol=new double[3];

            double Ptfn=0;
            double xkT=IBST_MAX*XK;

            double dipolzz=(DIPOL_D/(8*PI*XEO))*Math.pow(XE,2);
            double dipolxx=(DIPOL_D/(8*PI*XEO))*Math.pow(XE,2);
            double POT_MIN=1.0E8;


            double qpol,dqpolx,dqpoly,dqpolz;


            double dpolx=0.0,dpoly=0.0,dpolz=0.0;
            double xc,yc,zc;
            double xx_center,xx_center2,yy_center,yy_center2,zz_center,zz_center2;
            double rxyz_center,rxyz_center2,rxyz_center3,rxyz_center5,const_k,temp_pot;

            for(int isamp=0;isamp<3;isamp++)
            {
                for (int ibatom = 0; ibatom < 3; ibatom++)
                {
                    rx = 0;
                    ry = 0;
                    rz = 0;
                    e00 = 0;
                    de00x = 0;
                    de00y = 0;
                    de00z = 0;
                    sum1 = 0;
                    sum2 = 0;
                    sum3 = 0;
                    sum4 = 0;
                    sum5 = 0;
                    sum6 = 0;
                    qpol = 0;
                    dqpolx = 0;
                    dqpoly = 0;
                    dqpolz = 0;
                    xc=0;
                    yc=0;
                    zc=0;
                    for (Iterator i$ = Molecule.atoms().iterator(); i$.hasNext(); ) {
                        a = (IAtom) i$.next();

                        if (isamp==0) {
                            xc = (BOND / 2)*(2 * (ibatom+1) - 3);
                        }
                        if(isamp==1)
                        {
                            yc = (BOND / 2)*(2 * (ibatom+1) - 3);
                        }
                        if(isamp==2) {
                            zc = (BOND / 2)*(2 * (ibatom+1) - 3);
                        }

                        dpolx = dipolzz;
                        dpoly = dipolxx;
                        dpolz = dipolxx;


                        xx_center=x- a.getPoint3d().x;
                        xx=xx_center+xc;
                        xx_center2= Math.pow(xx_center,2);
                        xx2=Math.pow(xx,2);

                        yy_center=y-a.getPoint3d().y;
                        yy=yy_center+yc;
                        yy_center2= Math.pow(yy_center,2);
                        yy2=Math.pow(yy,2);

                        zz_center=z-a.getPoint3d().z;
                        zz=zz_center+zc;
                        zz_center2=Math.pow(zz_center,2);
                        zz2=Math.pow(zz,2);

                        rxyz_center2=xx_center2+yy_center2+zz_center2;
                        rxyz2=xx2+yy2+zz2;

                        rxyz_center=Math.sqrt(rxyz_center2);
                        rxyz=Math.sqrt(rxyz2);
                        if(rxyz < dmax) {
                            dmax = rxyz;
                        }
                        rxyz3=Math.pow(rxyz,3);
                        rxyz8=Math.pow(rxyz,8);
                        rxyz14=Math.pow(rxyz,14);

                        rxyz_center3=Math.pow(rxyz_center,3);
                        rxyz_center5=Math.pow(rxyz_center,5);
                        // LJ potential
                        e00+=((eolj[i] * 4)*((Math.pow(rolj[i],12)/Math.pow(rxyz,12))- (Math.pow(rolj[i],6)/Math.pow(rxyz,6))));
                        // LJ derivative
                        de00=(eolj[i] * 4)*(((Math.pow(rolj[i],6)*6)/rxyz8)- ((Math.pow(rolj[i],12)*12)/rxyz14));

                        de00x+=(de00*xx);
                        de00y+=(de00*yy);
                        de00z+=(de00*zz);


                        if (a.getCharge() != null) {
                            rxyz3i=a.getCharge()/rxyz_center3;
                            rxyz5i=-3*a.getCharge()/rxyz_center5;
                            rx+=(xx_center*rxyz3i);
                            ry+=(yy_center*rxyz3i);
                            rz+=(zz_center*rxyz3i);

                            // ion-induced dipole derivative
                            sum1+=(rxyz3i+(xx_center2*rxyz5i));
                            sum2+=(xx_center*yy_center*rxyz5i);
                            sum3+=(xx_center*zz_center*rxyz5i);
                            sum4+=(rxyz3i+(yy_center2*rxyz5i));
                            sum5+=(yy_center*zz_center*rxyz5i);
                            sum6+=(rxyz3i+(zz_center2*rxyz5i));

                            //ion-partial charge coulomb potential(quadrupole)
                            const_k=a.getCharge()*(Math.pow(Constants.XE,2))/(4*Constants.PI*Constants.XEO);
                            qpol+=(PC_CENTER*const_k/rxyz_center);
                            qpol+=(PC*const_k/rxyz);

                            // ion-partial charge coulomb derivative(quadrupole)
                            dqpolx-=((PC_CENTER*const_k/rxyz_center3)*(xx_center));
                            dqpoly-=((PC_CENTER*const_k/rxyz_center3)*(yy_center));
                            dqpolz-=((PC_CENTER*const_k/rxyz_center3)*(zz_center));
                            dqpolx-=((PC*const_k/rxyz3)*(xx));
                            dqpoly-=((PC*const_k/rxyz3)*(yy));
                            dqpolz-=((PC*const_k/rxyz3)*(zz));


                        }


                        i++;

                    }
                    i=0;

                    pottry[isamp][ibatom]=e00-0.5*(((dpolx*rx*rx)+(dpoly*ry*ry) +(dpolz*rz*rz)))+qpol;
                    dpotxtry[isamp][ibatom]=de00x-0.5*((dpolx*2*rx*sum1)+(dpoly*2*ry*sum2)+(dpolz*2*rz*sum3))+dqpolx;
                    dpotytry[isamp][ibatom]=de00y-0.5*((dpolx*2*rx*sum2)+(dpoly*2*ry*sum4)+(dpolz*2*rz*sum5))+dqpoly;
                    dpotztry[isamp][ibatom]=de00z-0.5*((dpolx*2*rx*sum3)+(dpoly*2*ry*sum5)+(dpolz*2*rz*sum6))+dqpolz;


                }
                pot_mol[isamp]=pottry[isamp][0]+pottry[isamp][1];
                if(POT_MIN > pot_mol[isamp]) {
                    POT_MIN = pot_mol[isamp];
                }
                dpotx_mol[isamp]=dpotxtry[isamp][0]+dpotxtry[isamp][1];
                dpoty_mol[isamp]=dpotytry[isamp][0]+dpotytry[isamp][1];
                dpotz_mol[isamp]=dpotztry[isamp][0]+dpotztry[isamp][1];

            }

            for(int isamp=0;isamp<3;isamp++) {
                Ptfn = Ptfn + Math.exp(-(pot_mol[isamp] - POT_MIN)/ xkT);
            }


            pot=0;
            dpotx=0;
            dpoty=0;
            dpotz=0;

            double weight;

            for(int isamp=0;isamp<3;isamp++) {

                temp_pot = pot_mol[isamp] - POT_MIN;
                weight = Math.exp(-temp_pot / xkT) / Ptfn;
                pot += weight * pot_mol[isamp];
                dpotx += weight * dpotx_mol[isamp];
                dpoty += weight * dpoty_mol[isamp];
                dpotz += weight * dpotz_mol[isamp];
            }



        }
        return;

    }
    /**
     * This function rotates the cluster to a specified orientation
     * @param phi This is the angle between the x axis and the N axis
     * @param theta  This is the second parameter to addNum method
     * @param gamma  This is the angle between the N axis and the X axis
     */
    public void rotate(double phi, double theta, double gamma)
    {
        IAtom a=null;
        double rxy,rzy,nphi,otheta,ophi,ogamma;
        IAtomContainer M=null;

        try {
            M = (IAtomContainer)originalMolecule.clone();
        }
        catch (CloneNotSupportedException e)
        {
            //logger.debug(e);
        }

        for(Iterator i$ = M.atoms().iterator(); i$.hasNext();) {
            a = (IAtom)i$.next();
            rxy=Math.sqrt(Math.pow(a.getPoint3d().x,2)+Math.pow(a.getPoint3d().y,2));
            if(rxy==0)
            {
                a.getPoint3d().x=Math.cos(NTHETA)*rxy;
                a.getPoint3d().y=Math.sin(NTHETA)*rxy;

            }
            else{
                otheta=Math.acos(a.getPoint3d().x/rxy);
                if(a.getPoint3d().y < 0)
                    otheta=(2*Constants.PI)-otheta;
                NTHETA=otheta+theta;

                a.getPoint3d().x=Math.cos(NTHETA)*rxy;
                a.getPoint3d().y=Math.sin(NTHETA)*rxy;
            }
        }
        for(Iterator j$ = M.atoms().iterator(); j$.hasNext();) {
            a = (IAtom)j$.next();
            rzy=Math.sqrt(Math.pow(a.getPoint3d().z,2)+Math.pow(a.getPoint3d().y,2));
            if(rzy==0)
            {
                a.getPoint3d().z=Math.cos(NPHI)*rzy;
                a.getPoint3d().y=Math.sin(NPHI)*rzy;
            }
            else{
                ophi=Math.acos(a.getPoint3d().z/rzy);
                if(a.getPoint3d().y < 0)
                    ophi=(2*Constants.PI)-ophi;
                nphi=ophi+phi;

                a.getPoint3d().z=Math.cos(nphi)*rzy;
                a.getPoint3d().y=Math.sin(nphi)*rzy;
            }
        }

        for(Iterator k$ = M.atoms().iterator(); k$.hasNext();) {
            a = (IAtom)k$.next();
            rxy=Math.sqrt(Math.pow(a.getPoint3d().x,2)+Math.pow(a.getPoint3d().y,2));
            if(rxy==0)
            {
                a.getPoint3d().x=Math.cos(NGAMMA)*rxy;
                a.getPoint3d().y=Math.sin(NGAMMA)*rxy;
            }
            else{
                ogamma=Math.acos(a.getPoint3d().x/rxy);
                if(a.getPoint3d().y < 0)
                    ogamma=(2*Constants.PI)-ogamma;
                NGAMMA=ogamma+gamma;

                a.getPoint3d().x=Math.cos(NGAMMA)*rxy;
                a.getPoint3d().y=Math.sin(NGAMMA)*rxy;
            }
        }
        try {
            Molecule = (IAtomContainer)M.clone();
        }
        catch (CloneNotSupportedException e)
        {
            //logger.debug(e);
        }

    }
    /**
     * This method is used rorate the cluster to a random orientation
     * @param T This is a random number for theta
     * @param P This is a random number for phi
     * @param G This is a random number for gamma
     */
    public void rantate(double T,double P,double G)
    {

        double theta=T*2*Constants.PI;
        double phi=Math.asin((P*2)-1)+(Constants.PI/2);
        double gamma=G*2*Constants.PI;

        this.rotate(phi,theta,gamma);
    }
    /**
     * This method is used to assign basic parameters to the cluster
     */
    public void configureCluster(char bufferGas)
    {

        int i=0;double weight=0;
        IsotopeFactory factory;


        try {
            factory = Isotopes.getInstance();
            factory.configureAtoms(Molecule);
        }
        catch (IOException e)
        {
            //logger.debug(e);
        }



        for (IAtom atom : Molecule.atoms()) {
            weight+=atom.getExactMass();

            if(bufferGas=='N')
            {
                eolj[i] = this.setLjParametersN2("EO", atom.getSymbol());
                rolj[i] = this.setLjParametersN2("RO", atom.getSymbol());
                rhs[i] = this.setLjParametersN2("RH", atom.getSymbol());
            }
            else {
                eolj[i] = this.setLjParameters("EO", atom.getSymbol());
                rolj[i] = this.setLjParameters("RO", atom.getSymbol());
                rhs[i] = this.setLjParameters("RH", atom.getSymbol());
            }
            i++;
        }


        Constants.M2=weight;
        Point3d centerOfMass = GeometryUtil.get3DCentreOfMass(Molecule);

        for (IAtom atom : Molecule.atoms()) {


            atom.getPoint3d().x=(atom.getPoint3d().x-centerOfMass.x)*1.E-10;
            atom.getPoint3d().y=(atom.getPoint3d().y-centerOfMass.y)*1.E-10;
            atom.getPoint3d().z= (atom.getPoint3d().z-centerOfMass.z)*1.E-10;



            // above is a test code for center mass

        }

        try {
            originalMolecule = (IAtomContainer)Molecule.clone();
        }
        catch (CloneNotSupportedException e)
        {

        }

        double[] roljMax=rolj.clone();
        Arrays.sort(roljMax);
        ROMAX =roljMax[roljMax.length-1];

        if(bufferGas=='N')
        {
            ROMAX=ROMAX+(1.1055E-10/2.0);
        }
        dmax=2*ROMAX;
    }
    /**
     * This method is used to find LJ parameters of each Atom in the cluster for hilium
     * @param Parameter This is the required LJ parameter
     * @param AtomType This is the Type of the atom
     * @return double LJ Parameter value
     */
    public double setLjParameters(String Parameter, String AtomType)
    {
        double ParameterVaule=0;

        if(Parameter=="EO")
        {
            switch(AtomType)
            {
                case "O" :
                    ParameterVaule=1.34E-3*Constants.XE;
                    break;
                case "N" :
                    ParameterVaule=1.34E-3*Constants.XE;
                    break;
                case "C" :
                    ParameterVaule=1.34E-3*Constants.XE;
                    break;
                case "H" :
                    ParameterVaule=0.65E-03*Constants.XE;
                    break;
                case "S" :
                    ParameterVaule=1.35E-3*Constants.XE;
                    break;

            }
        }
        else if(Parameter=="RO")
        {

            switch(AtomType)
            {
                case "O" :
                    ParameterVaule=3.043E0*1.0E-10;
                    break;
                case "N" :
                    ParameterVaule=3.043E0*1.0E-10;
                    break;
                case "C" :
                    ParameterVaule=3.043E0*1.0E-10;
                    break;
                case "H" :
                    ParameterVaule=2.38E0*1.0E-10;
                    break;
                case "S" :
                    ParameterVaule=3.5E0*1.0E-10;
                    break;
            }
        }
        else if(Parameter=="RH")
        {
            switch(AtomType)
            {
                case "O" :
                    ParameterVaule=2.7*1.0E-10;
                    break;
                case "N" :
                    ParameterVaule=2.7*1.0E-10;
                    break;
                case "C" :
                    ParameterVaule=2.7*1.0E-10;
                    break;
                case "H" :
                    ParameterVaule=2.2*1.0E-10;
                    break;
                case "S" :
                    ParameterVaule=3.5*1.0E-10;
                    break;

            }
        }
        return ParameterVaule;
    }
    /**
     * This method is used to find LJ parameters of each Atom in the cluster for hilium
     * @param Parameter This is the required LJ parameter
     * @param AtomType This is the Type of the atom
     * @return double LJ Parameter value
     */
    public double setLjParametersN2(String Parameter, String AtomType)
    {
        double ParameterVaule=0;

        if(Parameter=="EO")
        {
            switch(AtomType)
            {
                case "O" :
                    ParameterVaule=Math.sqrt(EOGAS*0.0558)*CONVE*Constants.XE;
                    break;
                case "N" :
                    ParameterVaule=Math.sqrt(EOGAS*0.0828)*CONVE*Constants.XE;
                    break;
                case "C" :
                    ParameterVaule=Math.sqrt(EOGAS*0.0977)*CONVE*Constants.XE;
                    break;
                case "H" :
                    ParameterVaule=Math.sqrt(EOGAS*0.0189)*CONVE*Constants.XE;
                    break;
                case "Na":
                    ParameterVaule=Math.sqrt(EOGAS*0.03000)*CONVE*Constants.XE;
                    break;
                case "Si":
                    ParameterVaule=Math.sqrt(EOGAS*0.4020)*CONVE*Constants.XE;
                    break;
                case "S":
                    ParameterVaule=Math.sqrt(EOGAS*0.2740)*CONVE*Constants.XE;
                    break;
                case "Fe":
                    ParameterVaule=Math.sqrt(EOGAS*0.0130)*CONVE*Constants.XE;
                    break;
                case "P":
                    ParameterVaule=Math.sqrt(EOGAS*0.305)*CONVE*Constants.XE;
                    break;
                case "F":
                    ParameterVaule=Math.sqrt(EOGAS*0.0465)*CONVE*Constants.XE;
                    break;



            }
        }
        else if(Parameter=="RO")
        {

            switch(AtomType)
            {
                case "O" :
                    ParameterVaule=Math.sqrt(ROGAS*3.2550)*CONVR*1.0E-10;
                    break;
                case "N" :
                    ParameterVaule=Math.sqrt(ROGAS*4.3920)*CONVR*1.0E-10;
                    break;
                case "C" :
                    ParameterVaule=Math.sqrt(ROGAS*3.5814)*CONVR*1.0E-10;
                    break;
                case "H" :
                    ParameterVaule=Math.sqrt(ROGAS*1.2409)*CONVR*1.0E-10;
                    break;
                case "Na":
                    ParameterVaule=Math.sqrt(ROGAS*2.9830)*CONVR*1.0E-10;
                    break;
                case "Si":
                    ParameterVaule=Math.sqrt(ROGAS*4.2950)*CONVR*1.0E-10;
                    break;
                case "S":
                    ParameterVaule=Math.sqrt(ROGAS*4.0350)*CONVR*1.0E-10;
                    break;
                case "Fe":
                    ParameterVaule=Math.sqrt(ROGAS*2.9120)*CONVR*1.0E-10;
                    break;
                case "P":
                    ParameterVaule=Math.sqrt(ROGAS*4.1470)*CONVR*1.0E-10;
                    break;
                case "F":
                    ParameterVaule=Math.sqrt(ROGAS*3.1285)*CONVR*1.0E-10;
                    break;
            }
        }
        else if(Parameter=="RH")
        {
            switch(AtomType)
            {
                case "O" :
                    ParameterVaule=2.7*1.0E-10;
                    break;
                case "N" :
                    ParameterVaule=2.7*1.0E-10;
                    break;
                case "C" :
                    ParameterVaule=2.7*1.0E-10;
                    break;
                case "H" :
                    ParameterVaule=2.2*1.0E-10;
                    break;
                case "Na":
                    ParameterVaule=2.853*1.0E-10;
                    break;
                case "Si":
                    ParameterVaule=2.95*1.0E-10;
                    break;
                case "S":
                    ParameterVaule=3.5*1.0E-10;
                    break;
                case "Fe":
                    ParameterVaule=3.5*1.0E-10;
                    break;
                case "P":
                    ParameterVaule=4.2*1.0E-10;
                    break;
                case "F":
                    ParameterVaule=2.7*1.0E-10;
                    break;

            }
        }
        return ParameterVaule;
    }
    public Object clone()throws CloneNotSupportedException{
        return super.clone();
    }
}
