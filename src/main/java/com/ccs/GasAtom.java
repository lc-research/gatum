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

import org.openscience.cdk.interfaces.IAtom;
import java.util.Iterator;
import static com.ccs.Constants.*;
/**
 * <h1>GasAtom!</h1>
 * This class defines the properties and functions of the buffer gas atom.
 * In addition to this, the class includes methods for a MD simulation.
 *
 * @author  M. Dias
 * @version 1.0
 */

public class GasAtom implements Cloneable{
    private char type='H';
    private double[] w = new double[6];
    private double[] dw = new double[6];
    private double dt;
    private boolean initiator=true;
    private double[][] tmp_array=new double[6][6];
    private double scatAngle;
    private double Mass;
    private final double SW1=0.00005;
    private final double SW2=0.005;
    private double N2_MASS=28.0134;
    private double HE_MASS=4.0026;
    private final int MAX_STEPS=30000;
    private final double DTSF1=0.5;
    private final double  DTSF2=0.1;

    public GasAtom(double[] w)
    {
        this.w=w.clone();
    }
    public GasAtom(char type)
    {
        Mass = (type=='N')? N2_MASS : HE_MASS;
        type=type;
    }

    public char getType()
    {
        return type;
    }

    public double getScatAngle(){
        return scatAngle;
    }

    public double getMass()
    {
        return this.Mass;
    }
    /**
     * This method is used to difine hamiltonian equations of motion.
     * @param Mol This is the molecule used in the MD simulation.
     */
    public void hamiltonian(Cluster Mol)
    {

        dw[0]=w[1]/Constants.MU;
        dw[2]=w[3]/Constants.MU;
        dw[4]=w[5]/Constants.MU;

        Mol.potential(w[0],w[2],w[4],type);

        dw[1]=- Mol.getDpotx();
        dw[3]=- Mol.getDpoty();
        dw[5]=- Mol.getDpotz();


    }
    /**
     * This method is used to find the trajectory of the buffer gas atom
     * @param Mol This is the molecule used for the MD simulation
     */
    public void diffeq(Cluster Mol)
    {

        if(initiator) {
            this.rungeKutta(Mol);
        }

        // Adams-moulton predictor-corrector
        double[] savw=w.clone();
        double[] savdw=dw.clone();


        for(int j=0;j<6;j++)
        {

            w[j] += ((savdw[j]+ (-0.111059153612 * tmp_array[0][j])+ (0.672667757774 * tmp_array[1][j])+ (-1.70633621697 * tmp_array[2][j])+ ( 2.33387888707* tmp_array[3][j])+ (-1.8524668225* tmp_array[4][j])) * (dt * 2.97013888888));
        }


        this.hamiltonian(Mol);


        for(int j=0;j<6;j++)
        {

            tmp_array[0][j]=tmp_array[1][j];
            tmp_array[1][j]=tmp_array[2][j];
            tmp_array[2][j]=tmp_array[3][j];
            tmp_array[3][j]=tmp_array[4][j];

            tmp_array[4][j]=savdw[j];
            w[j]=savw[j]+(dt * 0.990972222222)*(tmp_array[4][j]+(tmp_array[3][j]*-0.55921513665+(tmp_array[2][j]*0.337771548703+(tmp_array[1][j]*-0.121233356692+ (tmp_array[0][j]*0.0189208128941+0.332866152768*dw[j])))));
        }


        hamiltonian(Mol);


    }
    /**
     * This method is used to generate the first 5 positions of the buffer gas atom in its trajectory
     * @param Mol This is the molecule used for the MD simulation
     */
    public void rungeKutta(Cluster Mol)
    {
        double[] q= new double[6];

        dt = 0.5 * dt;

        double r;

        for(int j=0;j<5;j++) {

            hamiltonian(Mol);
            for (int i = 0; i < 6; i++) {
                dw[i] = dt * dw[i];
                r = 0.50 * (dw[i] - 2.0 * q[i]);
                w[i] += r;
                q[i] += 3.0 * r + (-0.5) * dw[i];
            }

            hamiltonian(Mol);
            for (int i = 0; i < 6; i++) {
                dw[i] = dt * dw[i];
                r = 0.292893218814 * (dw[i] - 1 * q[i]);
                w[i] += r;
                q[i] += 3.0 * r + (-0.292893218814) * dw[i];
            }

            hamiltonian(Mol);
            for (int i = 0; i < 6; i++) {
                dw[i] = dt * dw[i];
                r = 1.70710678118 * (dw[i] - 1 * q[i]);
                w[i] += r;
                q[i] += 3.0 * r + (-1.70710678118) * dw[i];
            }

            hamiltonian(Mol);
            for (int i = 0; i < 6; i++) {
                dw[i] = dt * dw[i];
                r = 0.1666666666667 * (dw[i] - 2 * q[i]);
                w[i] += r;
            }

            hamiltonian(Mol);

            tmp_array[j] = dw.clone();
        }

        initiator = false;
        dt = 2.0 * dt;


    }
    /**
     * This method is used to calculate the scattering angle of the buffer gas atom.
     * @param Mol This is the molecule involved in the binary collision
     */
    public void scatteringAngle(Cluster Mol,double velocity,double impactParameter) {


        double vy = -velocity;
        double vx = 0;
        double vz = 0;

        // determine time step

        double top = (velocity / 95.2381) - 0.5;
        if (velocity > 1000)
            top = 10;
        if (velocity > 2000)
            top = 10.0 - ((velocity - 2000) * 7.5E-3);
        if (velocity > 3000)
            top = 2.5;

        double dt1 = top * DTSF1 * 1.0E-11 / velocity;
        double dt2 = dt1 * DTSF2;

        // determine trajectory start position

        double e0 = 0.5 * Constants.MU * Math.pow(velocity, 2);
        double x = impactParameter;
        double z = 0;

        double ymin = 0;
        double ymax = 0;

        IAtom a;

        for (Iterator i$ = Mol.getMolecule().atoms().iterator(); i$.hasNext(); ) {
            a = (IAtom) i$.next();
            if (a.getPoint3d().y > ymax)
                ymax = a.getPoint3d().y;
            if (a.getPoint3d().y < ymin)
                ymin = a.getPoint3d().y;
        }


        ymax = ymax / 1.0E-10;
        ymin = ymin / 1.0E-10;
        int iymin = (int) ymin - 1;
        int iymax = (int) ymax + 1;
        int id2 = iymax;
        double y = (double) id2 * 1.0E-10;


        Mol.potential(x,y,z,type);



        if (Math.abs((Mol.getPot()) / e0) < SW1) {

            do {
                id2 = id2 - 1;
                y = (double) id2 * 1.0E-10;
                Mol.potential(x,y,z,type);
                if (id2 < iymin) {
                    return;
                }

                y = (double) id2 * 1.0E-10;
                Mol.potential(x,y,z,type);

            } while ((Math.abs(Mol.getPot() / e0)) < SW1);
        } else {
            do {
                id2 = id2 + 10;
                y = (double) id2 * 1.0E-10;
                Mol.potential(x,y,z,type);

            } while (Math.abs(Mol.getPot() / e0) > SW1);

            do {
                id2 = id2 - 1;
                y = (double) id2 * 1.0E-10;
                Mol.potential(x,y,z,type);
            } while (Math.abs(Mol.getPot() / e0) < SW1);

        }


        y = (double) id2 * 1.0E-10; // trajectory start position

        // initial coordinates and momenta
        double[] w = new double[7];
        w[0] = x;
        w[1] = vx * Constants.MU;
        w[2] = y;
        w[3] = vy * Constants.MU;
        w[4] = z;
        w[5] = vz * Constants.MU;

        // initialize the time derivatives of the coordinates and momenta
        GasAtom bufferGas=new GasAtom(w);
        bufferGas.dt=dt1;


        int noOfSteps = 0;
        int nw = 0;
        double redPot=0.0;

        do {
            do {
                bufferGas.diffeq(Mol);

                nw++;
                if ((noOfSteps + nw) > MAX_STEPS)// check if trajectory lost: too many steps is an indication of a lost trajectory
                {
                    return;
                }
            } while (Mol.getDmax()<Constants.ROMAX);

            redPot = Math.abs(Mol.getPot()) / e0;
            if ((redPot > SW2) && bufferGas.dt == dt1) {
                bufferGas.dt = dt2;
                bufferGas.initiator=false;
            }


            if ((redPot < SW2) && bufferGas.dt == dt2) {
                bufferGas.dt = dt1;
                bufferGas.initiator=false;
            }
        } while (((redPot > SW1)) || (nw < 50));


        scatAngle = Math.acos((bufferGas.dw[2] * (-velocity)) / (velocity * Math.sqrt(Math.pow(bufferGas.dw[0], 2) + Math.pow(bufferGas.dw[2], 2) + Math.pow(bufferGas.dw[4], 2))));


        if (!(bufferGas.dw[0] > 0)) {
            scatAngle=scatAngle*-1;
        }


    }
    public Object clone()throws CloneNotSupportedException{
        return super.clone();
    }
}
