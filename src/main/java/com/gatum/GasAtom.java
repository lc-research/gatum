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

package com.gatum;

import org.openscience.cdk.interfaces.IAtom;
import java.util.Iterator;

/**
 * <h1>GasAtom!</h1>
 * This class defines the properties and functions of the buffer gas atom.
 * In addition to this, the class includes methods for a MD simulation.
 *
 * @author  M. Dias
 * @version 1.0
 */

public class GasAtom implements Cloneable{
    double[] w = new double[6];
    double[] dw = new double[6];
    double dt;
    int l;
    double[][] tmp_array=new double[6][6];
    double scatAngle=0.0;
    double M1=4.0026E0;

    public GasAtom(double[] w)
    {
        this.w=w.clone();
    }
    public GasAtom(char type)
    {
        if(type=='N')
            M1=28.0134E0;

        Constants.bufferGas=type;
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

        Mol.potential(w[0],w[2],w[4]);

        dw[1]=- Mol.dpotx;
        dw[3]=- Mol.dpoty;
        dw[5]=- Mol.dpotz;


    }
    /**
     * This method is used to find the trajectory of the buffer gas atom
     * @param Mol This is the molecule used for the MD simulation
     */
    public void diffeq(Cluster Mol)
    {

        if(l==0) {
            this.rungeKutta(Mol);
        }

        // Adams-moulton predictor-corrector
        double[] savw=w.clone();
        double[] savdw=dw.clone();


        for(int j=0;j<6;j++)
        {

            w[j] += ((savdw[j]+ (-0.111059153612E0 * tmp_array[0][j])+ (0.672667757774E0 * tmp_array[1][j])+ (-1.70633621697E0 * tmp_array[2][j])+ ( 2.33387888707E0* tmp_array[3][j])+ (-1.8524668225E0* tmp_array[4][j])) * (dt * 2.97013888888E0));
        }


        this.hamiltonian(Mol);


        for(int j=0;j<6;j++)
        {

            tmp_array[0][j]=tmp_array[1][j];
            tmp_array[1][j]=tmp_array[2][j];
            tmp_array[2][j]=tmp_array[3][j];
            tmp_array[3][j]=tmp_array[4][j];

            tmp_array[4][j]=savdw[j];
            w[j]=savw[j]+(dt * 0.990972222222E0)*(tmp_array[4][j]+(tmp_array[3][j]*-0.55921513665E0+(tmp_array[2][j]*0.337771548703E0+(tmp_array[1][j]*-0.121233356692E0+ (tmp_array[0][j]*0.0189208128941E0+0.332866152768E0*dw[j])))));
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
                r = 0.50E0 * (dw[i] - 2.0E0 * q[i]);
                w[i] += r;
                q[i] += 3.0 * r + (-0.5E0) * dw[i];
            }

            hamiltonian(Mol);
            for (int i = 0; i < 6; i++) {
                dw[i] = dt * dw[i];
                r = 0.292893218814E0 * (dw[i] - 1.0E0 * q[i]);
                w[i] += r;
                q[i] += 3.0 * r + (-0.292893218814E0) * dw[i];
            }

            hamiltonian(Mol);
            for (int i = 0; i < 6; i++) {
                dw[i] = dt * dw[i];
                r = 1.70710678118E0 * (dw[i] - 1.0E0 * q[i]);
                w[i] += r;
                q[i] += 3.0 * r + (-1.70710678118E0) * dw[i];
            }

            hamiltonian(Mol);
            for (int i = 0; i < 6; i++) {
                dw[i] = dt * dw[i];
                r = 0.1666666666667E0 * (dw[i] - 2.0E0 * q[i]);
                w[i] += r;
            }

            hamiltonian(Mol);

            tmp_array[j] = dw.clone();
        }

        l = -1;
        dt = 2.0 * dt;


    }
    /**
     * This method is used to calculate the scattering angle of the buffer gas atom.
     * @param Mol This is the molecule involved in the binary collision
     * @param v This is the initial velocity of the buffer gas atom.
     * @param  b This is the impact parameter
     */
    public void scatteringAngle(Cluster Mol,double v,double b) {


        double vy = -v;
        double vx = 0.E0;
        double vz = 0.E0;

        // determine time step

        double top = (v / 95.2381E0) - 0.5E0;
        if (v > 1000.E0)
            top = 10.E0;
        if (v > 2000.E0)
            top = 10.0E0 - ((v - 2000.E0) * 7.5E-3);
        if (v > 3000.E0)
            top = 2.5E0;

        double dt1 = top * Constants.DTSF1 * 1.0E-11 / v;
        double dt2 = dt1 * Constants.DTSF2;

        // determine trajectory start position

        double e0 = 0.5E0 * Constants.MU * Math.pow(v, 2);
        double x = b;
        double z = 0.E0;

        double ymin = 0.E0;
        double ymax = 0.E0;

        IAtom a;

        for (Iterator i$ = Mol.Molecule.atoms().iterator(); i$.hasNext(); ) {
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


        Mol.potential(x,y,z);



        if (Math.abs((Mol.pot) / e0) < Constants.SW1) {

            do {
                id2 = id2 - 1;
                y = (double) id2 * 1.0E-10;
                Mol.potential(x,y,z);
                if (id2 < iymin) {
                    return;
                }

                y = (double) id2 * 1.0E-10;
                Mol.potential(x,y,z);

            } while ((Math.abs(Mol.pot / e0)) < Constants.SW1);
        } else {
            do {
                id2 = id2 + 10;
                y = (double) id2 * 1.0E-10;
                Mol.potential(x,y,z);

            } while (Math.abs(Mol.pot / e0) > Constants.SW1);

            do {
                id2 = id2 - 1;
                y = (double) id2 * 1.0E-10;
                Mol.potential(x,y,z);
            } while (Math.abs(Mol.pot / e0) < Constants.SW1);

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
        bufferGas.l = 0;
        double redPot=0.0;

        do {
            do {
                    bufferGas.diffeq(Mol);

                nw++;
                if ((noOfSteps + nw) > 30000)// check if trajectory lost: too many steps is an indication of a lost trajectory
                {
                    return;
                }
            } while (Mol.dmax<Constants.ROMAX);

            redPot = Math.abs(Mol.pot) / e0;
            if ((redPot > Constants.SW2) && bufferGas.dt == dt1) {
                bufferGas.dt = dt2;
                bufferGas.l=0;
            }


            if ((redPot < Constants.SW2) && bufferGas.dt == dt2) {
                bufferGas.dt = dt1;
                bufferGas.l=0;
            }
        } while (((redPot > Constants.SW1)) || (nw < 50));


        scatAngle = Math.acos((bufferGas.dw[2] * (-v)) / (v * Math.sqrt(Math.pow(bufferGas.dw[0], 2) + Math.pow(bufferGas.dw[2], 2) + Math.pow(bufferGas.dw[4], 2))));


        if (!(bufferGas.dw[0] > 0.E0)) {
            scatAngle=scatAngle*-1;
        }


    }
    public Object clone()throws CloneNotSupportedException{
        return super.clone();
    }
}
