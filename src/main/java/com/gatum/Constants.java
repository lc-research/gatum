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

/**
 * <h1>Constants!</h1>
 * The class holds all the constants used for CCS calculation
 *
 * @author M. Dias
 * @version 1.0
 */
public final class Constants {
    private Constants(){}

    public static char bufferGas='H';
    public static final int temperature=298;

    public static final double PI=3.1415926535897931;
    public static final double XEO=8.854187817E-12;
    public static final double XE=1.60217733E-19;
    public static final double XK=1.380658E-23;
    public static final double XN=6.0221367E23;
    public static double NTHETA=3.6297469670513871E-308;
    public static double NPHI=6.9533558064128858E-310;
    public static double NGAMMA=6.9533558064128858E-310;
    public static double DIPOL=(0.204956E-30/(2*4*PI*XEO))*XE*XE;
    public static double N2_MASS=28.0134;
    public static double HE_MASS=4.0026;
    public static double M2; //Mass of the molecule.
    public static double MU; //Mass constants.
    public static double RMAX=0;
    public static double OMEGA;

    public static final int INP=40; //Number of points in velocity integration.
    public static final int ITN=10; //Number of complete cycles for average mobility calculation.
    public static final int IMP=25; // Number of points in Monte Carlo integrations of impact parameter and orientation.
    public static final int IRN=1000;
    public static final int IBST_MAX=500;
    public static final int IMP_FAC=5;
    public static final int MAX_STEPS=30000;
    public static final double CMIN=0.0005;
    public static final double DBST=0.1;

    // Lennard-Jones scaling parameters
    public static final double EO=1.34E-03*XE;
    public static final double RO=3.043*1.0E-10;
    public static final double DTSF1=0.5;
    public static final double  DTSF2=0.1;
    public static final double SW1=0.00005;
    public static final double EOGAS=0.06900;
    public static final double ROGAS =3.6600;
    public static final double CONVE=4.2*0.01036427;
    public static final double CONVR=0.890898718;
    public static final double SW2=0.005;
    public static final double BOND=1.0976E-10;
    public static final double PC=-0.4825;
    public static final double PC_CENTER=-(PC);
    public static final double DIPOL_D=1.710E-30;
    public static final double POT_MIN=1.0E8;
    public static double ROMAX=0.0;

}

