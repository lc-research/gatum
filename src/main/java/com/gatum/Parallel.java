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
 * <h1>Parallel!</h1>
 * This class contains computationally expensive functions used to calculate CCS
 *
 *
 * @author  M. Dias
 * @version 1.0
 */
public class Parallel {

    public Parallel(int win)
    {
        window=win;
    }
    @Override
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    public int window;
    /**
     * This method is used to perform the MC integration
     * @param ion This is the molecule used for the CCS.
     * @param bufferGas This is the buffer gas atom used for collisions.
     * @param b2max This is an arry of impact parameters used to perform the MD simulation
     * @param pgst This is a parameter for the MD simulation
     * @param wgst This is a parameter for the MD simulation
     * @param random This is the random number sequence used in the MC integration
     * @param start The starting position of the parallelized work load
     * @param end The end position of the parallelized work load
     * @return int This returns sum of numA and numB.
     */
    public double mcIntegration(Cluster ion, GasAtom bufferGas, double [] b2max,double[] pgst,double[] wgst,double[] random, int start,int end) {

        int ig=0,r=start;
        if(r==0)
        {
            window=end - start;
        }
        else{
            r = (start / window) * window * Constants.IMP * 4;
        }
        ig=start%40;
        double v,b,mom11st=0.0;

        for(int j=start; j<end;j++)
        {
            v=Math.sqrt((Math.pow(pgst[ig],2)*Constants.EO)/(0.5E0* Constants.MU));

            for(int im=0;im<Constants.IMP;im++)
            {
                ion.rantate(random[r+1],random[r+2],random[r+3]);

                b=Constants.RO*Math.sqrt(random[r]*b2max[ig]);
                bufferGas.scatteringAngle(ion,v,b);

                mom11st+=((1.E0-Math.cos(bufferGas.scatAngle))*b2max[ig]*wgst[ig])/Constants.IMP;
                r+=4;

            }


            if(ig==39)
                ig=0;
            else
                ig++;
        }
        return mom11st;
    }
    /**
     * This method is used to calculate possible impact parameter values for the MD simulation
     * @param ion This is the molecule used for the MD simulation
     * @param bufferGas This is the buffer gas atom used for the MD simulation
     * @param pgst This is a parameter for the MD simulation
     * @param rmaxx This is the
     * @param start Starting position of the parallelized workload
     * @param end Ending postion of the parallelized workload
     */
    public double[] b2Max(Cluster ion, GasAtom bufferGas, double[] pgst,double rmaxx,int start,int end)
    {
        double[] b2max=new double[40];
        double dbst2 = 1.E0;
        double dbst22 = dbst2 / 10.E0;
        double cmin = 0.0005, gst2;
        double v, b;
        int ibst;
        double ang;
        double[] cosx = new double[500];


        for (int ig = start - 1; ig >=end; ig--) {
            gst2 = Math.pow(pgst[ig], 2);
            v = Math.sqrt((gst2 * Constants.EO) / (0.5E0 * Constants.MU));
            ibst = (int) (rmaxx / Constants.RO) - 6;

            if (ig < Constants.INP - 1)
                ibst = (int) (b2max[ig + 1] / dbst2) - 6;

            if (ibst < 0)
                ibst = 0;


            do {

                b = Constants.RO * Math.sqrt(dbst2 * (double) ibst);

                bufferGas.scatteringAngle(ion, v, b);
                ang = bufferGas.scatAngle;
                cosx[ibst] = 1.E0 - Math.cos(ang);

                if (ibst < 5) {
                    ibst = ibst + 1;
                }

                if (ibst >= 5) {
                    if (cosx[ibst] < cmin && cosx[ibst - 1] < cmin && cosx[ibst - 2] < cmin && cosx[ibst - 3] < cmin && cosx[ibst - 4] < cmin) {
                        b2max[ig] = (double) (ibst - 5) * dbst2;
                        do {
                            b2max[ig] = b2max[ig] + dbst22;
                            b = Constants.RO * Math.sqrt(b2max[ig]);
                            bufferGas.scatteringAngle(ion, v, b);
                            ang = bufferGas.scatAngle;
                        } while ((1.E0 - Math.cos(ang)) > cmin);
                        break;

                    } else {
                        ibst = ibst + 1;
                    }

                }

            } while (ibst < 500);
        }
        return b2max;
    }

}
