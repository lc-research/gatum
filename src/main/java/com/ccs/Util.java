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

import edu.cornell.lassp.houle.RngPack.Ranlux;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.IntStream;
import static com.ccs.Constants.*;
import static org.math.array.DoubleArray.*;
import static org.math.array.LinearAlgebra.*;
import static org.math.array.StatisticSample.*;
/**
 * <h1>Util!</h1>
 * This class provides utility functions to calculate CCS
 *
 *
 * @author  M. Dias
 * @version 1.0
 */

public class Util{
    private final int INP=40; //Number of points in velocity integration.
    private final int ITN=10; //Number of complete cycles for average mobility calculation.
    private final int IMP=25; // Number of points in Monte Carlo integrations of impact parameter and orientation.
    private final int IRN=1000;
    private final int IBST_MAX=500;
    private final int IMP_FAC=5;
    private final double CMIN=0.0005;
    private final double DBST=0.1;
    private double OMEGA;
    private final double RO=3.043*1.0E-10;


    public Util() {

    }

    /**
     * This method calls the required functions for the CSS calculation.
     * @param ion The calculations will be performed to calculate the CCS of this molecule.
     * @param bufferG This is the buffer gas atom of the CCS
     */
    public String collisionCrossSection(Cluster ion,char bufferG, int temperature)throws InterruptedException, ExecutionException {
        //determin maximum extent and orientation along x axis
        double RMAX=0;
        IAtom atom = null;
        double r;
        int ihold = 0, i = 0;

        final GasAtom bufferGas = new GasAtom(bufferG);
        ion.configureCluster(bufferGas.getType());


        for (Iterator i$ = ion.getOriginalMolecule().atoms().iterator(); i$.hasNext(); ) {

            atom = (IAtom) i$.next();

            r = Math.sqrt(Math.pow(atom.getPoint3d().x, 2) + Math.pow(atom.getPoint3d().y, 2) + Math.pow(atom.getPoint3d().z, 2));

            if (r > RMAX) {
                RMAX = r;
                ihold = i;
            }
            i++;

        }


        atom = ion.getOriginalMolecule().getAtom(ihold);
        double phi;

        phi = Math.acos(atom.getPoint3d().z / (Math.sqrt(Math.pow(atom.getPoint3d().z, 2) + Math.pow(atom.getPoint3d().y, 2))))+(Constants.PI / 2.0);
        if (atom.getPoint3d().y > 0)
            phi = (2 * Constants.PI) - phi;

        double theta = 0, gamma = 0;
        ion.rotate(phi, theta, gamma);


        atom = ion.getMolecule().getAtom(ihold);
        gamma = Math.acos(atom.getPoint3d().x / (Math.sqrt(Math.pow(atom.getPoint3d().x, 2) + Math.pow(atom.getPoint3d().y, 2))));

        if (atom.getPoint3d().y < 0)
            gamma = (2 * Constants.PI - gamma);

        gamma = (2 * Constants.PI) - gamma;

        ion.rotate(phi, theta, gamma);



        //      determine rmax, emax, and r00 along x, y, and z directions

        double x;

        for (int ir = 1; ir <= IRN; ir++) {
            x = RMAX + ROMAX - (ir * ((RMAX + ROMAX) / IRN));

            ion.potential(x, 0, 0,bufferGas.getType());

            if (ion.getPot() == 0) {
                break;
            }
            if (ion.getPot() < 0) {
                RMAX = x;
            }

        }


        // setup integration over gst
        Map rGst = this.integrationOverGst(temperature);
        final double[] pgst = (double[]) rGst.get("pgst");
        final double[] wgst = (double[]) rGst.get("wgst");


        final double[] b2max = new double[INP];


        MU=((bufferGas.getMass()*M2)/(bufferGas.getMass()+M2))/(XN*1.0E3);

        final Cluster[] ionParallel = new Cluster[INP*ITN];
        final GasAtom[] bufferGasParallel = new GasAtom[INP*ITN];
        for(int c=0;c<INP*ITN;c++){
            try{
                ionParallel[c]=(Cluster) ion.clone();
                bufferGasParallel[c]=(GasAtom)bufferGas.clone();
            }catch(CloneNotSupportedException e){}

        }

        final double EO_DIV_MU=EO/(MU/2);
        final int IBST = (int) (RMAX / RO) - 6;


        IntStream.iterate(INP-1,ig->ig-1).limit(INP).parallel().forEach(ig-> {

            double scatAngle,impactParamter;
            double[] cosx = new double[IBST_MAX];

            double velocity = Math.sqrt((Math.pow(pgst[ig], 2) * EO_DIV_MU));
            int ibst = IBST;

            if (ig < INP - 1)
                ibst = (int) (b2max[ig + 1] / DBST) - 6;

            if (ibst < 0)
                ibst = 0;


            do {

                impactParamter = RO * Math.sqrt(DBST * (double) ibst);

                bufferGas.scatteringAngle(ionParallel[ig], velocity, impactParamter);
                scatAngle = bufferGas.getScatAngle();
                cosx[ibst] = 1 - Math.cos(scatAngle);

                if (ibst < IMP_FAC) {
                    ibst++;
                }

                if (ibst >= IMP_FAC) {
                    if (cosx[ibst] < CMIN && cosx[ibst - 1] < CMIN && cosx[ibst - 2]
                            < CMIN && cosx[ibst - 3] < CMIN && cosx[ibst - 4] < CMIN) {
                        b2max[ig] = (double) (ibst - IMP_FAC) * DBST;
                        do {
                            b2max[ig] = b2max[ig] + DBST;
                            impactParamter = RO * Math.sqrt(b2max[ig]);
                            bufferGas.scatteringAngle(ionParallel[ig], velocity, impactParamter);
                            scatAngle = bufferGas.getScatAngle();
                        } while ((1 - Math.cos(scatAngle)) > CMIN);
                        break;

                    } else {
                        ibst++;
                    }

                }

            } while (ibst < IBST_MAX);
        } );



        double cs;

        Ranlux Sequence = new Ranlux(3,5013486);

        Sequence.raw();

        final double[] RandomN = new double[INP*ITN*IMP*4];
        for (int j = 0; j < INP*ITN*IMP*4; j++) {
            RandomN[j] = Sequence.raw();
        }

        OMEGA=0;
        IntStream.range(0,INP*ITN).parallel().forEach(j-> {
            double velocity,impactParameter;
            int ig=j%INP;
            int pointer=j;
            if(j>0)
                pointer=IMP*j*4;

            velocity=Math.sqrt(Math.pow(pgst[ig],2)*EO_DIV_MU);

            for(int im=0;im<IMP;im++)
            {
                ionParallel[j].rantate(RandomN[pointer+1],RandomN[pointer+2],RandomN[pointer+3]);

                impactParameter=RO*Math.sqrt(RandomN[pointer]*b2max[ig]);
                bufferGasParallel[j].scatteringAngle(ionParallel[j],
                        velocity,impactParameter);

                OMEGA +=((1-Math.cos(bufferGasParallel[j].getScatAngle()))*b2max[ig]*wgst[ig])/IMP;
                pointer+=4;

            }

        });


        cs = (OMEGA / (double) ITN) * Constants.PI * Math.pow(RO, 2)*Math.pow(10,20);

        return "The CCS value of the molecule " + ion.getMolecule().getProperty("cdk:Title") + " is " + String.format("%1.3E",cs)+'\n';



    }

    /**
     * This method is used to perform the integration over velocity
     */
    public Map integrationOverGst(int temperature)
    {
        double eo=1.34E-03*Constants.XE;
        double tst=Constants.XK*temperature/eo;
        double tst3=Math.pow(tst,3);

        double dgst=5.0E-7*6*Math.sqrt(tst);
        double gst=dgst;
        double sum=0, sum1=0, sum2=0;

        for(int i=0;i<=INP;i++) {

            sum1 = sum1 + Math.sqrt(i);
        }

        double hold1, hold2, gstt;
        double[] wgst = new double[INP];
        double[] pgst = new double[INP];

        for(int i=1;i<=INP;i++) {
            hold1 = Math.sqrt(i);
            hold2 = Math.sqrt(i-1);
            sum2 += hold2;
            wgst[i-1] = hold1 / sum1;

            gstt = tst3 * (sum2 + (hold1 / 2))/sum1;

            do {
                sum += (Math.exp(-1*Math.pow(gst,2)/ tst)* Math.pow(gst,5) * dgst);
                gst = gst + dgst;
                if (sum>gstt)
                    pgst[i-1] = gst - (dgst / 2);
            }while (sum<gstt);

        }

        Map  rGst= new HashMap();
        rGst.put("pgst",pgst);
        rGst.put("wgst",wgst);

        return rGst;
    }

    public String Response(String dataFile)throws InterruptedException, ExecutionException
    {

        DataFile inputFile=new DataFile("input_file/"+dataFile);
        IAtomContainer molecule[]= inputFile.retMolecule();
        String CCS="False";
        long start,end;

        for(int i=0;i<molecule.length;i++)
        {
            start = System.currentTimeMillis();
            Cluster ion = new Cluster(molecule[i]);
            CCS=this.collisionCrossSection(ion,'H',298);

            end = System.currentTimeMillis();
            float sec = (end - start)/ 1000F;

            System.out.println(CCS+" : "+sec+" seconds");

        }
        return CCS;
    }

}
