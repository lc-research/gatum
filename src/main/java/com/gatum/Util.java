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
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.*;
import java.util.concurrent.*;



/**
 * <h1>Util!</h1>
 * This class provides utility functions to calculate CCS
 *
 *
 * @author  M. Dias
 * @version 1.0
 */

public class Util{
    public Util() {

    }

    /**
     * This method calls the required functions for the CSS calculation.
     * @param ion The calculations will be performed to calculate the CCS of this molecule.
     * @param bufferG This is the buffer gas atom of the CCS
     */
    public String collisionCrossSection(Cluster ion,char bufferG)throws InterruptedException, ExecutionException {
        //determin maximum extent and orientation along x axis


        IAtom atom = null;
        double r, rmax=0.0;
        int ihold = 0, i = 0;

        final GasAtom bufferGas = new GasAtom(bufferG);
        ion.configureCluster();



        for (Iterator i$ = ion.originalMolecule.atoms().iterator(); i$.hasNext(); ) {

            atom = (IAtom) i$.next();

            r = Math.sqrt(Math.pow(atom.getPoint3d().x, 2) + Math.pow(atom.getPoint3d().y, 2) + Math.pow(atom.getPoint3d().z, 2));

            if (r > rmax) {
                rmax = r;
                ihold = i;
            }
            i++;

        }


        atom = ion.originalMolecule.getAtom(ihold);
        double phi;

        phi = Math.acos(atom.getPoint3d().z / (Math.sqrt(Math.pow(atom.getPoint3d().z, 2) + Math.pow(atom.getPoint3d().y, 2))))+(Constants.PI / 2.0);
        if (atom.getPoint3d().y > 0)
            phi = (2.0 * Constants.PI) - phi;

        double theta = 0.0, gamma = 0.0;
        ion.rotate(phi, theta, gamma);


        atom = ion.Molecule.getAtom(ihold);
        gamma = Math.acos(atom.getPoint3d().x / (Math.sqrt(Math.pow(atom.getPoint3d().x, 2) + Math.pow(atom.getPoint3d().y, 2))));

        if (atom.getPoint3d().y < 0)
            gamma = (2 * Constants.PI - gamma);

        gamma = (2 * Constants.PI) - gamma;

        ion.rotate(phi, theta, gamma);



        //      determine rmax, emax, and r00 along x, y, and z directions

        int irn = 1000;

        double x;

        for (int ir = 1; ir <= irn; ir++) {
            x = rmax + Constants.ROMAX - (ir * ((rmax + Constants.ROMAX) / irn));

            ion.potential(x, 0.E0, 0.E0);

            if (ion.pot == 0) {
                break;
            }
            if (ion.pot < 0.E0) {
                rmax = x;
            }

        }


        // setup integration over gst
        Map rGst = this.integrationOverGst();
        final double[] pgst = (double[]) rGst.get("pgst");
        final double[] wgst = (double[]) rGst.get("wgst");


        final double[] b2max = new double[40];

        int numOfCpus = Runtime.getRuntime().availableProcessors();
        int iterations = 40;
        int window=iterations/numOfCpus;
        final double rmaxxT=rmax;

        final int[] stepsB=new int[numOfCpus*2];
        int h=0;
        for(i=0;i<=iterations;)
        {
            stepsB[h]=i;
            i+=window;
            h++;
        }

        if(iterations%numOfCpus>0)
        {
            stepsB[numOfCpus]=iterations;
        }

        final Parallel[] integration = new Parallel[numOfCpus];
        final Cluster[] ion1 = new Cluster[numOfCpus];
        final GasAtom[] bufferGass = new GasAtom[numOfCpus];
        ExecutorService exec = Executors.newFixedThreadPool(numOfCpus);
        Collection<Callable<double[]>> tasksB2MaxMain = new ArrayList<Callable<double[]>>();

        Constants.MU=((bufferGas.M1*Constants.M2)/(bufferGas.M1+Constants.M2))/(Constants.XN*1.0E3);



        for(int c=0;c<numOfCpus;c++){
            final int m=c;
            try{
                integration[c]=new Parallel(window);
                ion1[c]=(Cluster) ion.clone();
                bufferGass[c]=(GasAtom)bufferGas.clone();

                tasksB2MaxMain.add((new Callable<double[]>() {
                    @Override
                    public double[] call() {
                        return integration[m].b2Max(ion1[m],bufferGass[m], pgst, rmaxxT, stepsB[m+1], stepsB[m]);
                    }
                }));
            }catch(CloneNotSupportedException e){}

        }

        List<Future<double[]>> futuresB2Max = exec.invokeAll(tasksB2MaxMain);

        for (int j = 0; j < numOfCpus; j++) {
            Future<double[]> B2MaxA = futuresB2Max.get(j);
            double[] b2Max_tmp1 = B2MaxA.get();
            for (i = 0; i < 40; i++) {
                if (b2Max_tmp1[i] > 0)
                    b2max[i] = b2Max_tmp1[i];
            }
        }


        double cs;
        double mom11st = 0.E0;


        Random.ranlux();


        final double[] RandomN = new double[40000];
        for (int j = 0; j < 40000; j++) {
            RandomN[j] = Random.ranlux();
        }

        iterations = Constants.INP * Constants.ITN;

        window=iterations/numOfCpus;
        final int[] steps=new int[numOfCpus*2];
        h=0;
        for(i=0;i<=iterations;)
        {
            steps[h]=i;
            i+=window;
            h++;
        }

        if(iterations%numOfCpus>0)
        {
            steps[numOfCpus]=iterations;
        }

        Collection<Callable<Double>> mcI = new ArrayList<Callable<Double>>();

            for(int c=0;c<numOfCpus;c++) {
                final int j=c;

                mcI.add((new Callable<Double>() {
                    @Override
                    public Double call() {
                        return integration[j].mcIntegration(ion1[j], bufferGass[j], b2max, pgst, wgst, RandomN, steps[j],steps[j+1]);
                    }
                }));
            }




            List<Future<Double>> mcIntegration = exec.invokeAll(mcI );
            exec.shutdown();


        for(i=0;i<numOfCpus;i++)
        {
            mom11st+=mcIntegration.get(i).get();
        }

            cs = (mom11st / (double) Constants.ITN) * Constants.PI * Math.pow(Constants.RO, 2)*Math.pow(10,20);
            //System.out.print("The CCS value of " + ion.Molecule.getProperty("cdk:Title") + " :" + Double.toString(cs));
        System.out.print(ion.Molecule.getProperty("cdk:Title") + " ," +  String.format("%1.3E",cs)+"\n");
        return "The CCS value of the molecule " + ion.Molecule.getProperty("cdk:Title") + " is " + String.format("%1.3E",cs)+'\n';



    }

    /**
     * This method is used to perform the integration over velocity
     */
    public Map integrationOverGst()
    {
        double eo=1.34E-03*Constants.XE;
        double tst=Constants.XK*Constants.temperature/eo;



        double tst3=Math.pow(tst,3);

        double dgst=5.0E-7*6.E0*Math.sqrt(tst);
        double gst=dgst;
        double sum=0.E0;
        double sum1=0.E0;
        double sum2=0.E0;

        for(int i=0;i<=Constants.INP;i++) {

            sum1 = sum1 + Math.sqrt(i);
        }

        double hold1=0.0, hold2=0.0, hold3=0.0, gstt=0.0;
        double[] wgst = new double[40];
        double[] pgst = new double[40];

        for(int i=1;i<=Constants.INP;i++) {
            hold1 = Math.sqrt(i);
            hold2 = Math.sqrt(i-1);
            sum2 = sum2 + hold2;
            wgst[i-1] = hold1 / sum1;

            gstt = tst3 * (sum2 + (hold1 / 2.E0))/sum1;

            do {
                sum = sum + (Math.exp(-1*Math.pow(gst,2)/ tst)* Math.pow(gst,5) * dgst);
                gst = gst + dgst;
                if (sum>gstt)
                    pgst[i-1] = gst - (dgst / 2.E0);
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
        IAtomContainer molecule2[]= inputFile.retMolecule();



        for(int i=0;i<molecule2.length;i++) {


            Random.rluxgo(3,5013486,0,0);

                Cluster ion =new Cluster(molecule2[i]);
  //              this.collisionCrossSection(ion,'H');
                 return this.collisionCrossSection(ion,'H');

        }
        return "False";
    }

}
