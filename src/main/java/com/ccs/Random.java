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

/**
 * <h1>Random!</h1>
 * This class provides functions to generate random numbers for the MC integration
 *
 *
 * @author  M. Dias
 * @version 1.0
 */
public final class Random {
    public Random() {
    }


    private static double[] SEEDS=new double[24];
    private static double CARRY=0;
    private static int J24=4,I24=18;
    private static double TWOM24=5.9604644775390625E-08;
    private static int[] NEXT=new int[24];
    private static double TWOM12=0.000244140625;
    private static int IN24=9;
    private static int NSKIP;
    private static int LENV=1;
    private static int LXDFLT=3;
    private static int LUXLEV;
    private static int MAXLEV=4;
    private static int JSEED;
    private static int JSDFLT=314159265;
    private static int INSEED;
    private static boolean NOTYET;
    private static double[] ISEEDS=new double[24];
    private static int ICONS=2147483563;
    private static double ITWO24=Math.pow(2,24);
    private static int IGIGA=1000000000;
    private static int KOUNT=0;
    private static int MKOUNT=0;

    /**
     * Subtract-and-borrow random number generator proposed by
     * Marsaglia and Zaman, implemented by F. James with the name
     * RCARRY in 1991, and later improved by Martin Luescher
     * in 1993 to produce "Luxury Pseudorandom Numbers".
     * @return static double This is the generated random number.
     */
    public static double ranlux()
    {
        double RVEC;
        double UNI = SEEDS[J24] - SEEDS[I24] - CARRY;

        if (UNI < 0) {
            UNI++;
            CARRY = TWOM24;
        } else {
            CARRY = 0;
        }

        SEEDS[I24] = UNI;
        I24 = NEXT[I24];
        J24 = NEXT[J24];

        RVEC = UNI;

        if (UNI < TWOM12) {
            RVEC += TWOM24 * SEEDS[J24];
            if (RVEC == 0) {
                RVEC = TWOM24 * TWOM24;
            }
        }

        IN24++;

        if (IN24 == 24) {
            IN24 = 0;


            KOUNT += NSKIP;
            KOUNT++;
            for (int ISK = 1; ISK <= NSKIP; ISK++) {
                UNI = SEEDS[J24] - SEEDS[I24] - CARRY;
                if (UNI < 0) {
                    UNI++;
                    CARRY = TWOM24;
                } else {
                    CARRY = 0;
                }

                SEEDS[I24] = UNI;
                I24 = NEXT[I24];
                J24 = NEXT[J24];

            }
        }


        KOUNT += LENV;
        if (KOUNT > IGIGA) {
            MKOUNT++;
            KOUNT -= IGIGA;
        }

        return RVEC;
    }

    /**
     * This method is used to initiate the random number generator
     * @param LUX  This defines the luxury levels
     * @param INS  This defines the input seed
     * @param K1,K2  These specify how many numbers were generated since the
     * initialization with LUX and INT
     */


    public static void rluxgo(int LUX,int INS,int K1,int K2) {
        int[] NDSKIP = {0, 24, 73, 199, 365};
        if (LUX < 0) {
            LUXLEV = LXDFLT;
        } else if (LUX < MAXLEV) {
            LUXLEV = LUX;
        } else if (LUX < 24 || LUX > 2000) {
            LUXLEV = MAXLEV;
        } else {
            LUXLEV = LUX;
            for (int ILX = 0; ILX <= MAXLEV; ILX++) {
                if (LUX == NDSKIP[ILX] + 24) {
                    LUXLEV = ILX;
                }
            }
        }

        if (LUXLEV < MAXLEV) {
            NSKIP = NDSKIP[LUXLEV];
        } else {
            NSKIP = LUXLEV - 24;
        }
        IN24 = 0;

        if (INS > 0) {
            JSEED = INS;
        } else {
            JSEED = JSDFLT;
        }


        INSEED = JSEED;
        NOTYET = false;
        TWOM24 = 1;
        int K;
        for (int I = 0; I < 24; I++) {
            TWOM24 = TWOM24 *0.5;
            K = JSEED / 53668;
            JSEED = 40014 * (JSEED - K * 53668) - K * 12211;
            if (JSEED < 0) {
                JSEED = JSEED + ICONS;
            }
            ISEEDS[I] =  JSEED % ITWO24;
        }

        TWOM12 = TWOM24 * 4096;

        for (int I = 0; I < 24; I++) {
            SEEDS[I] = ISEEDS[I] * TWOM24;
            NEXT[I] = I-1;
        }

        NEXT[0] = 23;
        I24 = 23;
        J24 = 9;
        CARRY = 0;
        if (SEEDS[23] == 0) {
            CARRY = TWOM24;
        }

        int INNER;
        double UNI;
        KOUNT = K1;
        MKOUNT = K2;
        if(K1 + K2 != 0)
        {
            for (int IOUTER = 1; IOUTER <= K2 + 1; IOUTER++) {
                INNER = IGIGA;
                if(IOUTER == K2 + 1)
                {
                    INNER = K1;
                }

                for (int ISK = 1; ISK <= INNER; ISK++) {
                    UNI = SEEDS[J24] - SEEDS[I24] - CARRY;

                    if (UNI < 0) {
                        UNI = UNI + 1.0;
                        CARRY = TWOM24;
                    } else {
                        CARRY = 0;
                    }

                    SEEDS[I24] = UNI;
                    I24 = NEXT[I24];
                    J24 = NEXT[J24];
                }

            }

            IN24 = KOUNT% (NSKIP + 24);
            int IZIP,IZIP2;
            if(MKOUNT> 0)
            {
                IZIP = IGIGA% (NSKIP + 24);
                IZIP2 = MKOUNT * IZIP + IN24;
                IN24 = IZIP2% (NSKIP + 24);
            }

            if(IN24> 23) {
                IN24 = 0;
            }
        }

    }

}
