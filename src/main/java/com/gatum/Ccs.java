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

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import java.util.concurrent.*;

/**
 * <h1>Gatum!</h1>
 * This application calculates the collision cross section of
 * a given molecule using the Trajectory method.
 * <p>
 *
 *
 * @author  M. Dias
 * @version 1.0
 */

//@SpringBootApplication
public class Ccs {
    public static void main(String[] args) throws InterruptedException, ExecutionException {

 //       SpringApplication.run(Ccs.class, args);
        Util ccs=new Util();
        ccs.Response("deprotonated_3D.sdf");


    }
}


