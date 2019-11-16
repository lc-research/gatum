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

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;

/**
 * <h1>DataFile</h1>
 * This class provides functions to read a chemical data file and build a CDK molecule
 *
 *
 * @author  M. Dias
 * @version 1.0
 */
public class DataFile {

    private IAtomContainer[] molecule;

    /**
     * This constructor method is used to read the chemical data file,
     * which is in SDF format.
     * @param fileName This is the file name of the chemical data file
     */
    DataFile(String fileName)
    {
        IteratingSDFReader reader=null;
        File sdfFile = new File(fileName);

        ArrayList<IAtomContainer> moleculeList = new ArrayList<IAtomContainer>();
        try {
            reader = new IteratingSDFReader(new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance());
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        while (reader.hasNext()) {
            moleculeList.add((IAtomContainer)reader.next());
        }

        molecule = moleculeList.toArray(new IAtomContainer[moleculeList.size()]);
    }

    /**
     * This method is used to return a copy of the molecule.
     * @return IAtomContainer This returns a copy of the molecule.
     */
    public IAtomContainer[] retMolecule()
    {
        return molecule;
    }


}
