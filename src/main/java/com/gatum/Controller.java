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

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.concurrent.ExecutionException;

import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;




@RestController
public class Controller {

    @RequestMapping(value="/ccs", method= RequestMethod.POST)
    public String ccs(@RequestParam("file") MultipartFile file)throws InterruptedException, ExecutionException {

        Util ccs=new Util();
        if (!file.isEmpty()) {

            try {
                byte[] bytes = file.getBytes();
                BufferedOutputStream stream =
                        new BufferedOutputStream(new FileOutputStream(new File("input_file/"+file.getOriginalFilename())));
                stream.write(bytes);
                stream.close();

                return ccs.Response(file.getOriginalFilename());
            } catch (Exception e) {
                return "Failed to upload " + file.getOriginalFilename() + " => " + e.getMessage();
            }
        }

        return "Error";
    }

}
