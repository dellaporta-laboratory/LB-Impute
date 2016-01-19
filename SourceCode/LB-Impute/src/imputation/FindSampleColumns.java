package imputation;

//Copyright 2015, Christopher Heffelfinger, Christopher Fragoso, Hongyu Zhao, Stephen Dellaporta
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FindSampleColumns { //This subroutine used to do more - now it just finds which row of the VCF file contains the sample names and which is the first row with genotype calls.
	
	public Map<Integer, String> findsamples (List<String> variants, Map<String, String> samplenames, int startcolumn){
		Map<Integer, String> returnvalues = new HashMap<Integer, String>();
		boolean foundnames = false;
		int i = 0;
		while(!foundnames){
			String[] varline = variants.get(i).split("\t");
			if(varline[0].equals("#CHROM")){
				foundnames = true;
				for(int k = startcolumn; k < varline.length; k++){
					if(samplenames.containsKey(varline[k])){
						returnvalues.put(k, varline[k]);
					}
				}
			}
			i++;
		}
		
		return returnvalues;
	}

}
