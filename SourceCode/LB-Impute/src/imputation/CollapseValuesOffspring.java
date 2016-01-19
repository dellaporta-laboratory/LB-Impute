package imputation;

//Copyright 2015, Christopher Heffelfinger, Christopher Fragoso, Hongyu Zhao, Stephen Dellaporta
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CollapseValuesOffspring {
	


	public List<List<String>> collapseprobs(double[][] probs){
		List<List<String>> returnvalues = new ArrayList<List<String>>();
		List<String> newprobs = new ArrayList<String>(); //newprobs are emission probabilities for the set of markers that are retained.
		List<String> probmap = new ArrayList<String>(); //the probmap is actually a map indicating whether a given marker in the original set has been kept or removed.
		boolean miss = true;
		for(int i = 0; i < probs.length; i++){
			miss = true;

			for(int k = 1; k < probs[i].length; k++){ //if a set of emission probabilities are all the same at a given marker, that indicates that it is missing or uninformative and should be removed
				if(probs[i][k] != probs[i][k - 1]){
					miss = false;
				}
			}
			if(miss){ 
				probmap.add((Integer.toString(i)) + "\t" + "Remove");
			}
			if(!miss){
				String addstring = Double.toString(probs[i][0]);
				for(int l = 1; l < probs[i].length; l++){
					addstring += "\t" + Double.toString(probs[i][l]);
				}
				probmap.add((Integer.toString(i)) + "\t" + "Keep");
				newprobs.add(addstring);
			}
		}
		returnvalues.add(newprobs); //both lists are returned.
		returnvalues.add(probmap);
		return returnvalues;
	}

}
