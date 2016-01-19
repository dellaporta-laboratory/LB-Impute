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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FindPath2 {
	
	private int markovorder;
	private int states;
	private double[][] probs;
	private boolean forward = true;
	private int[] dists;
	private int recomb;
	private int minmarkers = 0;
	private boolean drp = true;
 	
	public FindPath2 (int mystates, int mymarkovorder, List<String> probabilities, int[] distances, int myrecomb, String order, int myminmarkers, boolean mydrp){ //constructor. adjusts path to be imputed based on forward or reverse.
		markovorder = mymarkovorder;
		states = mystates;
		recomb = myrecomb;
		drp = mydrp;
		if(order.equals("forward")){forward = true;}
		if(order.equals("reverse")){forward = false;}
		probs = new double[probabilities.size()][states];
		dists = new int[distances.length];
		minmarkers = myminmarkers;
		if(forward){
			dists = distances;
			for(int i = 0; i < probabilities.size(); i++){
				String[] probparts = probabilities.get(i).split("\t");
				for(int j = 0; j < probparts.length; j++){
					probs[i][j] = Double.parseDouble(probparts[j]);
				}
			}
		}
		if(!forward){
			int h = 0;
			for(int i = distances.length - 1; i >= 0; i = i - 1){
				dists[h] = distances[i];
				h++;
			}
			h = 0;
			for(int i = probabilities.size() - 1; i >= 0; i = i - 1){
				String[] probparts = probabilities.get(i).split("\t");
				for(int j = 0; j < probparts.length; j++){
					probs[h][j] = Double.parseDouble(probparts[j]);
				}
				h++;
			}
		}
	}
	
	
	public List<List<String>> viterbi(){
		List<List<String>> returns = new ArrayList<List<String>>();
		List<String> returnvalues = new ArrayList<String>();
		boolean first = true;
		Map<String, Double> probpaths = new HashMap<String, Double>();
		Map<String, Double> oldprobpaths = new HashMap<String, Double>();
		List<String> probabilities = new ArrayList<String>();
		//System.out.println(probs.length);
		if(probs.length > minmarkers){ //if too few markers are present on a chromosome for a sample, returns a set of empty values.
			for(int i = 0; i < probs.length - 1; i++){ //iterates from first to last value
				if(!first){
					String lastvalue = returnvalues.get(returnvalues.size() - 1);
					probpaths.clear();
					probpaths.put(lastvalue + ",", 1.0);
					for(int k = i + 1; k < i + 1 + markovorder + 1 && k < probs.length; k++){ //builds the viterbi trellis for each value from i to i + markovorder + 1. So if the markovorder is 4, markers 0 to 5 are incorporated for a total trellis length of 6.
						oldprobpaths.clear();
						oldprobpaths.putAll(probpaths);
						probpaths.clear();
						for(int l = 0; l < probs[k].length; l++){
							for(String key: oldprobpaths.keySet()){ //for each path in the trellis, all possible subsequent paths are calculated and a new trellis regenerated. transition probabilities are determined based on distance. emission probabilities are preexisting based on allelic depth of coverage
								String newkey = null;
								double newscore = 0;
								int lastpos = Integer.parseInt(key.split(",")[key.split(",").length - 1]);
								if(lastpos == l && lastpos != probs[i].length - 1 && l != probs[i].length - 1){ //homozygous no recombination
									newscore = (0.5 * (1 + Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb))) * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								if(lastpos != l && lastpos != probs[i].length - 1 && l != probs[i].length - 1){ //from homozygous to homozygous state
									newscore = Math.pow(0.5 * (1 - Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb)), 2) * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								if(lastpos != l && lastpos != probs[i].length - 1 && l == probs[i].length - 1){ //from homozygous to heterozygous state
									newscore = (0.5 * (1 - Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb))) * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								if(lastpos == l && lastpos == probs[i].length - 1 && l == probs[i].length - 1){ // no recobination homozygous
									newscore = (0.5 * (1 + Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb)))  * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								if(lastpos != l && lastpos == probs[i].length - 1 && l != probs[i].length - 1){ // heterozygous to homozygous recombination
									newscore = (0.5 * (1 - Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb))) * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								probpaths.put(newkey, newscore);
							}		
						}
					}
					double maxval = 0;
					for(String key: probpaths.keySet()){
						double check = probpaths.get(key);
						if(check > maxval){maxval = check;}
					}
					String bestpath = new String();
					double totalprob = 0;
					for(String key: probpaths.keySet()){
						double check = probpaths.get(key);
						if(check == maxval){
							bestpath = key;
						}
						totalprob += check;
					}
					probabilities.add(Double.toString(maxval));
					String addvalue = bestpath.split(",")[1];
					returnvalues.add(addvalue);
				}
				if(first){ //this section is only to determine the first position in the chain. It is required to be a bit different as there is no true 0 marker to affect the probability of marker 1
					for(int j = 0; j < probs[0].length; j++){
						probpaths.put(Integer.toString(j) + ",", probs[0][j]);
					}
					for(int k = 1; k < 1 + markovorder + 1; k++){
						oldprobpaths.clear();
						oldprobpaths.putAll(probpaths);
						probpaths.clear();
						for(int l = 0; l < probs[k].length; l++){
							for(String key: oldprobpaths.keySet()){
								String newkey = null;
								double newscore = 0;
								int lastpos = Integer.parseInt(key.split(",")[key.split(",").length - 1]);
								if(lastpos == l && lastpos != probs[i].length - 1 && l != probs[i].length - 1){ //no recombination, homozygous
									newscore = (0.5 * (1 + Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb)))  * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								if(lastpos != l && lastpos != probs[i].length - 1 && l != probs[i].length - 1 && drp){ //recombination, homozygous to homozygous
									newscore = Math.pow(0.5 * (1 - Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb)), 2)  * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								if(lastpos != l && lastpos != probs[i].length - 1 && l != probs[i].length - 1 && !drp){ //recombination, homozygous to homozygous
									newscore = (0.5 * (1 - Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb)))  * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								if(lastpos != l && lastpos != probs[i].length - 1 && l == probs[i].length - 1){ //recombination, homozygous to heterozygous
									newscore = (0.5 * (1 - Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb)))  * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								if(lastpos == l && lastpos == probs[i].length - 1 && l == probs[i].length - 1){ //no recombination, heterozygous
									newscore = (0.5 * (1 + Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb)))  * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								if(lastpos != l && lastpos == probs[i].length - 1 && l != probs[i].length - 1){ // recombination, heterozygous to homozygous
									newscore = (0.5 * (1 + Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb)))  * probs[k][l] * oldprobpaths.get(key);
									newkey = key + Integer.toString(l) + ",";
								}
								//System.out.println(dists[k - 1] + "\t" + recomb + "\t" + -(double)dists[k - 1]/(double)recomb + "\t" + Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb) + "\t" + (0.5 * (1 - Math.pow(Math.E,-(double)dists[k - 1]/(double)recomb))));
								probpaths.put(newkey, newscore);
							}		
						}	
					}
					List<Double> sortvalues = new ArrayList<Double>();
					for(String key: probpaths.keySet()){ //once the trellis is complete, the probability of each path is added to a list.
						double check = probpaths.get(key);
						sortvalues.add(check);
					}
					Collections.sort(sortvalues); 
					Collections.reverse(sortvalues);
					double maxval = sortvalues.get(0); //the highest probability is retrieved 
					String bestpath = new String();
					for(String key: probpaths.keySet()){
						double check = probpaths.get(key); //the path with the highest probability is recovered
						if(check == maxval){
							bestpath = key;
						}
					}
					probabilities.add(Double.toString(maxval)); 
					String addvalue = bestpath.split(",")[0]; //the t + 1 value of the path with the highest probability is added to the final chain.
					returnvalues.add(addvalue);
					first = false;
					i = i - 1;
				}
			} 
		} else {
			for(int i = 0; i < probs.length; i++){ //if there are not enough markers to perform imputation, missing values are added for all positions.
				returnvalues.add("X");
			}
		}
		returns.add(returnvalues);
		returns.add(probabilities);
		return returns;
	}
	
}

