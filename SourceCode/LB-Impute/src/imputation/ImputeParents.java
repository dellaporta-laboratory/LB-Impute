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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ImputeParents {
	
	
	public int states;
	public int startcol;
	private int[][] parmap;
	private double[][] probpath;
	private int[] distmap;
	private int genotypecheck = 0;
	private int coveragecheck = 2;
	private int readcheck = 1;
	public int markovorder = 6;
	public double err = 0.01;
	public int recomb = 1000000;
	public int minsamples;
	public double minfraction;
	public boolean assumebiallele = false;
	public boolean resolveconflictcallswithmaxprob = false;
	public int minmarkers = 0;
	public double maxprob = 1;
	public double minprob = 0;
	public boolean drp = true;
	
	public ImputeParents (int mystates, int mystartcol, int varcount, int samplecount, int mygenotypecheck, int myreadcheck, int mycoveragecheck, double myerr, int myrecomb, int mymarkovorder, int myminsamples, double myminfraction, boolean myassumebiallele, boolean myresolve, int myminmarkers, double mygenotypeerror, boolean mydrp){
		states = mystates;
		startcol = mystartcol;
		parmap = new int[varcount][mystates];
		probpath = new double[varcount][mystates];
		distmap = new int[varcount];
		genotypecheck = mygenotypecheck;
		coveragecheck = mycoveragecheck;
		markovorder = mymarkovorder;
		minsamples = myminsamples;
		minfraction = myminfraction;
		assumebiallele = myassumebiallele;
		readcheck = myreadcheck;
		err = myerr;
		recomb = myrecomb;
		resolveconflictcallswithmaxprob = myresolve;
		minmarkers = myminmarkers;
		if(minmarkers < markovorder){
			minmarkers = markovorder;
		}
		maxprob = 1 - (2 * mygenotypeerror);
		minprob = mygenotypeerror;
		drp = mydrp;
	}
	
	public List<String> imputeparents(List<String> variants, Map<Integer, String> parents){
		List<String> returnvalues = new ArrayList<String>();
		retrieveparentals(variants, parents);
		int total = variants.get(0).split("\t").length;
		String[][] allsamplesmap = new String[variants.size()][total - startcol + 1];
		getdistances(variants);
		for(int i = startcol; i < total; i++){
			int samplecount = i - startcol + 1;
			int totalcount = total - startcol;
			System.out.println("Imputing sample " + samplecount + " of " + totalcount);
			List<String> samplestates = new ArrayList<String>();
			if(!parents.containsKey(i)){
				getprobabilities2(variants, i);
				samplestates = makestatemap(variants, i);
				for(int j = 0; j < samplestates.size(); j++){
					allsamplesmap[j][i - startcol] = samplestates.get(j);
				}
				samplestates.clear();
			} else {
				//System.out.println("Parental sample");
			}
		}
		returnvalues = imputeparentstates(variants, allsamplesmap, parents);
		return returnvalues;
	}
	
	public List<String> imputeparentstates (List<String> variants, String[][] allsamplemap, Map<Integer, String> parents){
		List<String> returnvalues = new ArrayList<String>();
		int[][] probabilities = new int[states][4];
		String parseline = variants.get(2).split("\t")[8];
		String[] parselineparts = parseline.split(":");
		int correct = 0;
		int incorrect = 0;
		for(int f = 0; f < parselineparts.length; f++){
			if(parselineparts[f].equals("GT")){genotypecheck = f;}
		}
		for(int f = 0; f < parselineparts.length; f++){
			if(parselineparts[f].equals("DP")){coveragecheck = f;}
		}
		for(int i = 0; i < variants.size(); i++){
			String[] vparts = variants.get(i).split("\t");
			int sample = 0;
			for(int q = 0; q < states; q++){
				for(int p = 0; p < 4; p++){
					probabilities[q][p] = 0;
				}
			}
			for(int k = startcol; k < vparts.length; k++){
				if(!vparts[k].equals("./.") && !parents.containsKey(k)){
					String[] vcomp = vparts[k].split(":");
					int coverage = Integer.parseInt(vcomp[coveragecheck]);
					if(coverage > 0 && !vcomp[genotypecheck].equals("./.")){
						String[] genotypes = vcomp[genotypecheck].split("/");
						if(!allsamplemap[i][sample].equals("X")){
							int state = Integer.parseInt(allsamplemap[i][sample]);
							if(vcomp[genotypecheck].equals("0/0")){probabilities[state][0]++;}
							if(vcomp[genotypecheck].equals("1/1")){probabilities[state][1]++;}
							if(vcomp[genotypecheck].equals("2/2")){probabilities[state][2]++;}
							if(!genotypes[0].equals(genotypes[1])){probabilities[state][3]++;}
						}
					}
				}
				sample++;
			}
			Map<Integer, String> results = new HashMap<Integer, String>();
			for(int l = 0; l < states; l++){
				double[] probs = new double[4];
				int total = 0;
				for(int m = 0; m < 4; m++){
					total = total + probabilities[l][m];
				}
				for(int m = 0; m < 4; m++){
					if(total > 0){
						probs[m] = (double)probabilities[l][m]/(double)total;
					}
					else {
						probs[m] = 0;
					}
				}
				double maxprob = 0;
				int maxgeno = 999;
				for(int m = 0; m < 4; m++){
					if(probs[m] > maxprob){
						maxprob = probs[m];
						maxgeno = m;
					}
				}
				
				if(total >= minsamples && maxprob >= minfraction){
					if(maxgeno == 0){results.put(l, "0/0");}
					if(maxgeno == 1){results.put(l, "1/1");}
					if(maxgeno == 2){results.put(l, "2/2");}
					if(maxgeno == 3){results.put(l, "X");}
				} else {
					results.put(l, "X");
				}
				//System.out.println("State\t" + l + "\t" + total + "\t" + maxprob + "\t" + maxgeno);
			}
			if(assumebiallele){
				int misscount = 0;
				int misspar = 0;
				String presentallele = "None";
				for(int l = 0; l < 2; l++){
					if(results.get(l).equals("X")){
						misscount++;
						misspar = l;
					}
					if(!results.get(l).equals("X")){
						presentallele = results.get(l);
					}
				}
				if(misscount == 1 && !presentallele.equals("None")){
					String replaceallele = new String();
					if(presentallele.equals("0/0")){
						replaceallele = "1/1";
					}
					if(presentallele.equals("1/1")){
						replaceallele = "0/0";
					}
					if(presentallele.equals("2/2")){
						replaceallele = "X";
					}
					results.put(misspar, replaceallele);
				}
			}
			int parent = 0;
			String newline = vparts[0];
			for(int k = 1; k < startcol; k++){
				newline += "\t" + vparts[k];
			}
			for(int k = startcol; k < vparts.length; k++){
				if(parents.containsKey(k)){
					if(!vparts[k].equals("./.")){
						String[] vcomp = vparts[k].split(":");
						int coverage = Integer.parseInt(vcomp[coveragecheck]);
						if(coverage > 0 && !vcomp[genotypecheck].equals("./.")){
							if(!results.get(parent).equals("X")){
								String actualcall = vcomp[genotypecheck];
								String predictedcall = results.get(parent);
								if(predictedcall.equals(actualcall)){correct++;}
								if(!predictedcall.equals(actualcall)){incorrect++;}
								//System.out.println(predictedcall + "\t" + actualcall);
							}
							newline += "\t" + vparts[k];
						} else {
							if(results.get(parent).equals("X")){ 
								newline += "\t" + vparts[k];
								} else {
									newline += "\t" + results.get(parent) + ":99,99:99:99:99";
								}
 						}
					} else {
						if(results.get(parent).equals("X")){
							newline += "\t" + vparts[k];
							} else {
								newline += "\t" + results.get(parent) + ":99,99:99:99:99";
							}
					}
				parent++;
				} else {
					newline += "\t" + vparts[k];
				}
			}
			//System.out.println(newline);
			returnvalues.add(newline);
		}
		int tot = correct + incorrect;
		//System.out.println(correct + "\t" + incorrect + "\t" + tot);
		return returnvalues;
	}

	
	public void getdistances(List<String> variants){
		for(int i = 0; i < variants.size(); i++){
			distmap[i] = Integer.parseInt(variants.get(i).split("\t")[1]);
		}
	}
	
	public List<String> makestatemap(List<String> variants, int sample){
		List<String> collapsemap = new ArrayList<String>();
		List<String> collapsevalues = new ArrayList<String>();
		List<List<String>> collapsereturns = new ArrayList<List<String>>();
		CollapseValuesOffspring collapse = new CollapseValuesOffspring();
		collapsereturns = collapse.collapseprobs(probpath);
		collapsevalues = collapsereturns.get(0);
		collapsemap = collapsereturns.get(1);
		List<String> combinedstates = new ArrayList<String>();
		if(collapsevalues.size() > 0){ // if no called markers are present on chromosome, this step ensures the chromosome is omitted.
			ParseDistances parsedists = new ParseDistances(collapsemap, distmap, collapsevalues.size() - 1);
			int[] distances = new int[collapsevalues.size() - 1];
			distances = parsedists.finddistances();
			FindPath2 forward = new FindPath2(states, markovorder, collapsevalues, distances, recomb, "forward", minmarkers, drp);
			FindPath2 reverse = new FindPath2(states, markovorder, collapsevalues, distances, recomb, "reverse", minmarkers, drp);
			List<String> forwardstates = new ArrayList<String>();
			List<List<String>> viterbireturnsforward = new ArrayList<List<String>>();
			List<List<String>> viterbireturnsreverse = new ArrayList<List<String>>();
			viterbireturnsforward = forward.viterbi();
			List<String> reversestates = new ArrayList<String>();
			viterbireturnsreverse = reverse.viterbi();
			forwardstates = viterbireturnsforward.get(0);
			reversestates = viterbireturnsreverse.get(0);
			List<String> forwardprobs = new ArrayList<String>();
			List<String> reverseprobs = new ArrayList<String>();
			forwardprobs = viterbireturnsforward.get(1);
			reverseprobs = viterbireturnsreverse.get(1);
			if(!resolveconflictcallswithmaxprob){
				combinedstates = combinestates(forwardstates, reversestates);
			}
			if(resolveconflictcallswithmaxprob){
				combinedstates = combinestatesresolve(forwardstates, reversestates, forwardprobs, reverseprobs);
			}
			//String printline = new String();
			//for(String c: combinedstates){
			//	printline += c;
			//}
			//System.out.println(printline);
		}
		List<String> finalmap = new ArrayList<String>();
		finalmap = deconvolute(combinedstates, collapsemap);
		return finalmap;
	}
	
	public List<String> deconvolute(List<String> combinedstates, List<String> collapsemap){
		List<String> returnvalues = new ArrayList<String>();
		String oldstate = "X";
		boolean remove = true;
		int statepos = 0;
		int removecount = 0;
		for(int i = 0; i < collapsemap.size(); i++){
			String call = collapsemap.get(i).split("\t")[1];
			if(call.equals("Keep")){
				String combinedstate = combinedstates.get(statepos);
				if(remove && removecount > 0){
					if(oldstate.equals(combinedstates.get(statepos))){
						for(int k = 0; k < removecount; k++){
							returnvalues.add(oldstate);
						}
					} else {
						for(int k = 0; k < removecount; k++){
							returnvalues.add("X");
						}
					}
				}
				oldstate = combinedstates.get(statepos);
				returnvalues.add(oldstate);
				statepos++;
				removecount = 0;
				remove = false;
			}
			if(call.equals("Remove")){
				removecount++;
				remove = true;
			}
		}
		if(remove){
			for(int k = 0; k < removecount; k++){
				returnvalues.add("X");
			}
		}
		return returnvalues;
	}

	public List<String> combinestatesresolve(List<String> forwardstates, List<String> reversestates, List<String> forwardprobs, List<String> reverseprobs){
		List<String> returnvalues = new ArrayList<String>();
		for(int i = 0; i < forwardstates.size(); i++){
			int j = reversestates.size() - 1 - i;
			if(forwardstates.get(i).equals(reversestates.get(j))){
				returnvalues.add(forwardstates.get(i));
			} else { 
				double forwardprob = Double.parseDouble(forwardprobs.get(i));
				double reverseprob = Double.parseDouble(reverseprobs.get(j));
				if(forwardprob > reverseprob){returnvalues.add(forwardstates.get(i));}
				if(reverseprob > forwardprob){returnvalues.add(reversestates.get(j));}
				if(forwardprob == reverseprob){returnvalues.add("X");}
			}
		}
		return returnvalues;
	}
	
	public List<String> combinestates(List<String> forwardstates, List<String> reversestates){
		List<String> returnvalues = new ArrayList<String>();
		for(int i = 0; i < forwardstates.size(); i++){
			int j = reversestates.size() - 1 - i;
			if(forwardstates.get(i).equals(reversestates.get(j))){
				returnvalues.add(forwardstates.get(i));
			} else { 
				returnvalues.add("X");
			}
		}
		return returnvalues;
	}
	
	
	public void getprobabilities2(List<String> variants, int sample){
		for(int k = 0; k < variants.size(); k++){
			String checkvar = variants.get(k).split("\t")[sample];
			if(checkvar.equals("./.")){
				for(int j = 0; j < states; j++){
					probpath[k][j] = 1;
				}
			}
			if(!checkvar.equals("./.")){
				if(checkvar.split(":")[genotypecheck].equals("./.") | checkvar.split(":")[coveragecheck].equals("0")){
					for(int j = 0; j < states; j++){
						probpath[k][j] = 1;
					}
				} else {
					int refallele = Integer.parseInt(checkvar.split(":")[readcheck].split(",")[0]);
					int nonrefallele = Integer.parseInt(checkvar.split(":")[readcheck].split(",")[1]);
					double check1 = Math.pow((1 - err), nonrefallele) * Math.pow(err, refallele);
					double check2 = Math.pow((1 - err), refallele) * Math.pow(err, nonrefallele);
					double check3 = Math.pow(0.5, refallele) * Math.pow(0.5, nonrefallele);
					double maxcheck = 0;
					if(check1 >= check2 && check1 >= check3){maxcheck = check1;}
					if(check2 >= check1 && check2 >= check3){maxcheck = check2;}
					if(check3 >= check1 && check3 >= check2){maxcheck = check3;}
					double rawvalue = 0;
					double normvalue = 0;
					for(int j = 0; j < states - 1; j++){
						if(parmap[k][j] == 1){
							rawvalue = Math.pow((1 - err), nonrefallele) * Math.pow(err, refallele);
							normvalue = rawvalue/maxcheck * maxprob + minprob;
							probpath[k][j] = normvalue;
						}
						if(parmap[k][j] == 0){
							rawvalue = Math.pow((1 - err), refallele) * Math.pow(err, nonrefallele);
							normvalue = rawvalue/maxcheck * maxprob + minprob;
							probpath[k][j] = normvalue;
						}
						if(parmap[k][j] == 999){
							if(check1 > check2){rawvalue = check1;}
							if(check2 >= check1){rawvalue = check2;}
							normvalue = rawvalue/maxcheck * maxprob + minprob;
							probpath[k][j] = normvalue;
						}
					}
					rawvalue = Math.pow(0.5, refallele) * Math.pow(0.5, nonrefallele);
					normvalue = rawvalue/maxcheck * maxprob + minprob;
					probpath[k][states - 1] = normvalue;
				}
			}
		}
	}
	
	public void retrieveparentals(List<String> variants, Map<Integer, String> parents){
		for(int i = 0; i < variants.size(); i++){
			String[] vparts = variants.get(i).split("\t");
			int parent = 0;
			for(int k = startcol; k < vparts.length; k++){
				if(parents.containsKey(k)){
					if(!vparts[k].equals("./.")){
						if(!vparts[k].split(":")[genotypecheck].equals("./.") && !vparts[k].split(":")[coveragecheck].equals("0")){
							parmap[i][parent] = Integer.parseInt(vparts[k].split(":")[genotypecheck].split("/")[1]);
						} else {
							parmap[i][parent] = 999;
						}
					} else {
						parmap[i][parent] = 999;
					}
					parent++;
				}
			}	
		}
	}
	

	
}
