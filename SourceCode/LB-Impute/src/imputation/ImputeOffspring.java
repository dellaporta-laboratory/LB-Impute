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
import java.util.List;
import java.util.Map;



public class ImputeOffspring {
	// These values are just to indicate defaults. They will be changed when the class is instantized.
	public int states;
	public int startcol;
	private int[][] parmap;
	private double[][] probpath;
	private int[] distmap;
	private int genotypecheck = 0;
	private int coveragecheck = 2;
	private int readcheck = 1;
	public int markovorder = 6;
	public double err = 0.05;
	public int recomb = 1000000;
	public boolean resolveconflictcallswithmaxprob = false;
	public int minmarkers = 0;
	public double maxprob = 0.90;
	public double minprob = 0.05;
	public boolean drp = true;
	public boolean keeporiginal = false;
	
	public ImputeOffspring (int mystates, int mystartcol, int varcount, int samplecount, int mygenotypecheck, int myreadcheck, int mycoveragecheck, double myerr, int myrecomb, int mymarkovorder, boolean myresolve, int myminmarkers, double mygenotypeerror, boolean mydrp, boolean mykeeporiginal){
		states = mystates;
		startcol = mystartcol;
		parmap = new int[varcount][mystates];
		probpath = new double[varcount][mystates];
		distmap = new int[varcount];
		genotypecheck = mygenotypecheck;
		coveragecheck = mycoveragecheck;
		markovorder = mymarkovorder;
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
		keeporiginal = mykeeporiginal;
	}
	
	public List<String> imputeoffspring(List<String> variants, Map<Integer, String> parents){ // primary imputation subroutine.
		List<String> returnvalues = new ArrayList<String>();
		retrieveparentals(variants, parents);
		int total = variants.get(0).split("\t").length;
		String[][] allsamplesmap = new String[variants.size()][total - startcol + 1];
		getdistances(variants);
		for(int i = startcol; i < total; i++){
			int samplecount = i - startcol + 1;
			int totalcount = total - startcol;
			List<String> samplestates = new ArrayList<String>();
			if(!parents.containsKey(i)){
				getprobabilities2(variants, i); 
				samplestates = makestatemap(variants, i);
				System.out.println("Imputing sample " + samplecount + " of " + totalcount);
				for(int j = 0; j < samplestates.size(); j++){
					allsamplesmap[j][i - startcol] = samplestates.get(j);
				}
				samplestates.clear();
			} else {
				//System.out.println("Parental sample");
			}
		}
		returnvalues = imputeoffspringstates(variants, allsamplesmap, parents);
		return returnvalues;
	}
	
	public List<String> imputeoffspringstates (List<String> variants, String[][] allsamplemap, Map<Integer, String> parents){
		List<String> returnvalues = new ArrayList<String>();
		String parseline = variants.get(2).split("\t")[8];
		String[] parselineparts = parseline.split(":");	
		for(int f = 0; f < parselineparts.length; f++){
			if(parselineparts[f].equals("GT")){genotypecheck = f;}
		}
		for(int f = 0; f < parselineparts.length; f++){
			if(parselineparts[f].equals("DP")){coveragecheck = f;}
		}
		
		for(int i = 0; i < variants.size(); i++){
			int sample = 0;
			String[] vparts = variants.get(i).split("\t");
			String newline = vparts[0];
			for(int k = 1; k < startcol; k++){
				newline += "\t" + vparts[k];
			}
			for(int k = startcol; k < vparts.length; k++){
				//System.out.println(k);
				//for(int p: parents.keySet()){
					//System.out.println(p);
				//}
				if(!parents.containsKey(k)){
					String samplecall = allsamplemap[i][sample];
					//System.out.println(i + "\t" + sample + "\t" + samplecall);
					if(samplecall.equals("X")){
						if(!keeporiginal){newline += "\t" + "./.";}
						if(keeporiginal){newline += "\t" + variants.get(i).split("\t")[k];}
					}
					if(!samplecall.equals("X")){
						int state = Integer.parseInt(allsamplemap[i][sample]);
						String call = new String();
						if(state == states - 1){
							call = "0/1";
						}
						if(state != states - 1){
							call = parmap[i][state] + "/" + parmap[i][state];
						}
						if(!call.equals("999/999")){
							newline += "\t" + call + ":99,99:99:99:99";
						} else {
							if(!keeporiginal){newline += "\t" + "./.";}
							if(keeporiginal){newline += "\t" + variants.get(i).split("\t")[k];}
						}
					}
				}
				if(parents.containsKey(k)){
					newline += "\t" + vparts[k];
				}
				sample++;
			}
			returnvalues.add(newline);
		}
		return returnvalues;		
	}

	
	public void getdistances(List<String> variants){ //calculates distances between markers
		for(int i = 0; i < variants.size(); i++){
			distmap[i] = Integer.parseInt(variants.get(i).split("\t")[1]);
		}
	}
	
	public List<String> makestatemap(List<String> variants, int sample){ //make state map refers to creating a path of predicted states (IBD genotypes) for a given offspring
		List<String> collapsemap = new ArrayList<String>();
		List<String> collapsevalues = new ArrayList<String>();
		List<List<String>> collapsereturns = new ArrayList<List<String>>();
		CollapseValuesOffspring collapse = new CollapseValuesOffspring(); //collapsing values refers to removal of values corresponding to missing markers from the path to be imputed. 
		collapsereturns = collapse.collapseprobs(probpath);
		collapsevalues = collapsereturns.get(0); //these are the retained emission probabilities, which will be sent to the viterbi algorithm
		collapsemap = collapsereturns.get(1); //this keeps track of which markers were retained and which were removed and will need to be imputed later based on viterbi results.
		List<String> combinedstates = new ArrayList<String>();
		if(collapsevalues.size() > 0){ // omits chromosome if no called markers are present
			ParseDistances parsedists = new ParseDistances(collapsemap, distmap, collapsevalues.size() - 1);
			int[] distances = new int[collapsevalues.size() - 1]; //creates an array to hold distances between retained markers. Size of this array is one less than the total number of retained markers.
			distances = parsedists.finddistances(); //calculates distances between retained markers, since the original set of distances was for the total set of markers. The viterbi algorithm should only calculate transition probabilities based on distances between retained though.
			FindPath2 forward = new FindPath2(states, markovorder, collapsevalues, distances, recomb, "forward", minmarkers, drp); //instantizes the forward viterbi path.
			FindPath2 reverse = new FindPath2(states, markovorder, collapsevalues, distances, recomb, "reverse", minmarkers, drp); //reverse viterbi path - the array is read last to first value.
			List<String> forwardstates = new ArrayList<String>();
			List<List<String>> viterbireturnsforward = new ArrayList<List<String>>();
			List<List<String>> viterbireturnsreverse = new ArrayList<List<String>>();
			viterbireturnsforward = forward.viterbi(); //forward path
			List<String> reversestates = new ArrayList<String>();
			viterbireturnsreverse = reverse.viterbi(); // reverse path
			forwardstates = viterbireturnsforward.get(0);
			reversestates = viterbireturnsreverse.get(0);
			List<String> forwardprobs = new ArrayList<String>();
			List<String> reverseprobs = new ArrayList<String>();
			forwardprobs = viterbireturnsforward.get(1);
			reverseprobs = viterbireturnsreverse.get(1);
			if(!resolveconflictcallswithmaxprob){ // this results in the final statemap having the consensus set of calls of the forward and reverse algorithms
				combinedstates = combinestates(forwardstates, reversestates);
			}
			if(resolveconflictcallswithmaxprob){ // this results in the final statemap having the consensus + higher probability set of calls from the forward and reverse algorithms
				combinedstates = combinestatesresolve(forwardstates, reversestates, forwardprobs, reverseprobs);
			}
			String printline = new String();
			for(String c: combinedstates){
				printline += c;
			}
			System.out.println(printline); // just for validation of final imputed genotypes
		} 
		List<String> finalmap = new ArrayList<String>();
		finalmap = deconvolute(combinedstates, collapsemap); // resolves missing data using proximal markers, or leaves data missing if proximal markers conflict
		return finalmap;
	}
	
	public List<String> deconvolute(List<String> combinedstates, List<String> collapsemap){ //after imputation has finished, values that were removed from the imputed path due to being missing need to be imputed. The immediately proximal calls are used for this. In the case where the proximal calls are discordant, the missing call is left missing.
		List<String> returnvalues = new ArrayList<String>();
		String oldstate = "X";
		boolean remove = true;
		int statepos = 0;
		int removecount = 0;
		for(int i = 0; i < collapsemap.size(); i++){
			String call = collapsemap.get(i).split("\t")[1];
			if(call.equals("Keep")){ // Keep means that a genotype call is present.
				if(remove && removecount > 0){
					if(oldstate.equals(combinedstates.get(statepos))){ //checks to see if the upstream and downstream calls from missing values are concordant
						for(int k = 0; k < removecount; k++){ // if concordant, missing values are set to that state.
							returnvalues.add(oldstate);
						}
					} else { // if not concordant, missing values are left missing.
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
			if(call.equals("Remove")){ // counts the number of removed calls so, in the event that the immediate upstream and downstream calls are concordant, the missing values are set to this state
				removecount++;
				remove = true;
			}
		}
		if(remove){
			for(int k = 0; k < removecount; k++){ // Calls at the end of the chromosome are left missing.
				returnvalues.add("X");
			}
		}
		return returnvalues;
	}

	public List<String> combinestatesresolve(List<String> forwardstates, List<String> reversestates, List<String> forwardprobs, List<String> reverseprobs){ //determines the consensus state based on forward and reverse paths. if the paths conflict, the state with the higher probability is added
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
		
	
	
	public List<String> combinestates(List<String> forwardstates, List<String> reversestates){ //this subroutine determines the consensus state based on forward and reverse paths. if the paths conflict, the state is left empty.
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
	
	
	public void getprobabilities2(List<String> variants, int sample){ // this calculates the probabilities of each state based on allelic depth of coverage
		for(int k = 0; k < variants.size(); k++){
			String checkvar = variants.get(k).split("\t")[sample];
			if(checkvar.equals("./.")){ // when a value is missing in that offspring, all states are set to an equal probability (1) at that position.
				for(int j = 0; j < states; j++){
					probpath[k][j] = 1; 
				}
			}
			if(!checkvar.equals("./.")){ //if value is missing in offspring, all states are set to equal probability.
				if(checkvar.split(":")[genotypecheck].equals("./.") | checkvar.split(":")[coveragecheck].equals("0")){
					for(int j = 0; j < states; j++){
						probpath[k][j] = 1;
					}
				} else {
					int refallele = Integer.parseInt(checkvar.split(":")[readcheck].split(",")[0]); //gets the reference allele count
					int nonrefallele = Integer.parseInt(checkvar.split(":")[readcheck].split(",")[1]); //gets the non-reference allele count
					double check1 = Math.pow((1 - err), nonrefallele) * Math.pow(err, refallele); //calculates the three state probabilities. Does not assign them to a particular state yet.
					double check2 = Math.pow((1 - err), refallele) * Math.pow(err, nonrefallele);
					double check3 = Math.pow(0.5, refallele) * Math.pow(0.5, nonrefallele);
					double maxcheck = 0;
					if(check1 >= check2 && check1 >= check3){maxcheck = check1;} //finds the highest probability value.
					if(check2 >= check1 && check2 >= check3){maxcheck = check2;}
					if(check3 >= check1 && check3 >= check2){maxcheck = check3;}
					double rawvalue = 0;
					double normvalue = 0;
					for(int j = 0; j < states - 1; j++){
						if(parmap[k][j] == 1){ //finds the non-reference parent (state)
							rawvalue = Math.pow((1 - err), nonrefallele) * Math.pow(err, refallele); //determines the initial probability of this state.
							normvalue = rawvalue/maxcheck * maxprob + minprob; //normalizes the state to 1, multiplies it by the maximum possible probability - min prob, then adds the min prob.
							probpath[k][j] = normvalue; //adds this normalized value to the final path
						}
						if(parmap[k][j] == 0){ // same thing, just for reference state.
							rawvalue = Math.pow((1 - err), refallele) * Math.pow(err, nonrefallele);
							normvalue = rawvalue/maxcheck * maxprob + minprob; 
							probpath[k][j] = normvalue;
						}
						if(parmap[k][j] == 999){ //if parent is missing, it sets it to the maximum homozygous probability.
							if(check1 > check2){rawvalue = check1;}
							if(check2 >= check1){rawvalue = check2;}
							normvalue = rawvalue/maxcheck * maxprob + minprob;
							probpath[k][j] = normvalue;
						}
					}
					rawvalue = Math.pow(0.5, refallele) * Math.pow(0.5, nonrefallele); //determines heterozygous state probability and normalizes it.
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
