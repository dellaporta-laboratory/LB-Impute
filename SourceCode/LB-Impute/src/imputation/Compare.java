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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class Compare {
	
	public int startcol = 9;
	public String parentlist;
	public String mainfile;
	public String missfile;
	public String imputefile;
	private int genotypecheck = 0;
	public boolean compareparents = true;
	
	public void compare(String[] args){
		parseargs(args);
		
		List<String> missvar = new ArrayList<String>();
		List<String> imputevar = new ArrayList<String>();
		List<String> variants = new ArrayList<String>();
		Map<Integer, String> parents = new HashMap<Integer, String>();
		Map<String, String> parentnames = new HashMap<String, String>();
		parentnames = parseparents();
		FindSampleColumns findsamples = new FindSampleColumns();
		
		try{
			missvar = loadFile(missfile);
			imputevar = loadFile(imputefile);
			variants = loadFile(mainfile);
		} catch (Exception e){
			if (e instanceof IOException){
				throw new RuntimeException(e);
			}
		}
		parents = findsamples.findsamples(variants, parentnames, startcol);
		GetVariantNames getnames = new GetVariantNames();
		Map<Integer, String> skiplist = new HashMap<Integer, String>();
		HashMap<Integer, List<String>> getNameReturn = new HashMap<Integer, List<String>>();
		getNameReturn = getnames.getName(variants, skiplist, startcol);
		int mainstartline = Integer.parseInt(getNameReturn.get(1).get(0));
		getNameReturn = getnames.getName(missvar, skiplist, startcol);
		int missstartline = Integer.parseInt(getNameReturn.get(1).get(0));
		getNameReturn = getnames.getName(imputevar, skiplist, startcol);
		int imputestartline = Integer.parseInt(getNameReturn.get(1).get(0));
		for(int i = 0; i < mainstartline; i++){
			variants.remove(0);
		}
		for(int i = 0; i < missstartline; i++){
			missvar.remove(0);
		}
		for(int i = 0; i < imputestartline; i++){
			imputevar.remove(0);
		}
		String parseline = variants.get(2).split("\t")[8];
		String[] parselineparts = parseline.split(":");
		for(int f = 0; f < parselineparts.length; f++){
			if(parselineparts[f].equals("GT")){genotypecheck = f;}
		}
		comparemain(variants, missvar, imputevar, parents);
	}
	
	public void comparemain (List<String> originalvariants, List<String> subtractedvariants, List<String> imputedvariants, Map<Integer, String> parents){
		int totaloriginalcount = 0;
		int totaloriginalimpute = 0;
		int totaloriginalimputecorrect = 0;
		int misserroneouscount = 0;
		int misshetcount = 0;
		int misshomcount = 0;
		int misserroneousimpute = 0;
		int misserroneousimputecorrect= 0;
		int hetimpute = 0;
		int hetcorrect = 0;
		int homimpute = 0;
		int homcorrect =0;
		int originalhettotal = 0;
		int originalhetimpute = 0;
		int originalhomtotal = 0;
		int originalhomimpute = 0;
		int imputehetcorrect = 0;
		int imputehomcorrect = 0;
		int erroneouscount = 0;
		int erroneouscorrect = 0;
		
		
		if(compareparents){
			for(int i = 0; i < originalvariants.size(); i++){
				String[] oparts = originalvariants.get(i).split("\t");
				String[] sparts = subtractedvariants.get(i).split("\t");
				String[] iparts = imputedvariants.get(i).split("\t");
				for(int k = startcol; k < oparts.length; k++){
					String ocall = oparts[k].split(":")[genotypecheck];
					String icall = new String();
					if(iparts[k].split(":").length > 1){
						icall = iparts[k].split(":")[genotypecheck];
					} else {
						icall = iparts[k];
					}
					String scall = sparts[k].split(":")[genotypecheck];
					if(parents.containsKey(k) && !oparts[k].equals("./.")){
						totaloriginalcount++;
						if(ocall.equals("0/1") | ocall.equals("0/2") | ocall.equals("1/2")){
							originalhettotal++;
							if(!icall.equals("./.")){originalhetimpute++;}
							if(ocall.equals(icall)){imputehetcorrect++;}
						}
						if(ocall.equals("0/0") | ocall.equals("1/1") | ocall.equals("2/2")){
							originalhomtotal++;
							if(!icall.equals("./.")){originalhomimpute++;}
							if(ocall.equals(icall)){imputehomcorrect++;}
						}
						if(!icall.equals("./.")){totaloriginalimpute++;}
						if(ocall.equals(icall)){totaloriginalimputecorrect++;}
						if(icall.equals("0/0") | icall.equals("1/1") | icall.equals("2/2")){
							homimpute++;
							if(ocall.equals(icall)){homcorrect++;}
						}
						if(icall.equals("0/1") | icall.equals("0/2") | icall.equals("1/2")){
							hetimpute++;
							if(ocall.equals(icall)){hetcorrect++;}
						}
						if(!ocall.equals(scall)){
							misserroneouscount++;
							if(ocall.equals("0/1") | ocall.equals("0/2") | ocall.equals("1/2")){misshetcount++;}
							if(ocall.equals("1/1") | ocall.equals("0/0") | ocall.equals("2/2")){misshomcount++;}
							if(!icall.equals("./.")){misserroneousimpute++;}
							if(ocall.equals(icall)){misserroneousimputecorrect++;}
						}
					}
				}
			}
		}
		if(!compareparents){
			for(int i = 0; i < originalvariants.size(); i++){
				String[] oparts = originalvariants.get(i).split("\t");
				String[] sparts = subtractedvariants.get(i).split("\t");
				String[] iparts = imputedvariants.get(i).split("\t");
				for(int k = startcol; k < oparts.length; k++){
					String ocall = oparts[k].split(":")[genotypecheck];
					String icall = new String();
					if(iparts[k].split(":").length > 1){
						icall = iparts[k].split(":")[genotypecheck];
					} else {
						icall = iparts[k];
					}
					String scall = sparts[k].split(":")[genotypecheck];
					if(!parents.containsKey(k) && !oparts[k].equals("./.")){
						totaloriginalcount++;
						if(ocall.equals("0/1") | ocall.equals("0/2") | ocall.equals("1/2")){
							originalhettotal++;
							if(!icall.equals("./.")){originalhetimpute++;}
							if(ocall.equals(icall)){imputehetcorrect++;}
						}
						if(ocall.equals("0/0") | ocall.equals("1/1") | ocall.equals("2/2")){
							originalhomtotal++;
							if(!icall.equals("./.")){originalhomimpute++;}
							if(ocall.equals(icall)){imputehomcorrect++;}
						}
						if(!icall.equals("./.")){totaloriginalimpute++;}
						if(ocall.equals(icall)){totaloriginalimputecorrect++;}
						if(icall.equals("0/0") | icall.equals("1/1") | icall.equals("2/2")){
							homimpute++;
							if(ocall.equals(icall)){homcorrect++;}
						}
						if(icall.equals("0/1")){
							hetimpute++;
							if(ocall.equals(icall)){hetcorrect++;}
						}
						if(!ocall.equals(scall)){
							misserroneouscount++;
							if(!scall.equals("./.")){
								erroneouscount++;
								if(icall.equals(ocall)){erroneouscorrect++;}
							}
							if(ocall.equals("0/1") | ocall.equals("0/2") | ocall.equals("1/2")){misshetcount++;}
							if(ocall.equals("1/1") | ocall.equals("0/0") | ocall.equals("2/2")){misshomcount++;}
							if(!icall.equals("./.")){misserroneousimpute++;}
							if(ocall.equals(icall)){misserroneousimputecorrect++;
							}
						}
					}
				}
			}
		}
		System.out.println(imputefile + "\t" + mainfile + "\t" + missfile);
		System.out.println("Total original markers, original markers called\t" + totaloriginalcount + "\t" + totaloriginalimpute);
		System.out.println("Total original imputed, total original correct\t" + totaloriginalimpute + "\t" + totaloriginalimputecorrect);
		System.out.println("Total missing or erroneous, imputed\t" + misserroneouscount + "\t" + misserroneousimpute);
		System.out.println("Total missing or erroneous and imputed, correct\t" + misserroneousimpute + "\t" + misserroneousimputecorrect);
		//System.out.println("Total erroneous, correct\t" + erroneouscount + "\t" + erroneouscorrect);
		System.out.println("Total original homozygous, imputed\t" + originalhomtotal + "\t" + originalhomimpute);
		System.out.println("Total original heterozygous, imputed\t" + originalhettotal + "\t" + originalhetimpute);
		System.out.println("Total original homozygous and imputed, correct\t" + originalhomimpute + "\t" + imputehomcorrect);
		System.out.println("Total original heterozygous and imputed, correct\t" + originalhetimpute + "\t" + imputehetcorrect);
		System.out.println("Values imputed as homozygous, correct\t" + homimpute + "\t" + homcorrect);
		System.out.println("Missing values imputed as heterozygous, correct\t" + hetimpute + "\t" + hetcorrect);
		
	}
	
	
	
	public Map<String, String> parseparents(){
		Map<String, String> returnvalues = new HashMap<String, String>();
		String[] parentvals = parentlist.split(",");
		for(String p: parentvals){
			returnvalues.put(p, "parent");
		}
		return returnvalues;
	}
	
	public ArrayList<String> loadFile(String fileName) throws IOException { // this subroutine loads the files
		
		FileReader fileReader = new FileReader(fileName);
		BufferedReader bufferedReader = new BufferedReader(fileReader);
		String line = null;		
		ArrayList<String> fileLines = new ArrayList<String>();		
		try {
			while ((line = bufferedReader.readLine()) != null) {
				fileLines.add(line);
			}
		}	finally {
			bufferedReader.close();
		}
		return fileLines;
	}
	
	
	public void parseargs(String[] args){
		for(int i = 0; i < args.length; i++){
			boolean foundarg = false;
			if(args[i].equals("-startcol")){
				i++;
				startcol = Integer.parseInt(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-parents")){
				i++;
				parentlist = args[i];
				foundarg = true;
			}
			if(args[i].equals("-mainfile")){
				i++;
				mainfile = args[i];
				foundarg = true;
			}
			if(args[i].equals("-imputefile")){
				i++;
				imputefile = args[i];
				foundarg = true;
			}
			if(args[i].equals("-missfile")){
				i++;
				missfile = args[i];
				foundarg = true;
			}
			if(args[i].equals("-compareoffspring")){
				compareparents = false;
				foundarg = true;
			}
			if(args[i].equals("-method")){ // method options are compare, randomremove, and impute
				i++;
				foundarg = true;
			}
			if(args[i].equals("-help")){
				foundarg = true;
			}
			if(!foundarg){
				System.out.println(args[i] + " is not a valid option. If " + args[i] + " is the name of a parental sample column, be sure that there is no space between the comma and the sample column name.");
				System.out.println();
				System.out.println("For more information, please use the -help option");
				System.exit(0);
			}
		}
	}
	
}
