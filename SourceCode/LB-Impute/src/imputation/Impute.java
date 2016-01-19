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
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Impute {

	public String variantsfile; //name of file to be imputed
	public String outputfile; //name of the imputed output vcf file.
	public double err = 0.05; //Default per variant error rate. Modifiable with -err flag.
	public int recomb = 1000000; //Default distance until probability of transition to a different state (recombination event) becomes 0.5. Modifiable with -recombdist flag.
	public int markovorder = 6; //Number of markers to "peer forward" when determining the probability of a transition. Actual trellis window is markovorder + 2 (0 marker, 1 marker, and then a number of markers equal to markerorder). Changeable with -window flag
	public boolean offspringimpute = false; //Can impute offspring or parents. Default is parents. -offspringimpute sets this to true.
	public boolean parentimpute = true; //Default imputation of parental genotypes. -parentimpute is left as a flag for clarity, though not necessary.
	public boolean resolveconflicts = false; // Resolving conflicts will, for a marker with differing calls on the forward and reverse path, use the state with the higher probability. -resolveconflicts sets this to true.
	public String parentlist; //parentlist is the list of parental samples. 
	public int genotypecheck = 0; //VCF files use a colon deliminated set of values to provide genotype information for each sample for each marker. This value says which of these values contains the genotype call. This value will be assigned automatically.
	public int coveragecheck = 2; //VCF files use a colon deliminated set of values to provide genotype information for each sample for each marker. This value says which of these values contains the coverage call. This value will be assigned automatically.
	public int startcol = 9; // startcol is the first column of the VCF information to contain genotype calls for samples.
	public int states; // the total number of possible states. should be 3 (parent 1, parent 2, heterozygous). This valuable was originally variable for imputation of MAGIC populations.
	public boolean assumebiallele = false; // When genotypes are assumed, if one parent is successfully imputed and the other is not, the parent which was not initially imputed will receive the call opposite the one that was imputed. -assumebiallelic will set this to true.
	public int minsamples = 5; // This is the minimum number of samples that must have a successful state call at a given site for a parental imputation to be made. Can be changed with -minsamples flag.
	public double minfraction = 0.5; // The minimum fraction of genotype calls assigned to a parent for a given marker for that genotype to be assigned to the parent. Can be changed with -minfraction flag.
	public int minmarkers = 0; // The minimum number of non-missing genotypes on a chromosome for imputation to take place for a given sample. If the sample does not have a sufficient number of markers on a given chromosome, all markers on that chromosome will be set to missing. Will default to the markov order in later steps. Can be changed with -minmarkers flag.
	public double genotypeerror = 0.05; // The expected genotyping error independent of coverage. For instance, misalignments or paralogous artifacts. Maximum emission probability of any state is 1 - genotyping error. Set with -genotypeerr flag.
	public boolean drp = true;
	public boolean keeporiginal = false;
	
	public void impute(String[] args){
		Long starttime = System.currentTimeMillis();
		parseargs(args);
		List<String> variants = new ArrayList<String>();
		// loads vcf file into variants arraylist.
		try{
			variants = loadFile(variantsfile);
		} catch (Exception e){
			if (e instanceof IOException){
				throw new RuntimeException(e);
			}
		}
		// This section identifies the first line of genotype calls and which columns are offspring/which columns are parents by looking at vcf name row.
		Map<Integer, String> parents = new HashMap<Integer, String>();
		Map<Integer, String> skiplist = new HashMap<Integer, String>(); //depriciated
		Map<String, String> parentnames = new HashMap<String, String>();
		parentnames = parseparents();
		GetVariantNames getnames = new GetVariantNames();
		FindSampleColumns findsamples = new FindSampleColumns();
		parents = findsamples.findsamples(variants, parentnames, startcol);
		HashMap<Integer, List<String>> getNameReturn = new HashMap<Integer, List<String>>();
		getNameReturn = getnames.getName(variants, skiplist, startcol);
		int startline = Integer.parseInt(getNameReturn.get(1).get(0));
		
		//removes header, including sample names, from vcf arraylist. Restored in the final output file.
		List<String> header = new ArrayList<String>();
		for(int f = 0; f < startline; f++){
			header.add(variants.get(0));
			variants.remove(0);
		}
		
		//identifies where the genotype call, depth of coverage, and allele specific coverage values are in the colon delimated genotype-calls for each sample.
		String parseline = variants.get(2).split("\t")[8];
		String[] parselineparts = parseline.split(":");
		int readcheck = 1;
		for(int f = 0; f < parselineparts.length; f++){
			if(parselineparts[f].equals("GT")){genotypecheck = f;} // GT = genotype, DP = depth, AD = allele depth. Each of these values identifies which colon separated position these values will be in.
		}
		for(int f = 0; f < parselineparts.length; f++){
			if(parselineparts[f].equals("DP")){coveragecheck = f;}
		}
		for(int f = 0; f < parselineparts.length; f++){
			if(parselineparts[f].equals("AD")){readcheck = f;}
		}
		
		//parent 1 homozyogus, parent 2 homozygous, heterozygous. Possible states for samples. Should be 3.
		states = parents.size() + 1;
		
		//For imputation, the vcf file must now be split by chromosome, and each chromosome imputed independently. To do this, the chromosome for each variant is checked. So long as the chromsome matches the old one, the variant is added to currentchrvariants. When the chromosome changes, currentchrvariants is imputed, returning currentchrimputed, which is added to imputedvariants. Currentchrvariants is then cleared and set to the first variant with the new chromosome.
		String oldchromosome = variants.get(0).split("\t")[0];
		List<String> currentchrvariants = new ArrayList<String>();
		List<String> currentchrimputed = new ArrayList<String>();
		List<String> imputedvariants = header; //imputed variants is the new vcf file with variants imputed. The header of the old vcf file is added to this new one.
		int samplecount = 0; //the number of offspring being imputed. parents are not imputed directly. This value is not zero, but will be determined later.
		for(int i = 0; i < variants.size(); i++){
			if (variants.get(i).split("\t")[0].equals(oldchromosome)){ //if the chromosome does not change, the variant being checked is added to currentchrimputed.
				currentchrvariants.add(variants.get(i));
			}  else { //if the chromsome between two variants does change, then the chromosome is assumed to have ended and imputation takes place.
				samplecount = currentchrvariants.get(0).split("\t").length - parents.size() - startcol + 1;
				int varcount = currentchrvariants.size(); //the total number of markers for the current chromosome being imputed.
				System.out.println("Imputing sequence " + oldchromosome);
				if(offspringimpute){ //offspringimpute is to impute missing offspring.
					ImputeOffspring imputeoffspring = new ImputeOffspring(states, startcol, varcount, samplecount, genotypecheck, readcheck, coveragecheck, err, recomb, markovorder, resolveconflicts, minmarkers, genotypeerror, drp, keeporiginal);
					currentchrimputed = imputeoffspring.imputeoffspring(currentchrvariants, parents);
				}
				if(parentimpute){ //parentimpute is to impute missing parents. 
					ImputeParents imputeparents = new ImputeParents(states, startcol, varcount, samplecount, genotypecheck, readcheck, coveragecheck, err, recomb, markovorder, minsamples, minfraction, assumebiallele, resolveconflicts, minmarkers, genotypeerror, drp);
					currentchrimputed = imputeparents.imputeparents(currentchrvariants, parents);
				}
				imputedvariants.addAll(currentchrimputed);
				oldchromosome = variants.get(i).split("\t")[0];
				currentchrvariants.clear();
				currentchrvariants.add(variants.get(i));
			}
		}
		
		//One final imputation is performed on the final chromosome.
		System.out.println("Imputing sequence " + oldchromosome);
		if(offspringimpute){
			samplecount = currentchrvariants.get(0).split("\t").length - parents.size() - startcol + 1; 
			int varcount = currentchrvariants.size();
			ImputeOffspring imputeoffspring = new ImputeOffspring(states, startcol, varcount, samplecount, genotypecheck, readcheck, coveragecheck, err, recomb, markovorder, resolveconflicts, minmarkers, genotypeerror, drp, keeporiginal);
			currentchrimputed = imputeoffspring.imputeoffspring(currentchrvariants, parents);
		}
		if(parentimpute){
			samplecount = currentchrvariants.get(0).split("\t").length - parents.size() - startcol + 1;
			int varcount = currentchrvariants.size();
			ImputeParents imputeparents = new ImputeParents(states, startcol, varcount, samplecount, genotypecheck, readcheck, coveragecheck, err, recomb, markovorder, minsamples, minfraction, assumebiallele, resolveconflicts, minmarkers, genotypeerror, drp);
			currentchrimputed = imputeparents.imputeparents(currentchrvariants, parents);
		}
		imputedvariants.addAll(currentchrimputed);
		Long runtime = System.currentTimeMillis() - starttime;
		System.out.println(runtime);
		try{
			output(imputedvariants, outputfile);
		} catch (Exception e){
			if(e instanceof IOException){
				throw new RuntimeException(e);
			}
		}
	}
	
	
	public void parseargs(String[] args){ //argument parsing. See global variables at the top for descriptions of what each flag does.
		for(int i = 0; i < args.length; i++){
			boolean foundarg = false;
			if(args[i].equals("-f")){
				i++;
				variantsfile = args[i];
				foundarg = true;
			}
			if(args[i].equals("-o")){
				i++;
				outputfile = args[i];
				foundarg = true;
			}
			if(args[i].equals("-readerr")){
				i++;
				err = Double.parseDouble(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-genotypeerr")){
				i++;
				genotypeerror = Double.parseDouble(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-recombdist")){
				i++;
				recomb = Integer.parseInt(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-window")){
				i++;
				markovorder = Integer.parseInt(args[i]) - 1;
				foundarg = true;
				minmarkers = Integer.parseInt(args[i]) + 1;
			}
			if(args[i].equals("-parentimpute")){
				parentimpute = true;
				offspringimpute = false;
				foundarg = true;
			}
			if(args[i].equals("-offspringimpute")){
				parentimpute = false;
				offspringimpute = true;
				foundarg = true;
			}
			if(args[i].equals("-parents")){
				i++;
				parentlist = args[i];
				foundarg = true;
			}
			if(args[i].equals("-startcol")){
				i++;
				startcol = Integer.parseInt(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-assumebiallelic")){
				assumebiallele = true;
				foundarg = true;
			}
			if(args[i].equals("-minsamples")){
				i++;
				minsamples = Integer.parseInt(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-minfraction")){
				i++;
				minfraction = Double.parseDouble(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-resolveconflicts")){
				resolveconflicts = true;
				foundarg = true;
			}
			if(args[i].equals("-dr")){
				drp = false;
				foundarg = true;
			}
			if(args[i].equals("-keep")){
				keeporiginal = true;
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
	
	public Map<String, String> parseparents(){ //parses the parent names from the initial global String parentlist into a hash. This hash is fed into the FindSamples Class to find which columns contain the parental genotypes in the vcf file.
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
	
	public void output(List<String> outputarray, String filename) throws IOException{ //output subroutine.
		BufferedWriter writer = null;
		writer = new BufferedWriter(new FileWriter(filename));
		for(String line: outputarray){
			writer.write(line + "\n");
		}
		writer.flush();
		writer.close();
	}
	
}
