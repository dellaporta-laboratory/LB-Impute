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

public class ImputeMain {
	
	public String method;
	public boolean printhelp = false;
	
	public static void main(String[] args){
		if(args.length == 0){
			System.out.println("LB-Impute ver 1.0 written by Christopher Heffelfinger and Christopher Fragoso");
			System.out.println("Try -help for more information or check out the manual");
			System.out.println("Source code is available upon request.");
			System.exit(0);
		}
		ImputeMain imputemain = new ImputeMain();
		imputemain.start(args);
	}

	public void start(String[] args){
		parseargs(args);
		if(printhelp){
			System.out.println("LB-Impute is a program for imputing missing data and correcting erroneous data in biallelic populations. Command line options are as follows:");
			System.out.println("-method	<impute, compare, randomremove>	Determines what operation LB-Impute will undertake. REQUIRED.\n");
			System.out.println();
			System.out.println();
			System.out.println("IMPUTE OPTIONS\n");	
			System.out.println("-f	<filename>	VCF input file name. REQUIRED.");
			System.out.println("-o	<filename>	Output file name. REQUIRED.");
			System.out.println("-parents	<parentnames>	Names of parental samples separated by comma. NO SPACE BETWEEN COMMA AND NAMES. Names with spaces or commas will not work properly. REQUIRED.");
			System.out.println("-offspringimpute		LB-Impute will impute offspring.");
			System.out.println("-parentimpute		Default. LB-Impute will impute parents.");
			System.out.println("-readerr	<double>	Default value is 0.05. Probability that any given read is erroneous.");
			System.out.println("-genotypeerr	<double>	Default value is 0.05. Probability that any given genotype call is erroneous.");
			System.out.println("-recombdist	<integer>	Default value is 10,000,000. Expected distance of 50cM.");
			System.out.println("-resolveconflicts		Default is off. When on, it will use the path with the highest probability to resolve a conflict rather than leaving a conflicting genotype empty. May reduce accuracy.");
			System.out.println("-dr 		Default is off. When on, a homozygous to homozygous recombination event (double recombination) will have the same probability as a single event. This may be useful for inbred populations such as NAMs.\n");
			System.out.println("-minsamples	<integer>	Default is 5. Minimum number of imputed samples required to infer a parental genotype. Has no function when imputing offspring.\n");
			System.out.println("-minfraction	<double>	Default is 0.5. For a genotype at a given locus to be assigned to a parent, that genotype must be at least this fraction of the total genotypes of offspring assigned to that parent. This value should be somewhat low for low coverage populations, as false homozygotes will be confounding.");
			System.out.println("-window	<integer>	Default is 7. Window length of the Vitterbi trellis. Longer windows produce more accurate results, but have longer runtimes. It is recommended to keep this between 5 and 7, though 2 is the minimum.");
			System.out.println();	
			System.out.println();
			System.out.println();
			System.out.println("RANDOMREMOVE OPTIONS\n");
			System.out.println("Random remove removes calls at a specific level or range of coverage to serve as a validation set. When used along with compare, it can validate imputation on a dataset. Prints new dataset to stdout.");
			System.out.println();
			System.out.println("-f	<filename>	Name of VCF file. REQUIRED.");
			System.out.println("-parents	<parentnames>	Names of parental samples separated by comma. NO SPACE BETWEEN COMMA AND NAMES. Names with spaces or commas will not work properly. REQUIRED.");
			System.out.println("-removefraction	<double>	Fraction of calls fitting criteria to be removed.");
			System.out.println("-minhomcov	<int>	Minimum coverage of homozygous calls to be removed. Inclusive. REQUIRED.");
			System.out.println("-maxhomcov	<int>	Maximum coverage of homozygous calls to be removed. Inclusive. REQUIRED.");
			System.out.println("-minhetcov	<int>	Minimum coverage of heterozygous calls to be removed. Inclusive. REQUIRED.");
			System.out.println("-maxhetcov	<int>	Maximum coverage of heterozygous calls to be removed. Inclusive. REQUIRED.");
			System.out.println();
			System.out.println();
			System.out.println("COMPARE OPTIONS");
			System.out.println("Compare compares the original file with an imputed file and a file with a validation set removed then outputs statistics to stdout.");		
			System.out.println("-mainfile	<filename>	Name of original file. REQUIRED.");
			System.out.println("-imputefile	<filename>	Name of imputed file. REQUIRED.");
			System.out.println("-missfile	<filename>	Name of file with validation set removed. REQUIRED.");
			System.out.println("-offspringcompare		Compares offspring rather than parental calls.");
			System.out.println("-parents	<parentnames>	Names of parental samples separated by comma. NO SPACE BETWEEN COMMA AND NAMES. Names with spaces or commas will not work properly. REQUIRED.");
			System.exit(0);
		}
		if(method.equals("compare")){ // this is the algorithm for accuracy metrics. Requires a file with correct genotypes, erroneous or missing genotypes, and imputed genotypes.
			Compare compare = new Compare();
			compare.compare(args);
		}
		if(method.equals("randomremove")){ // I don't think this is used anymore.
			RemoveValues randomremove = new RemoveValues();
			randomremove.remove(args);
		}
		if(method.equals("impute")){ // this is the primary algorithm for imputation.
			Impute impute = new Impute();
			impute.impute(args);
		}
	}
	
	public void parseargs(String[] args){
		boolean methodentered = false;
		for(int i = 0; i < args.length; i++){
			if(args[i].equals("-method")){ // method options are compare, randomremove, and impute
				i++;
				method = args[i];
				methodentered = true;
			}
			if(args[i].equals("-help")){
				printhelp = true;
			}
		}
		if(!methodentered){
			System.out.println("Please select a method to use via -method. See -help or the manual for more information");
			System.exit(0);
		}
	}
}
