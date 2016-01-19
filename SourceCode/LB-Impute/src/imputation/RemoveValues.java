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
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class RemoveValues {
	
	public String vcffile;
	public int minhomcov = 0;
	public int maxhomcov = 1000;
	public int minhetcov = 0;
	public int maxhetcov = 1000;
	public double removefraction = 0;
	public String parentlist;
	public int startcol = 9;
	public int removecount = 0;
	public int keepcount = 0;
	public int genotypecheck = 0;
	public int coveragecheck = 2;
	
	
	public static void main(String[] args){
		
		RemoveValues remove = new RemoveValues();
		
		remove.remove(args);
	}
	
	public void remove(String[] args){
		
		parseargs(args);
		List<String> variants = new ArrayList<String>();
		
		try{
			variants = loadFile(vcffile);
		} catch(Exception e){
			if(e instanceof IOException){
				throw new RuntimeException(e);
			}
		}
		
		Map<Integer, String> parents = new HashMap<Integer, String>();
		Map<Integer, String> skiplist = new HashMap<Integer, String>(); //depriciated
		Map<String, String> parentnames = new HashMap<String, String>();
		parentnames = parseparents();
		FindSampleColumns findsamples = new FindSampleColumns();
		parents = findsamples.findsamples(variants, parentnames, startcol);
		
		List<String> newvcf = new ArrayList<String>();
		boolean foundstart = false;
		int line = 0;
		int startline = 0;
		while(!foundstart){
			if(variants.get(line).substring(0, 1).equals("#")){
				newvcf.add(variants.get(line));
			}
			if(!variants.get(line).substring(0,1).equals("#")){
				foundstart = true;
				startline = line;
			}
			line++;
		}
		String parseline = variants.get(startline).split("\t")[8];
		String[] parselineparts = parseline.split(":");
		for(int f = 0; f < parselineparts.length; f++){
			if(parselineparts[f].equals("GT")){genotypecheck = f;}
		}
		for(int f = 0; f < parselineparts.length; f++){
			if(parselineparts[f].equals("DP")){coveragecheck = f;}
		}
		SecureRandom random = new SecureRandom();
		for(int i = startline; i < variants.size(); i++){
			boolean removeline = false;
			String newline = variants.get(i).split("\t")[0];
			String[] vparts = variants.get(i).split("\t");
			for(int k = 1; k < startcol; k++){
				newline += "\t" + vparts[k];
			}
			for(int k = startcol; k < vparts.length; k++){
				if(vparts[k].equals("./.") | parents.containsKey(k)){
					newline += "\t" + vparts[k];
				}
				if(!vparts[k].equals("./.") && !parents.containsKey(k)){
					String call = vparts[k].split(":")[genotypecheck];
					int coverage = Integer.parseInt(vparts[k].split(":")[coveragecheck]);
					if(call.equals("./.")){
						newline += "\t./.";
					}
					if(!call.equals("./.")){
						Double r = random.nextDouble();
						if(r > removefraction){
							newline += "\t" + vparts[k];
							keepcount++;
						} else {
							if(!call.split("/")[0].equals(call.split("/")[1]) && coverage >= minhetcov && coverage <= maxhetcov){
								newline += "\t./.";
								removecount++;
							}
							if(call.split("/")[0].equals(call.split("/")[1]) && coverage >= minhomcov && coverage <= maxhomcov){
								newline += "\t./.";
								removecount++;
							}
							if(!call.split("/")[0].equals(call.split("/")[1]) && (coverage < minhetcov | coverage > maxhetcov)){
								newline += "\t" + vparts[k];
								keepcount++;
							}
							if(call.split("/")[0].equals(call.split("/")[1]) && (coverage < minhomcov | coverage > maxhomcov)){
								newline += "\t" + vparts[k];
								keepcount++;
							}
						}
					}
				}
			}
			if(!removeline){newvcf.add(newline);}
		}
		try{
			output(newvcf, "newvcf.vcf");
		} catch (Exception e){
			if(e instanceof IOException){
				throw new RuntimeException(e);
			}
		}
		System.out.println("Removed: " + removecount + "\tkept: " + keepcount);
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
	
	public void parseargs(String[] args){
		
		for(int i = 0; i < args.length; i++){
			boolean foundarg = false;
			if(args[i].equals("-f")){
				i++;
				vcffile = args[i];
				foundarg = true;
			}
			if(args[i].equals("-minhomcov")){
				i++;
				minhomcov = Integer.parseInt(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-maxhomcov")){
				i++;
				maxhomcov = Integer.parseInt(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-minhetcov")){
				i++;
				minhetcov = Integer.parseInt(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-maxhetcov")){
				i++;
				maxhetcov = Integer.parseInt(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-removefraction")){
				i++;
				removefraction = Double.parseDouble(args[i]);
				foundarg = true;
			}
			if(args[i].equals("-parents")){
				i++;
				parentlist = args[i];
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
