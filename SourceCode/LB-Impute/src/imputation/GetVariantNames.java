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
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GetVariantNames {
	
	public HashMap<Integer, List<String>> getName (List<String> unfilteredvariants, Map<Integer, String> skip, int namestart) {
		HashMap<Integer, List<String>> returnval = new HashMap<Integer, List<String>>();
		int f = 0;
		int startline = 0;
		int nameline = 0;
		boolean foundstart = false;
		while(!foundstart) {
			String[] variantlineparts = unfilteredvariants.get(f).split("\t");
			if(variantlineparts[0].equals("#CHROM")){
				startline = f + 1;
				nameline = f;
				foundstart = true;
			}
			f++;
			if(f > unfilteredvariants.size()){
				System.out.println("Unable to find line beginning with #CHROM, which is the standard string in vcf files indicating line with variant names. Exiting");
				System.exit(0);
			}
		}
		String[] namelineparts = unfilteredvariants.get(nameline).split("\t");
		List<String> samplenames = new ArrayList<String>();
		for(int i = namestart; i < namelineparts.length; i++){
			samplenames.add(namelineparts[i]);
		}
		List<String> startlineval = new ArrayList<String>();
		startlineval.add(Integer.toString(startline));
		returnval.put(1, startlineval);
		returnval.put(2, samplenames);
		return returnval;
	}
}
