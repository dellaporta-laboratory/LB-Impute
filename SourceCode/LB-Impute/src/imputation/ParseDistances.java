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

import java.util.List;

public class ParseDistances {
	
	private List<String> collapsemap;
	private int[] distmap;
	private int collapsecount;
	
	public ParseDistances(List<String> mycollapsemap, int[] mydistmap, int mycollapsecount){ //requires the list of markers that have been retained or removed and the distances between them
		collapsemap = mycollapsemap;
		distmap = mydistmap;
		collapsecount = mycollapsecount;
	}
	
	public int[] finddistances(){ 
		int k = 0;
		int[] returnvalues = new int[collapsecount];
		boolean firststart = false;
		int oldvalue = 0;
		for(int i = 0; i < distmap.length; i++){
			String collapsecheck = collapsemap.get(i).split("\t")[1];
			//System.out.println(collapsecheck);
			if(firststart){ //if a marker is kept, it looks up the positions of this marker and the kept marker before it. calculates the distance between the two and passes it to the array.
				if(collapsecheck.equals("Keep")){
					int addvalue = distmap[i] - oldvalue;
					//System.out.println(collapsecount + "\t" + addvalue);
					returnvalues[k] = addvalue;
					k++;
					oldvalue = distmap[i];
				}
			}
			if(!firststart && collapsecheck.equals("Keep")){ //finds the first kept marker and sets it at position 0
				firststart = true;
				oldvalue = distmap[i];
			}
			
		}
		
		
		return returnvalues;
	}

}
