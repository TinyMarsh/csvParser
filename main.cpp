#include <fstream>
#include <iostream>

#include "main.h"

/*
This is a program which contains a parser for combiPeptData, calculatedRTs and peakAreas csv files outputted by Pescal++

Though really, this program fetches peptides which shared between Mascot and X!Tandem search results and outputs the calculated retention times and peak areas of these peptides
*/

int main(){

	//We will fetch all the sequence, modification position, and charge for each peptide and store them in vectors

	//...said vectors...
	std::vector<std::string> mascotSequences;
	std::vector<std::string> mascotModPos;
	std::vector<std::string> mascotCharges;

	//this function fills said vectors
	getMascot(mascotSequences, mascotModPos, mascotCharges);

	//more of said vectors, but these will hold the X!Tandem peptides
	std::vector<std::string> xTandemSequences;
	std::vector<std::string> xTandemModPos;
	std::vector<std::string> xTandemCharges;

	//this is another function filling said vectors, but this time for X!Tandem results
	getXTandem(xTandemSequences, xTandemModPos, xTandemCharges);

	//We'll store the peptide db IDs of peptides which are identical between each identification software in this handy vector of type pair
	std::vector<std::pair<int, int> > matchedIndexes; 

	//Now we will fill the above handy vector. For every mascot sequence, modpos, and charge, check if there exists an xTandem equivalent. If there is, store the db ID from Mascot and X!Tandem in an int pair
	//While we're here, we might as well output all these matched peptides to a file "sharedPeptides.csv"
	std::ofstream file("sharedPeptides.csv");
	int xTandemIndex=0;
	for(int i=0, size=mascotSequences.size(); i<size; ++i){
		int compare = mascotSequences[i].compare(xTandemSequences[xTandemIndex]);
		while(compare>0 && xTandemIndex<xTandemSequences.size()){
			++xTandemIndex;
			compare = mascotSequences[i].compare(xTandemSequences[xTandemIndex]);
		}
		bool found=0;
		int tempIndex=xTandemIndex;
		while(compare==0&&!found){
			if(mascotModPos[i]==xTandemModPos[tempIndex] && mascotCharges[i]==xTandemCharges[tempIndex]){
				found=1;
				std::pair<int, int> indexes;
				indexes.first = i;
				indexes.second = tempIndex;
				matchedIndexes.push_back(indexes);

				//this is the bit where we output the matched peptides to the file
				file << mascotSequences[i] << "," << mascotModPos[i] << "," << mascotCharges[i] << "\n";
			}
			if(tempIndex<xTandemSequences.size()){
				++tempIndex;
				compare = mascotSequences[i].compare(xTandemSequences[tempIndex]);
			} else{
				break;
			}
		}
	}

	//output all the matched db IDs to a file "indexes.csv"
	std::ofstream indexFile("indexes.csv");
	for(int i=0, size=matchedIndexes.size(); i<size; ++i){
		indexFile << matchedIndexes[i].first << "," << matchedIndexes[i].second << "\n";
	}


	//Now we have a list of indexes of shared peptides, let's get the quant for each
	std::vector<double> mascotQuants;
	std::vector<double> xTandemQuants;

	getMascotQuants(mascotQuants);
	getXTandemQuants(xTandemQuants);

	//output the quant for each matched peptide to a file "comparedQuants.csv"
	std::ofstream quantFile("comparedQuants.csv");
	for(int i=0, size=matchedIndexes.size(); i<size; ++i){
		if(!(mascotQuants[matchedIndexes[i].first]==0 || xTandemQuants[matchedIndexes[i].second]==0)){
			quantFile << matchedIndexes[i].first << "," << mascotQuants[matchedIndexes[i].first] << "," << xTandemQuants[matchedIndexes[i].second] << "\n";
		}
	}

	//Let's get the retention times of each matched peptide
	std::vector<double> mascotRTs;
	std::vector<double> xTandemRTs;

	getMascotRTs(mascotRTs);
	getXTandemRTs(xTandemRTs);

	//And store these retention times to file "comparedRTs.csv"
	std::ofstream RTFile("comparedRTs.csv");
	for(int i=0, size=matchedIndexes.size(); i<size; ++i){
		RTFile << matchedIndexes[i].first << "," << mascotRTs[matchedIndexes[i].first] << "," << xTandemRTs[matchedIndexes[i].second] << "\n";
	}

}

void getMascot(std::vector<std::string> &mascotSequences, std::vector<std::string> &mascotModPos, std::vector<std::string> &mascotCharges){

	std::ifstream file("F:\\data\\CTAM\\analysis\\mascotPosOfMod\\combiPeptData.csv");

	std::string value;

	//first getline just moves the cursor past the first row, since we don't want the headers
	std::getline(file, value, '\n');

	for(int peptide=0; peptide<45187; ++peptide){

		//skip the first five quatation marks so that we get the sequence
		for(int i=0; i<5; ++i){
			std::getline(file, value, '\"');
		}

		//get the sequence
		std::getline(file, value, '\"');
		mascotSequences.push_back(value);

		//skip the next 7 quatation marks to get the charge
		for(int i=0; i<7; ++i){
			std::getline(file, value, '\"');
		}

		//get the charge
		std::getline(file, value, '\"');
		mascotCharges.push_back(value);

		//skip the next quatation mark to get the modPos string
		std::getline(file, value, '\"');

		//get the modPos string
		std::getline(file, value, '\"');
		mascotModPos.push_back(value);

		std::getline(file, value, '\n');
	}
}

void getXTandem(std::vector<std::string> &xTandemSequences, std::vector<std::string> &xTandemModPos, std::vector<std::string> &xTandemCharges){


	std::ifstream file("F:\\data\\CTAM\\analysis\\peptideShakerPosOfMod\\combiPeptData.csv");

	std::string value;

	//first getline just moves the cursor past the first row, since we don't want the headers
	std::getline(file, value, '\n');

	for(int peptide=0; peptide<59410; ++peptide){

		//skip the first five quatation marks so that we get the sequence
		for(int i=0; i<5; ++i){
			std::getline(file, value, '\"');
		}

		//get the sequence
		std::getline(file, value, '\"');
		xTandemSequences.push_back(value);

		//skip the next 7 quatation marks to get the charge
		for(int i=0; i<7; ++i){
			std::getline(file, value, '\"');
		}

		//get the charge
		std::getline(file, value, '\"');
		xTandemCharges.push_back(value);

		//skip the next quatation mark to get the modPos string
		std::getline(file, value, '\"');

		//get the modPos string
		std::getline(file, value, '\"');

		std::size_t found=0;
		while(found!=std::string::npos){
			found = value.find("4", found+1);
			if(found!=std::string::npos){
				value[found] = '3';
			}
		}

		xTandemModPos.push_back(value);

		std::getline(file, value, '\n');
	}
}

void getMascotQuants(std::vector<double> &mascotQuants){
	std::ifstream file("F:\\data\\CTAM\\analysis\\mascotPosOfMod\\peakAreas.csv");

	std::string value;

	//not interested in first line
	std::getline(file, value, '\n');

	for(int i=0; i<45187; ++i){
		for(int j=0; j<352; ++j){
			std::getline(file, value, '\"');
			std::getline(file, value, '\"');
		}

		//get quant value
		std::getline(file, value, '\"');
		std::getline(file, value, '\"');

		mascotQuants.push_back(std::stod(value));

		std::getline(file, value, '\n');
	}
}

void getXTandemQuants(std::vector<double> &xTandemQuants){
	std::ifstream file("F:\\data\\CTAM\\analysis\\peptideShakerPosOfMod\\peakAreas.csv");

	std::string value;

	//not interested in first line
	//std::getline(file, value, '\n');

	for(int i=0; i<59410; ++i){
		for(int j=0; j<70; ++j){
			std::getline(file, value, '\"');
			std::getline(file, value, '\"');
		}

		//get quant value
		std::getline(file, value, '\"');
		std::getline(file, value, '\"');

		xTandemQuants.push_back(std::stod(value));

		std::getline(file, value, '\n');
	}
}

void getMascotRTs(std::vector<double> &mascotRTs){
	std::ifstream file("F:\\data\\CTAM\\analysis\\mascotPosOfMod\\calculatedRTs.csv");

	std::string value;

	//not interested in first line
	std::getline(file, value, '\n');

	for(int i=0; i<45187; ++i){
		for(int j=0; j<352; ++j){
			std::getline(file, value, '\"');
			std::getline(file, value, '\"');
		}

		//get quant value
		std::getline(file, value, '\"');
		std::getline(file, value, '\"');

		mascotRTs.push_back(std::stod(value));

		std::getline(file, value, '\n');
	}
}

void getXTandemRTs(std::vector<double> &xTandemRTs){
	std::ifstream file("F:\\data\\CTAM\\analysis\\peptideShakerPosOfMod\\calculatedRTs.csv");

	std::string value;

	//not interested in first line
	std::getline(file, value, '\n');

	for(int i=0; i<59410; ++i){
		for(int j=0; j<257; ++j){
			std::getline(file, value, '\"');
			std::getline(file, value, '\"');
		}

		//get quant value
		std::getline(file, value, '\"');
		std::getline(file, value, '\"');

		xTandemRTs.push_back(std::stod(value));

		std::getline(file, value, '\n');
	}
}

template<class Iter, class T>
Iter vectorBinarySearch(Iter begin, Iter end, T &val){
    // Finds the lower bound in at most log(last - first) + 1 comparisons
    Iter i = std::lower_bound(begin, end, val);

    if (i != end && !(val < *i))
        return i; // found
    else
    	std::cout << "error: vector binary search found no matches!\n";
        return end; // not found
}