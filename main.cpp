#include <fstream>
#include <iostream>

#include "main.h"

int main(){

	//We will fetch all the sequence, modification position, and charge for each peptide and store them in vectors

	std::vector<std::string> mascotSequences;
	std::vector<std::string> mascotModPos;
	std::vector<std::string> mascotCharges;

	getMascot(mascotSequences, mascotModPos, mascotCharges);

	std::vector<std::string> xTandemSequences;
	std::vector<std::string> xTandemModPos;
	std::vector<std::string> xTandemCharges;

	getXTandem(xTandemSequences, xTandemModPos, xTandemCharges);

	//We'll store the peptide db IDs of peptides which are identical between each identification software
	std::vector<std::pair<int, int> > matchedIndexes; 

	//for every mascot sequence, modpos, and charge, check if there exists an xTandem equivalent. If there is, store the db ID from Mascot and X!Tandem in an int pair
	int xTandemIndex=0;
	for(int i=0, size=mascotSequences.size(); i<size; ++i){
		int compare = mascotSequences[i].compare(xTandemSequences[xTandemIndex]);
		while(compare==0){
			if(mascotModPos[i]==xTandemModPos[xTandemIndex] && mascotCharges[i]==xTandemCharges[xTandemIndex]){
				std::pair<int, int> indexes;
				indexes.first = i;
				indexes.second = xTandemIndex;
				matchedIndexes.push_back(indexes);
			}
			if(xTandemIndex<xTandemSequences.size()){
				++xTandemIndex;
				compare = mascotSequences[i].compare(xTandemSequences[xTandemIndex]);
			} else{
				break;
			}
		}

		while(compare>0 && xTandemIndex<xTandemSequences.size()){
			++xTandemIndex;
			compare = mascotSequences[i].compare(xTandemSequences[xTandemIndex]);
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
		if(!(mascotQuants[matchedIndexes[i].first]<10000000 || xTandemQuants[matchedIndexes[i].second]<1000)){
			quantFile << mascotQuants[matchedIndexes[i].first] << "," << xTandemQuants[matchedIndexes[i].second] << "\n";
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
		RTFile << mascotRTs[matchedIndexes[i].first] << "," << xTandemRTs[matchedIndexes[i].second] << "\n";
	}

}

void getMascot(std::vector<std::string> &mascotSequences, std::vector<std::string> &mascotModPos, std::vector<std::string> &mascotCharges){

	std::ifstream file("C:\\data\\CTAM\\Pescal++\\mascot\\combiPeptData.csv");

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


	std::ifstream file("C:\\data\\CTAM\\Pescal++\\xtandem\\combiPeptData.csv");

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
	std::ifstream file("C:\\data\\CTAM\\Pescal++\\mascot\\peakAreas.csv");

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
	std::ifstream file("C:\\data\\CTAM\\Pescal++\\xtandem\\peakAreas.csv");

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
	std::ifstream file("C:\\data\\CTAM\\Pescal++\\mascot\\calculatedRTs.csv");

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
	std::ifstream file("C:\\data\\CTAM\\Pescal++\\xtandem\\calculatedRTs.csv");

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