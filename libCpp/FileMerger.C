#include "TROOT.h"
#include "TH1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TKey.h"

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>

void FileMerger(const unsigned int bins, const char* filename, const std::vector<std::string>& names) {
	if (bins != names.size())
        throw std::runtime_error("The number of RooFit input files does not correspond to the number of bins.");
	TFile *outputfile=new TFile(filename, "RECREATE");
	outputfile->cd();
	for (unsigned int i = 0; i != names.size(); i++) {
		TFile *file = new TFile(names[i].c_str());
		if (file->IsZombie()) {
			std::string string(names[i]);
			string += std::string(" is corrupted.");
			throw std::runtime_error(string.c_str());
		}
		TNamed *obj;
		TKey *key;
		TIter next( file->GetListOfKeys());
		while ((key = (TKey *) next())) {
			obj = (TNamed*) file->Get(key->GetName()); // copy object to memory
			obj->SetName(key->GetName());
			outputfile->cd();
			obj->Write(key->GetName());
            // std::cout << "Saving object " << key->GetName() << std::endl;
		}
		file->Close();
	}
	outputfile->Close();
}

int main(int argc, char **argv) {
    
	std::vector<std::string> names;
	for (int i = 3; i != argc; i++) {
		names.push_back(std::string(argv[i]));
	}
	try {
        std::istringstream myStream(argv[1]);
        unsigned int nBins = 0;
        myStream >> nBins;
        const char* outfilename = argv[2];
		FileMerger(nBins, outfilename, names);
	}
	catch (std::runtime_error& exception) {
		std::cerr << exception.what() << "\n";
		return -1;
	}
	return 0;
}
