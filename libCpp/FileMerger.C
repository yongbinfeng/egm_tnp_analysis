#include "TROOT.h"
#include "TH1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TKey.h"

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

void FileMerger(const char* NumberOfBins, const char* filename, std::vector<std::string> names) {
	istringstream myStream(NumberOfBins);
	unsigned int bins = 0;
	myStream >> bins;
	if (bins!=names.size()) throw std::runtime_error("The number of RooFit output files does not correspond to the number of bins.");
	TFile *outputfile=new TFile(filename,"RECREATE");
	outputfile->cd();
	for (unsigned int i=0; i!=names.size(); i++) {
		TFile *file=new TFile(names[i].c_str());
		if (file->IsZombie()) {
			std::string string(names[i]);
			string+=std::string(" is corrupted.");
			throw std::runtime_error(string.c_str());
		}
		TObject *obj;
		TKey *key;
		TIter next( file->GetListOfKeys());
		while ((key = (TKey *) next())) {
			obj = file->Get(key->GetName()); // copy object to memory
			outputfile->cd();
			obj->Write();
		}
		file->Close();
	}
	outputfile->Close();
}

int main(int argc, char **argv) {
	std::vector<std::string> names;
	for (int i=3; i!=argc; i++) {
		names.push_back(std::string(argv[i]));
	}
	try {
		FileMerger(argv[1],argv[2],names);
	}
	catch (std::runtime_error& exception) {
		std::cerr<<exception.what()<<"\n";
		return -1;
	}
	return 0;
}
