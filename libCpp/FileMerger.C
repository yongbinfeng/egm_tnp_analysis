#include "TROOT.h"
#include "TH1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TKey.h"

#include <vector>
#include <string>
#include <unordered_map>

void FileMerger(const char* filename,std::vector<std::string> names) {
	TFile *outputfile=new TFile(filename,"RECREATE");
	outputfile->cd();
	for (unsigned int i=0; i!=names.size(); i++) {
		TFile *file=new TFile(names[i].c_str());
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
	for (int i=2; i!=argc; i++) {
		names.push_back(std::string(argv[i]));
	}
	FileMerger(argv[1],names);
	return 0;
}
