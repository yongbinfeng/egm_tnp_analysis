#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "TH1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"

/// include pdfs
#include "RooCBExGaussShape.h"
#include "RooCMSShape.h"

#include <vector>
#include <string>
#ifdef __CINT__
#pragma link C++ class std::vector<std::string>+;
#endif

using namespace RooFit;
using namespace std;

class scaleFitter {
public:
  scaleFitter( TFile *file, std::string histname, int replica  );
  scaleFitter( TH1 *hPass, std::string histname, int replica  );
  ~scaleFitter(void) {if( _work != 0 ) delete _work; }
  void setZLineShape(TH1 *hZPass);
  void setWorkspace(std::vector<std::string>);
  void setOutputFile(TFile *fOut ) {_fOut = fOut;}
  void fits(bool mcTruth,std::string title = "");
  void useMinos(bool minos = true) {_useMinos = minos;}
  void textParForCanvas(RooFitResult *resP, TPad *p);
  
  void fixSigmaFtoSigmaP(bool fix=true) { _fixSigmaFtoSigmaP= fix;}

  void setFitRange(double xMin,double xMax) { _xFitMin = xMin; _xFitMax = xMax; }
private:
  RooWorkspace *_work;
  std::string _histname_base;
  TFile *_fOut;
  double _nTotP;
  double _bias;
  bool _useMinos;
  bool _fixSigmaFtoSigmaP;
  double _xFitMin,_xFitMax;
  int _nMCReplica;
};

scaleFitter::scaleFitter(TFile *filein, std::string histname, int replica   ) : _useMinos(false),_fixSigmaFtoSigmaP(false) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  _nMCReplica = replica;

  TString fullname = replica>0 ? TString::Format("%s_Stat%d",histname.c_str(), replica) :  TString::Format("%s",histname.c_str());
  _histname_base = fullname;  
  TH1 *hPass = (TH1*) filein->Get(fullname.Data());
  _nTotP = hPass->Integral();
  /// MC histos are done between 50-130 to do the convolution properly
  /// but when doing MC fit in 60-120, need to zero bins outside the range
  for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ )
   if(  hPass->GetXaxis()->GetBinCenter(ib) <= 60 || hPass->GetXaxis()->GetBinCenter(ib) >= 120 ) {
     hPass->SetBinContent(ib,0);
   }
  
  _work = new RooWorkspace("w") ;
  _work->factory("x[50,130]");

  RooDataHist rooPass("hPass","hPass",*_work->var("x"),hPass);
  _work->import(rooPass) ;
  _xFitMin = 60;
  _xFitMax = 120;
}

scaleFitter::scaleFitter(TH1 *hPass, std::string histname, int replica  ) : _useMinos(false),_fixSigmaFtoSigmaP(false) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  _histname_base = histname;
  _nMCReplica = replica;

  _nTotP = hPass->Integral();  
  /// MC histos are done between 50-130 to do the convolution properly
  /// but when doing MC fit in 60-120, need to zero bins outside the range
  for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ )
    if(  hPass->GetXaxis()->GetBinCenter(ib) <= 60 || hPass->GetXaxis()->GetBinCenter(ib) >= 120 ) {
      hPass->SetBinContent(ib,0);
    }
  
  _work = new RooWorkspace("w") ;
  _work->factory("x[50,130]");
  
  RooDataHist rooPass("hPass","hPass",*_work->var("x"),hPass);
  _work->import(rooPass) ;
  _xFitMin = 60;
  _xFitMax = 120;
  
}


void scaleFitter::setZLineShape(TH1 *hZPass) {
  RooDataHist rooPass("hGenZPass","hGenZPass",*_work->var("x"),hZPass);
  _work->import(rooPass) ;
}

void scaleFitter::setWorkspace(std::vector<std::string> workspace) {
  for( unsigned icom = 0 ; icom < workspace.size(); ++icom ) {
    _work->factory(workspace[icom].c_str());
  }

  _work->factory("HistPdf::sigPhysPass(x,hGenZPass)");
  _work->factory("FCONV::sigPass(x, sigPhysPass , sigResPass)");
  _work->factory(TString::Format("nSigP[%f,0.5,%f]",_nTotP*0.9,_nTotP*1.5));
  _work->factory(TString::Format("nBkgP[%f,0.5,%f]",_nTotP*0.1,_nTotP*1.5));
  _work->factory("SUM::pdfPass(nSigP*sigPass,nBkgP*bkgPass)");
  _work->Print();			         
}

void scaleFitter::fits(bool mcTruth,string title) {

  cout << " title : " << title << endl;

  
  RooAbsPdf *pdfPass = _work->pdf("pdfPass");

  if( mcTruth ) {
    _work->var("nBkgP")->setVal(0); 
    _work->var("nBkgP")->setConstant();
    if( _work->var("sosP")   ) { 
      _work->var("sosP")->setVal(0);
      _work->var("sosP")->setConstant(); 
    }
    if( _work->var("acmsP")  ) _work->var("acmsP")->setConstant();
    if( _work->var("betaP")  ) _work->var("betaP")->setConstant();
    if( _work->var("gammaP") ) _work->var("gammaP")->setConstant();
  }

  /// FC: seems to be better to change the actual range than using a fitRange in the fit itself (???)
  /// FC: I don't know why but the integral is done over the full range in the fit not on the reduced range
  _work->var("x")->setRange(_xFitMin,_xFitMax);
  _work->var("x")->setRange("fitMassRange",_xFitMin,_xFitMax);
  RooFitResult* resPass = pdfPass->fitTo(*_work->data("hPass"),Minos(_useMinos),SumW2Error(kTRUE),Save(),Range("fitMassRange"));
  //RooFitResult* resPass = pdfPass->fitTo(*_work->data("hPass"),Minos(_useMinos),SumW2Error(kTRUE),Save());

  RooPlot *pPass = _work->var("x")->frame(60,120);
  pPass->SetTitle("passing probe");
  
  _work->data("hPass") ->plotOn( pPass );
  _work->pdf("pdfPass")->plotOn( pPass, LineColor(kRed) );
  _work->pdf("pdfPass")->plotOn( pPass, Components("bkgPass"),LineColor(kBlue),LineStyle(kDashed));
  _work->data("hPass") ->plotOn( pPass );
  
  TCanvas c("c","c",1100,450);
  c.Divide(3,1);
  TPad *padText = (TPad*)c.GetPad(1);
  textParForCanvas( resPass, padText );
  c.cd(2); pPass->Draw();

  _fOut->cd();
  c.Write(TString::Format("%s_Canv_Stat%d",_histname_base.c_str(),_nMCReplica),TObject::kOverwrite);
  resPass->Write(TString::Format("%s_resP_Stat%d",_histname_base.c_str(),_nMCReplica),TObject::kOverwrite);

  
}





/////// Stupid parameter dumper /////////
void scaleFitter::textParForCanvas(RooFitResult *resP, TPad *p) {

  double eff = -1;
  double e_eff = 0;

  RooRealVar *nSigP  = _work->var("nSigP");
  RooRealVar *meanP  = _work->var("meanP");
  RooRealVar *sigmaP = _work->var("sigmaP");
  
  double nP   = nSigP->getVal();
  double e_nP = nSigP->getError();

  TPaveText *text1 = new TPaveText(0,0.8,1,1);
  text1->SetFillColor(0);
  text1->SetBorderSize(0);
  text1->SetTextAlign(12);

  text1->AddText(TString::Format("* fit status: %d",resP->status()));

//  text->AddText("* Passing parameters");
  TPaveText *text = new TPaveText(0,0,1,0.8);
  text->SetFillColor(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->AddText("    --- parmeters " );
  RooArgList listParFinalP = resP->floatParsFinal();
  for( int ip = 0; ip < listParFinalP.getSize(); ip++ ) {
    TString vName = listParFinalP[ip].GetName();
    text->AddText(TString::Format("   - %s \t= %1.3f #pm %1.3f",
				  vName.Data(),
				  _work->var(vName)->getVal(),
				  _work->var(vName)->getError() ) );
  }


  p->cd();
  text1->Draw();
  text->Draw();
}
