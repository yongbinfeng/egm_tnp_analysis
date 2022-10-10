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

class tnpFitter {
public:
  tnpFitter( TFile *file, std::string histname, int massbins, float massmin, float massmax  );
  tnpFitter( TH1 *hPass, TH1 *hFail, std::string histname, int massbins, float massmin, float massmax  );
  ~tnpFitter(void) {if( _work != 0 ) delete _work; }
  void setZLineShapes(TH1 *hZPass, TH1 *hZFail );
  void setWorkspace(std::vector<std::string>, bool, bool, bool);
  //python3void setOutputFile( TFile *fOut ) {_fOut = fOut;}
  //void setOutputFile(TString fname ) {_fOut = new TFile(fname, "UPDATE");}
  void setOutputFile(TString fname ) {_fOut = new TFile(fname, "recreate"); } 
  int fits(std::string title = "");
  void useMinos(bool minos = true) {_useMinos = minos;}
  void textParForCanvas(RooFitResult *resP, RooFitResult *resF, TPad *p);
  void fixSigmaFtoSigmaP(bool fix=true) { _fixSigmaFtoSigmaP= fix;}
  void setFitRange(double xMin,double xMax) { _xFitMin = xMin; _xFitMax = xMax; }
  void setPassStrategy(int strategy) {_strategyPassFit = std::clamp(strategy, 0, 2); } 
  void setFailStrategy(int strategy) {_strategyFailFit = std::clamp(strategy, 0, 2); } 
  void setPrintLevel(int level) {_printLevel = std::clamp(level, -1, 9); } 
  void setMaxSignalFractionFail(double max) { _maxSignalFractionFail = max; }
    
private:
  RooWorkspace *_work;
  std::string _histname_base;
  TFile *_fOut;
  double _nTotP, _nTotF;
  bool _useMinos;
  bool _fixSigmaFtoSigmaP;
  double _xFitMin,_xFitMax;
  int _strategyPassFit = 1;
  int _strategyFailFit = 1;
  int _printLevel = 3;
  double _maxSignalFractionFail = -1; // not used by default if negative
};

tnpFitter::tnpFitter(TFile *filein, std::string histname, int massbins, float massmin, float massmax ) : _useMinos(false),_fixSigmaFtoSigmaP(false) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  _histname_base = histname;  

  TH1 *hPass = (TH1*) filein->Get(TString::Format("%s_Pass",histname.c_str()).Data());
  TH1 *hFail = (TH1*) filein->Get(TString::Format("%s_Fail",histname.c_str()).Data());
  _nTotP = hPass->Integral();
  _nTotF = hFail->Integral();
  /// MC histos are done between 50-130 to do the convolution properly
  /// but when doing MC fit in 60-120, need to zero bins outside the range
  for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ )
   if(  hPass->GetXaxis()->GetBinCenter(ib) <= massmin || hPass->GetXaxis()->GetBinCenter(ib) >= massmax ) {
     hPass->SetBinContent(ib,0);
     hFail->SetBinContent(ib,0);
   }
  
  _work = new RooWorkspace("w") ;
  //_work->factory("x[50,130]");
  _work->factory(TString::Format("x[%f,%f]",massmin, massmax));

  RooDataHist rooPass("hPass","hPass",*_work->var("x"),hPass);
  RooDataHist rooFail("hFail","hFail",*_work->var("x"),hFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;
  _xFitMin = massmin;
  _xFitMax = massmax;
}

tnpFitter::tnpFitter(TH1 *hPass, TH1 *hFail, std::string histname, int massbins, float massmin, float massmax ) : _useMinos(false),_fixSigmaFtoSigmaP(false) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  _histname_base = histname;
  
  _nTotP = hPass->Integral();
  _nTotF = hFail->Integral();
  /// MC histos are done between 50-130 to do the convolution properly
  /// but when doing MC fit in 60-120, need to zero bins outside the range
  for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ )
    if(  hPass->GetXaxis()->GetBinCenter(ib) <= massmin || hPass->GetXaxis()->GetBinCenter(ib) >= massmax ) {
      hPass->SetBinContent(ib,0);
      hFail->SetBinContent(ib,0);
    }
 
  _work = new RooWorkspace("w") ;
  //_work->factory("x[50,130]");
  _work->factory(TString::Format("x[%f,%f]",massmin, massmax));


  RooDataHist rooPass("hPass","hPass",*_work->var("x"),hPass);
  RooDataHist rooFail("hFail","hFail",*_work->var("x"),hFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;
  _xFitMin = massmin;
  _xFitMax = massmax;
  
}


void tnpFitter::setZLineShapes(TH1 *hZPass, TH1 *hZFail ) {
  RooDataHist rooPass("hGenZPass","hGenZPass",*_work->var("x"),hZPass);
  RooDataHist rooFail("hGenZFail","hGenZFail",*_work->var("x"),hZFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;  
}

void tnpFitter::setWorkspace(std::vector<std::string> workspace, bool isMCfit = false, bool analyticPhysicsShape = false, bool modelFSR = false) {
  for( unsigned icom = 0 ; icom < workspace.size(); ++icom ) {
    _work->factory(workspace[icom].c_str());
  }

  if (not analyticPhysicsShape) {
      _work->factory("HistPdf::sigPhysPass(x,hGenZPass)");
      _work->factory("HistPdf::sigPhysFail(x,hGenZFail)");
  }
  
  _work->factory("FCONV::sigPass(x, sigPhysPass , sigResPass)");
  _work->factory(TString::Format("nSigP[%f,0.5,%f]",_nTotP*0.9,_nTotP*1.5));
  _work->factory(TString::Format("nBkgP[%f,0.5,%f]",_nTotP*0.1,_nTotP*1.5));
  _work->factory("SUM::pdfPass(nSigP*sigPass,nBkgP*bkgPass)");

  if (modelFSR) {

      _work->factory("FCONV::sigMainFail(x, sigPhysFail , sigResFail)");
      _work->factory("SUM::sigFail(fracMainF[0.95,0.8,1.0]*sigMainFail, sigFsrFail)");
      if (isMCfit) {
          _work->factory(TString::Format("nSigF[%f,%f,%f]",_nTotF*0.9,_nTotF*0.85,_nTotF*1.5));
          _work->factory(TString::Format("nBkgF[%f,0.5,%f]",_nTotF*0.1,_nTotF*0.15));
      } else{ 
          _work->factory(TString::Format("nSigF[%f,0.5,%f]",_nTotF*0.9,_nTotF*1.5));
          _work->factory(TString::Format("nBkgF[%f,0.5,%f]",_nTotF*0.1,_nTotF*1.5));
      }
      _work->factory("SUM::pdfFail(nSigF*sigFail,nBkgF*bkgFail)");
      
  } else {

      _work->factory("FCONV::sigFail(x, sigPhysFail , sigResFail)");
      if (isMCfit) {
          _work->factory(TString::Format("nSigF[%f,%f,%f]",_nTotF*0.9,_nTotF*0.85,_nTotF*1.5));
          _work->factory(TString::Format("nBkgF[%f,0.5,%f]",_nTotF*0.1,_nTotF*0.15));
      } else{ 
          if (_work->var("maxFracSigF") != nullptr) {
              // std::cout << "Test setting signal fraction" << std::endl;
              double maxFracSigF = _work->var("maxFracSigF")->getVal(); 
              double minFracBkgF = 1.0 - maxFracSigF;
              double halfFracSigF = maxFracSigF/2.0;
              _work->factory(TString::Format("nSigF[%f,0.5,%f]", _nTotF*halfFracSigF, _nTotF*maxFracSigF));
              _work->factory(TString::Format("nBkgF[%f,%f,%f]", _nTotF*(minFracBkgF+halfFracSigF), _nTotF*minFracBkgF,_nTotF*1.5));          
          } else {
              _work->factory(TString::Format("nSigF[%f,0.5,%f]",_nTotF*0.9,_nTotF*1.5));
              _work->factory(TString::Format("nBkgF[%f,0.5,%f]",_nTotF*0.1,_nTotF*1.5));          
          }
      }
      _work->factory("SUM::pdfFail(nSigF*sigFail,nBkgF*bkgFail)");
      
  }

  //_work->Print(); // FIXME: might want to comment this one to avoid unnecessary output text

}

int tnpFitter::fits(string title) {

  cout << " this is the title : " << title << endl;

  
  RooAbsPdf *pdfPass = _work->pdf("pdfPass");
  RooAbsPdf *pdfFail = _work->pdf("pdfFail");

  /// FC: seems to be better to change the actual range than using a fitRange in the fit itself (???)
  /// FC: I don't know why but the integral is done over the full range in the fit not on the reduced range
  _work->var("x")->setRange(_xFitMin,_xFitMax);
  _work->var("x")->setRange("fitMassRange",_xFitMin,_xFitMax);
  RooFitResult* resPass = pdfPass->fitTo(*_work->data("hPass"),Minos(_useMinos),SumW2Error(kTRUE),Save(),Range("fitMassRange"),Strategy(_strategyPassFit), PrintLevel(_printLevel));
  //RooFitResult* resPass = pdfPass->fitTo(*_work->data("hPass"),Minos(_useMinos),SumW2Error(kTRUE),Save());
  if( _fixSigmaFtoSigmaP ) {
    _work->var("sigmaF")->setVal( _work->var("sigmaP")->getVal() );
    _work->var("sigmaF")->setConstant();
  }

  // _work->var("sigmaF")->setVal(_work->var("sigmaP")->getVal());
  // _work->var("sigmaF")->setRange(0.8* _work->var("sigmaP")->getVal(), 3.0* _work->var("sigmaP")->getVal());
  RooFitResult* resFail = pdfFail->fitTo(*_work->data("hFail"),Minos(_useMinos),SumW2Error(kTRUE),Save(),Range("fitMassRange"),Strategy(_strategyFailFit), PrintLevel(_printLevel));

  RooPlot *pPass = _work->var("x")->frame(_xFitMin,_xFitMax); // always plot 50 - 130
  RooPlot *pFail = _work->var("x")->frame(_xFitMin,_xFitMax);
  pPass->SetTitle("passing probe");
  pFail->SetTitle("failing probe");

  _work->data("hPass") ->plotOn( pPass );
  _work->pdf("pdfPass")->plotOn( pPass, LineColor(kRed) );
  _work->pdf("pdfPass")->plotOn( pPass, Components("bkgPass"),LineColor(kBlue),LineStyle(kDashed));
  _work->data("hPass") ->plotOn( pPass );
  
  _work->data("hFail") ->plotOn( pFail );
  _work->pdf("pdfFail")->plotOn( pFail, LineColor(kRed) );
  _work->pdf("pdfFail")->plotOn( pFail, Components("bkgFail"),LineColor(kBlue),LineStyle(kDashed));
  _work->data("hFail") ->plotOn( pFail );

  TCanvas * c = new TCanvas(TString::Format("%s_Canv",_histname_base.c_str()) ,TString::Format("%s_Canv",_histname_base.c_str()) ,1100,450);
  c->Divide(3,1);
  TPad *padText = (TPad*)c->GetPad(1);
  textParForCanvas( resPass,resFail, padText );
  c->cd(2); pPass->Draw();
  c->cd(3); pFail->Draw();

  _fOut->cd();
  c->Write(TString::Format("%s_Canv",_histname_base.c_str()),TObject::kOverwrite);
  resPass->Write(TString::Format("%s_resP",_histname_base.c_str()),TObject::kOverwrite);
  resFail->Write(TString::Format("%s_resF",_histname_base.c_str()),TObject::kOverwrite);
  _fOut->Close();
  
  return 1;

}





/////// Stupid parameter dumper /////////
void tnpFitter::textParForCanvas(RooFitResult *resP, RooFitResult *resF,TPad *p) {

  double eff = -1;
  double e_eff = 0;

  RooRealVar *nSigP = _work->var("nSigP");
  RooRealVar *nSigF = _work->var("nSigF");
  
  double nP   = nSigP->getVal();
  double e_nP = nSigP->getError();
  double nF   = nSigF->getVal();
  double e_nF = nSigF->getError();
  double nTot = nP+nF;
  eff = nP / (nP+nF);
  e_eff = 1./(nTot*nTot) * sqrt( nP*nP* e_nF*e_nF + nF*nF * e_nP*e_nP );  // this is linear error propagation assuming uncorrelated nP and nF, but might not be correct when efficiency is close to 1

  TPaveText *text1 = new TPaveText(0,0.84,1,1);
  text1->SetFillColor(0);
  text1->SetBorderSize(0);
  text1->SetTextAlign(12);

  text1->AddText(TString::Format("* fit status pass: %d, fail : %d",resP->status(),resF->status()));
  text1->AddText(TString::Format("* eff = %1.4f #pm %1.4f",eff,e_eff));

  //  text->SetTextSize(0.06);

//  text->AddText("* Passing parameters");
  TPaveText *text = new TPaveText(0,0,1,0.84);
  text->SetFillColor(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->AddText("    --- parameters " );
  RooArgList listParFinalP = resP->floatParsFinal();
  for( int ip = 0; ip < listParFinalP.getSize(); ip++ ) {
    TString vName = listParFinalP[ip].GetName();
    text->AddText(TString::Format("   - %s \t= %1.3f #pm %1.3f",
				  vName.Data(),
				  _work->var(vName)->getVal(),
				  _work->var(vName)->getError() ) );
  }

//  text->AddText("* Failing parameters");
  RooArgList listParFinalF = resF->floatParsFinal();
  for( int ip = 0; ip < listParFinalF.getSize(); ip++ ) {
    TString vName = listParFinalF[ip].GetName();
    text->AddText(TString::Format("   - %s \t= %1.3f #pm %1.3f",
				  vName.Data(),
				  _work->var(vName)->getVal(),
				  _work->var(vName)->getError() ) );
  }

  p->cd();
  text1->Draw();
  text->Draw();
}
