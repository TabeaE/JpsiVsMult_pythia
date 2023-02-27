#include <iostream>
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TMath.h"
#include "TPaveText.h"
#include "defineEstimator.C"
#include "defineVariation.C"
#include "defineAliases.C"
using namespace std;
/**
 * 
 * The macro that plots J/psi production vs multiplicity
 * In addition:
 *  - Ratio to linear increase
 *  - pT, y distribution
 *  - fractions vs. multiplicity, pT and y
 * 
 * It uses functions defined in 
 *  - defineEstimator.C
 *  - defineAliases.C
 *  - defineVariation.C
 * 
 * 
 * Arguments:
 *
 *  - plotWhat: int, what should be plotted (defined in defineVariation.C). Options:
 *              1: plot inclusive, prompt and non-prompt J/psi separately
 *              2: plot inclusive J/psi in pT bins (< 4, 4 - 8, > 8 GeV/c)
 *              3: plot non-prompt J/psi in pT bins
 *              4: plot J/psi in different decay channels (dilepton and others)
 *              5: plot different charmonium states (J/psi, Psi(2s))
 *  - selector: int, which multiplicity estimator to use for the event selection
 *              (defined in defineEstimator.C). Options:
 *              e.g. 1 is nch in |eta| < 1
 *                  13 is nch in V0C in towards region
 *  - estimator: int, which multiplicity estimator to use for the x axis
 *               same meaning of the number as for "selector"
 *               so, it is possible to select events based on the V0M multiplicity, but plot as a function
 *               of nch in |eta| < 1
 *  - setSettings: string, name of input root file (in format onium_setSettings.root)
 *  - setEventClass: int, select event class
 *                   0: Not-diffractive events
 *                   1: All inelastic events
 *                   2: Diffractive events
 *                   3: Not single-diffractive events
 *  - setEventTrigger: int, select event trigger
 *                     0: all events
 *                     1: V0AND
 *                     2: INEL>0
 *  - signalRapidity: int, cut on J/psi rapidity
 *        1: |y| < 0.9
 *        2:  2.5 < |y| <  4.0
 *        3:  2.5 <  y  <  4.0
 *        4: -2.5 >  y  > -4.0
 *  - maxY: double, maximum value of y axis
 *  - folder: string, path to the input tree file
 * 
 * 
 */

void plotRelOrigin(
  int     plotWhat        = 1,
  int     selector        = 4,
  int     estimator       = 1,
  TString setSettings     = "full",
  int     setEventClass   = 0,
  int     setEventTrigger = 1,
  int     signalRapidity  = 1,
  double  maxY            = 22,
  TString folder          = "."
);

struct Variation;

void crash(TString message);
void drawInCanvas(TCanvas* c, TObject* o, TString option="");
TString normalized(TString str);
TCanvas* makeCanvas(TString name, bool logx=false, bool logy=false);
TLegend* makeLegend(TString header, int xPos=0, int yPos=0);
void plot(TString commonCut, TString commonName, TString header, TString inclusiveName,
          TString variationsName, vector<Variation> variations, double maxY, double maxX);
void save(TCanvas* c, TString folder, TString filename);
void styleFunction(TF1* h, Color_t color, int line, TString  title, TString  xtitle, TString ytitle,
                   double xMin=0, double xMax=0, double yMin=0, double yMax=0);
void stylePlot(TGraphErrors* h, Color_t color, int style, int line, TString  title, TString  xtitle,
               TString ytitle, double xMin=0, double xMax=0, double yMin=0, double yMax=0);


TTree* oniumTree;

double nEventsMB;
double meanMultMB;
const int maxPt   = 16;
const int nPtBins = 32;

const int maxRapidity   = 8;
const int nRapidityBins = 80;


TString nJpsi;
TString estimatorStringForAxisTitle;
TString selectorStringForAxisTitle;


TString settings;
TString estimatorString = "";
TString selectorString  = "";
int eventClass;
int eventTrigger;
bool noMPI;
bool noCR;

// Multiplicity bins are first multiplied by the factor and then shifted by -0.5
// pPb 5TeV, global tracks
const unsigned int nMultiplicityBins = 10;
double multiplicityBins[nMultiplicityBins] = {
         1., 20., 27., 34.,  44.,
        54., 64., 79., 94., 119.
      };

// pp 13TeV, SPD tracklets
// const unsigned int nMultiplicityBins = 24;
// double multiplicityBins[nMultiplicityBins] = {
//          0.,  4.,   8.,  12., 16.,
//         20., 24.,  28.,  32., 36.,
//         40., 44.,  48.,  52., 56.,
//         60., 64.,  68.,  72., 76.,
//         84., 92., 100., 116.
//       };
      
double meanMultMBPerBin[nMultiplicityBins];
double nEventsPerBin[nMultiplicityBins];
double errPerBin[nMultiplicityBins];


void plotRelOrigin(int plotWhat, int selector, int estimator, TString setSettings, int setEventClass,
                   bool setEventTrigger, int signalRapidity, double maxY, TString folder) {

  gEnv->SetValue("Hist.Precision.1D", "double");

  // Initialize values
  settings     = setSettings;
  eventClass   = setEventClass;
  eventTrigger = setEventTrigger;

  noMPI = strstr(settings.Data(), "noMPI");
  noCR  = strstr(settings.Data(), "noCR");

  nJpsi = "#it{N}_{J/#it{#psi}}";
  nJpsi = Form("d%s/d#it{y}", nJpsi.Data());
  nJpsi = normalized(nJpsi);
  

  // Set plot options
  gStyle->SetOptStat(0);
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleX(0.58);
  

  // Define multiplicity estimator
  TString regionString_estimator;
  TString regionString_selector;
  Double_t factor;
  
  defineEstimator(estimator, estimatorString, estimatorStringForAxisTitle, regionString_estimator, factor);
  defineEstimator(selector,  selectorString,  selectorStringForAxisTitle,  regionString_selector,  factor);
  
  cout << "estimator: " << estimatorString << endl;
  cout << "selector: "  << selectorString  << endl;
  
  if(noCR)  factor *= 3.;
  if(noMPI) factor /= 2.;
  
  if(estimator!=6) estimatorStringForAxisTitle = Form("d%s/d#it{#eta}", estimatorStringForAxisTitle.Data());
  estimatorStringForAxisTitle = normalized(estimatorStringForAxisTitle); 
  estimatorStringForAxisTitle = Form("%s(%s, selection: %s)", estimatorStringForAxisTitle.Data(),
                                     regionString_estimator.Data(), regionString_selector.Data());
  
  for(unsigned int i=0; i<nMultiplicityBins; ++i) {
    multiplicityBins[i] *= factor;
    multiplicityBins[i] -= 0.5;
  }


  // Read input file
  TString fn = "onium";
  TFile* f = TFile::Open(Form("%s/%s_%s.root", folder.Data(), fn.Data(), settings.Data()));
  cout << "opening " <<  Form("%s/%s_%s.root", folder.Data(), fn.Data(), settings.Data()) << endl;

  // Get tree
  TTree* eventTree = (TTree*) f->Get("eventTree");
  oniumTree = (TTree*) f->Get("oniumTree");
  
  
  // Define aliases
  defineAliases(eventTree, oniumTree);
  eventTree->SetAlias("mul", estimatorString.Data());
  eventTree->SetAlias("sel", selectorString.Data());
  oniumTree->SetAlias("sel", selectorString.Data());
  

  // Apply event cuts (select event trigger and class)
  TString eventCut;
  switch(eventClass) {
    case 0:  // Not-diffractive inelastic events
      eventCut =  (eventTrigger==0) ? "inelastic && !diffractive"
                : (eventTrigger==1) ? "inelastic && !diffractive && V0AND"
                : (eventTrigger==2) ? "inelastic && !diffractive && inel0";
      break;
    case 1:  // All inelastic events
      eventCut =  (eventTrigger==0) ? "inelastic"
                : (eventTrigger==1) ? "inelastic && V0AND"
                : (eventTrigger==2) ? "inelastic && inel0";
      break;
    case 2:  // Diffractive events
      eventCut =  (eventTrigger==0) ? "diffractive"
                : (eventTrigger==1) ? "V0AND && diffractive"
                : (eventTrigger==2) ? "inel0 && diffractive";
      break;
    case 3:  // Not single-diffractive events
      eventCut =  (eventTrigger==0) ? "NSD"
                : (eventTrigger==1) ? "V0AND && NSD"
                : (eventTrigger==2) ? "inel0 && NSD";
      break;
  }
  
  eventTree->SetAlias("eventCut", eventCut.Data());
  oniumTree->SetAlias("eventCut", eventCut.Data());
  
  
  // IMPORTANT
  // Draw the 2D multiplicity distribution for MB events, i.e. the multiplicity used for selecting events
  // vs. the one used as x axis
  TH2D* hMult = new TH2D("hMult", "hMult", 1000, -0.5, 999.5, 1000, -0.5, 999.5);
  hMult->Sumw2();
  hMult->SetDirectory(gDirectory);
  eventTree->Draw("mul:sel>>hMult", "eventCut", "goff");
  hMult->GetXaxis()->SetRange(0, hMult->GetNbinsX());

  meanMultMB = hMult->GetMean(2); 
  nEventsMB  = hMult->Integral();
  cout << "nEventsMB " << nEventsMB << endl << "mean mult: " << meanMultMB << endl;
  if(!nEventsMB)  crash("nEventsMB == 0!!!");
  if(!meanMultMB) crash("meanMultMB == 0!!!");

  for(unsigned int i=0; i<nMultiplicityBins-1; ++i) {
      hMult->GetXaxis()->SetRangeUser(multiplicityBins[i], multiplicityBins[i+1]);
      meanMultMBPerBin[i] = hMult->GetMean(2);
      int binMin = hMult->GetXaxis()->FindBin(multiplicityBins[i] );
      int binMax = hMult->GetXaxis()->FindBin(multiplicityBins[i+1])-1;
      nEventsPerBin[i] = hMult->IntegralAndError(binMin, binMax, 0, 1000, errPerBin[i]);
  }


  // Define Jpsi rapidity range
  TString commonCut, commonName, variationsName, header, inclusiveName;
  vector<Variation> variations;
  defineVariation(plotWhat, &commonName, &variationsName, &commonCut, &header, &inclusiveName, &variations);
  double maxX = 7;

  if(signalRapidity==1) {
    commonCut  += "&&midRapidity";
    header     += " (|#it{y}|<0.9)";
    commonName += "Mid";
  }
  else if(signalRapidity==2) {
    commonCut  += "&&muonRapidity";
    header     += " (2.5<|#it{y}|<4.0)";
    commonName += "Muons";
  }
  else if(signalRapidity==3) {
    commonCut  += "&&muonRapidityPos";
    header     += " (2.5<#it{y}<4.0)";
    commonName += "MuonsFwd";
  }
  else if(signalRapidity==4) {
    commonCut  += "&&muonRapidityNeg";
    header     += " (-2.5>#it{y}>-4.0)";
    commonName += "MuonsBkwd";
  }


  // Plot
  plot(commonCut, commonName, header, inclusiveName, variationsName, variations, maxY, maxX);

}


void plot(TString commonCut, TString commonName, TString header, TString inclusiveName,
          TString variationsName, vector<Variation> variations, double maxY, double maxX) {

  cout << commonName << " vs. " << variationsName << endl;
  
  TString title;
  TString eventClassString = "";
  switch(eventClass) {
    case 0:  // Not-diffractive events
      eventClassString = "nondiffractive";
      title           += "nondiffractive events";
      break;
    case 1:  // All inelastic events
      eventClassString = "inelastic";
      title           += "inelastic events";
      break;
    case 2:  // Diffractive events
      eventClassString = "diffractive";
      title           += "diffractive events";
      break;
    case 3:  // Not single-diffractive events
      eventClassString = "NSD";
      title           += "NSD events";
      break;
  }
  
  TString eventTriggerString = "";
  switch(eventClass) {
    case 1:  // V0AND
      eventClassString = "V0AND";
      title           += "V0AND";
      break;
    case 2:  // INEL>0
      eventClassString = "INEL0";
      title           += "INEL>0";
      break;
  }
  
  if(noCR)  title += ", CR off";
  if(noMPI) title += ", MPI off";
  
  TString folder1    = Form("%s%s%s", settings.Data(), eventClassString.Data(), eventTriggerString.Data());
  TString identifier = Form("%s_%s_%s", folder1.Data(), selectorString.Data(), estimatorString.Data());
  
  gSystem->mkdir(Form("output/%s", folder1.Data()));
  gSystem->mkdir(Form("output/%s/%s_%s", folder1.Data(),selectorString.Data(), estimatorString.Data()));
  TString subfolder = Form("%sVS%s", commonName.Data(), variationsName.Data());
  TString folder    = Form("output/%s/%s_%s/%s", folder1.Data(), selectorString.Data(),
                           estimatorString.Data(), subfolder.Data());
  gSystem->mkdir(folder.Data()); 
  
  TString fn_add = Form("%s_%s", identifier.Data(), subfolder.Data());
  
  
  TCanvas* cYieldVsMult    = makeCanvas("yieldVsMult");
  TCanvas* cRatioToLin     = makeCanvas("ratioToLin");
  
  TCanvas* cYieldVsPt      = makeCanvas("yieldVsPt", false, true);
  TCanvas* cYieldVsY       = makeCanvas("yieldVsY");
  
  TCanvas* cFractionVsMult = makeCanvas("fractionVsMult");
  TCanvas* cFractionVsPt   = makeCanvas("fractionVsPt");
  TCanvas* cFractionVsY    = makeCanvas("fractionVsY");
  
  
  TLegend* lmy = makeLegend(header, 0, 0);
  TLegend* lmf = makeLegend(header, 0, 0);
  TLegend* lpy = makeLegend(header, 1, 0);
  TLegend* lpf = makeLegend(header, 1, 0);
  TLegend* lyy = makeLegend(header, 2, 1);
  TLegend* lyf = makeLegend(header, 2, 0);
  
  
  TH1D* hJpsiVsMultTotal    = 0x0;
  TH1D* hJpsiVsRawMultTotal = 0x0;
  TH1D* hJpsiVsPtTotal      = 0x0;
  TH1D* hJpsiVsYTotal       = 0x0;
  
  
  TF1* linear = new TF1("linear", "x" , 0., maxX);
  styleFunction(linear, kBlack, 2, title, estimatorStringForAxisTitle, nJpsi, 0., maxX, 0., maxY);
  drawInCanvas(cYieldVsMult, linear, "L");
  
  TF1* flat = new TF1("flat", "1", 0., maxX);
  styleFunction(flat, kBlack, 2, title, estimatorStringForAxisTitle, "ratio to lin", 0., maxX, 0., 2.);
  drawInCanvas(cRatioToLin, flat, "L");

  TF1* null = new TF1("null", "0", 0., maxX);
  styleFunction(null, kWhite, 2, title, estimatorStringForAxisTitle, "Fraction of total yield",0.,maxX,0.,1.);
  drawInCanvas(cFractionVsMult, null, "L");
  
  TF1* onePt  = new TF1("onePt",  "0", 0., maxPt);
  TF1* nullPt = new TF1("nullPt", "0", 0., maxPt);
  styleFunction(nullPt, kWhite, 2, title, "#it{p}_{T} (GeV/#it{c})","Fraction of total yield",0.,maxPt,0.,1.);
  drawInCanvas(cFractionVsPt, nullPt, "L");
  
  TF1* one = new TF1("one", "1", 0., maxPt);
  styleFunction(onePt, kWhite, 2, title,"#it{p}_{T} (GeV/#it{c})", "yield( a.u.)", 0., maxPt, 0.01, 2. );
  drawInCanvas(cYieldVsPt, onePt, "L");
  
  TF1* nullY = new TF1("nullY", "0", -1.*maxRapidity, maxRapidity);
  TF1* oneY  = new TF1("oneY",  "1", -1.*maxRapidity, maxRapidity);
  styleFunction(nullY, kWhite, 2, title, "y", "Fraction of total yield", -1.*maxRapidity, maxRapidity,0.,1.3);
  drawInCanvas(cFractionVsY, nullY, "L");
  
  styleFunction(oneY, kWhite, 2, title, "y", "yield (a.u.)", -1.*maxRapidity, maxRapidity, 0., 1.1);
  drawInCanvas(cYieldVsY, oneY, "L");

  TFile out(Form("%s/output.root",folder.Data()), "RECREATE");
  
  for(int iCut=-1; iCut<(int)variations.size(); ++iCut) {
      
    bool cutBin     = iCut > -1;
    
    TString cut     = cutBin ? variations[iCut].cut   : "1" ;
    TString name    = cutBin ? variations[iCut].name  : inclusiveName;
    Color_t color   = cutBin ? variations[iCut].color : 1;
    int style       = cutBin ? variations[iCut].style : 20;
    int lineStyle   = 1;
    TString optionY = "PLsame";
    TString optionF = "PLsame";
    
    cout << "  " << name << endl;

    
    TGraphErrors* gFractionVsMult = new TGraphErrors();
    TGraphErrors* gFractionVsPt   = new TGraphErrors();
    TGraphErrors* gFractionVsY    = new TGraphErrors();
    
    TGraphErrors* gYieldVsMult    = new TGraphErrors();
    TGraphErrors* gRatioToLin     = new TGraphErrors();
    TGraphErrors* gRatioToQuad    = new TGraphErrors();
    TGraphErrors* gYieldVsPt      = new TGraphErrors();
    TGraphErrors* gYieldVsY       = new TGraphErrors();
    
    gDirectory->Delete("hJpsi");
    TH3D* hJpsi = new TH3D("hJpsi", "hJpsi", nRapidityBins, -1.*maxRapidity, maxRapidity,
                            200,  0.0, 100.0,
                           1000, -0.5, 999.5);
    hJpsi->SetDirectory(gDirectory);
    
    // IMPORTANT
    // Draw the multiplicity vs pt vs rapidity distribution of selected quarkonia into a 3D histogram
    oniumTree->Draw("sel:onium.pt:onium.y>>hJpsi", Form("eventCut && (%s) && (%s)",commonCut.Data(), cut.Data()), "goff");
  
  
  // Multiplicity
    TH1D* hJpsiVsMult = hJpsi->ProjectionZ();
    double meanYield  = hJpsiVsMult->Integral() / nEventsMB;
    hJpsiVsMult = (TH1D*) hJpsiVsMult->Rebin(nMultiplicityBins-1, "jpsiNewBins", multiplicityBins);
    
    if(!cutBin) {
      hJpsiVsMultTotal = (TH1D*) hJpsiVsMult->Clone();
      hJpsiVsMultTotal->SetDirectory(0);
    }
    
    if(!meanYield) {
      cout << "meanYield == 0 !!!" << endl;
      continue;
    }
    
    for(unsigned int iMultBin=0; iMultBin<nMultiplicityBins-1; ++iMultBin) {
      
      if(!nEventsPerBin[iMultBin]) {
        cout << (Form("nEventsPerBin[%d] == 0!!!", iMultBin)) << endl;
        nEventsPerBin[iMultBin] = 1;
      }

      // calculate relative multiplicity
      double mult = meanMultMBPerBin[iMultBin] / meanMultMB;
      
      // error calculation.
      // assumption: mean yield has no statistical uncertainty
      double nJpsi = hJpsiVsMult->GetBinContent(iMultBin+1);
      double eJpsi = hJpsiVsMult->GetBinError(iMultBin+1);
      double nEv   = nEventsPerBin[iMultBin];
      double eEv   = errPerBin[iMultBin];
      
      double yield      = nJpsi / nEv ;
      double yieldError = TMath::Sqrt(TMath::Power(eJpsi/nEv,2) + TMath::Power(eEv*nJpsi/nEv/nEv,2));
      
      
      gYieldVsMult->SetPoint(      iMultBin, mult, yield/meanYield      );
      gYieldVsMult->SetPointError( iMultBin, 0.,   yieldError/meanYield );
      
      if(mult){
        gRatioToLin->SetPoint(      iMultBin, mult, yield/meanYield/mult      );
        gRatioToLin->SetPointError( iMultBin, 0.,   yieldError/meanYield/mult );
        
        gRatioToQuad->SetPoint(      iMultBin, mult, yield/meanYield/mult/mult      );
        gRatioToQuad->SetPointError( iMultBin, 0.,   yieldError/meanYield/mult/mult );
      }
      
      if(cutBin){
        double total = hJpsiVsMultTotal->GetBinContent(iMultBin+1);
        double fraction      = total ? (hJpsiVsMult->GetBinContent(iMultBin+1) / total) : 0;
        double fractionError = total ? (hJpsiVsMult->GetBinError(iMultBin+1)   / total) : 1;
        
        gFractionVsMult->SetPoint(iMultBin,      mult, fraction);
        gFractionVsMult->SetPointError(iMultBin, 0.,   fractionError);
      }
    }  // end loop over mutliplicity bins
    
    if(cutBin) {
      stylePlot(gFractionVsMult, color, style, lineStyle, title, estimatorStringForAxisTitle,
                "Fraction of total yield", 0., maxX, 0., 1.);
      drawInCanvas(cFractionVsMult, gFractionVsMult, optionF);
      lmf->AddEntry(gFractionVsMult, name.Data(), "LP");
    }
    
    stylePlot(gYieldVsMult, color, style, lineStyle, title, estimatorStringForAxisTitle,
              nJpsi, 0., maxX, 0., maxY);
    drawInCanvas(cYieldVsMult, gYieldVsMult, optionY);
    
    stylePlot(gRatioToLin, color, style, lineStyle, title, estimatorStringForAxisTitle,
              "ratio to lin", 0., maxX, 0., 2.);
    drawInCanvas(cRatioToLin, gRatioToLin, optionY);
    lmy->AddEntry(gYieldVsMult, name.Data(), "LP");
    
    
    // Transverse momentum
    TH1D* hJpsiVsPt = hJpsi->ProjectionY();
    
    if(!cutBin) {
      hJpsiVsPtTotal = (TH1D*) hJpsiVsPt->Clone();
      hJpsiVsPtTotal->SetDirectory(0);
    }
    int iPoint = 0;
    for(unsigned int iBin=0; iBin<nPtBins; ++iBin) {
      double maxYieldVsPt = hJpsiVsPt->GetMaximum();
      double yield        = maxYieldVsPt ? hJpsiVsPt->GetBinContent(iBin+1) / maxYieldVsPt : 0.;
      double yieldError   = maxYieldVsPt ? hJpsiVsPt->GetBinError(iBin+1)   / maxYieldVsPt : 0.;
      double pt           = hJpsiVsPt->GetBinCenter(iBin+1);
      
      if(yield) {
        gYieldVsPt->SetPoint(      iPoint, pt, yield);
        gYieldVsPt->SetPointError( iPoint, 0., yieldError);
      
        if(cutBin) {
          double total = hJpsiVsPtTotal->GetBinContent(iBin+1);
          double fraction      = total ? hJpsiVsPt->GetBinContent(iBin+1) / total : 0 ;
          double fractionError = total ? hJpsiVsPt->GetBinError(iBin+1)   / total  : 1;
        
          gFractionVsPt->SetPoint(      iPoint, pt, fraction);
          gFractionVsPt->SetPointError( iPoint, 0., fractionError);
        }
        ++iPoint;
      }
    }  // end loop over pT bins
    
    stylePlot(gYieldVsPt, color, style, lineStyle, title, "#it{p}_{T} (GeV/c)", "yield (a.u.)", 0.,maxPt);
    drawInCanvas(cYieldVsPt, gYieldVsPt, optionY);
    lpy->AddEntry(gYieldVsPt, name.Data(), "LP");
    
    if(cutBin) {
      stylePlot(gFractionVsPt, color, style, lineStyle, title, "#it{p}_{T} (GeV/c)",
                "Fraction of total yield", 0., maxPt, 0., 1.);
      drawInCanvas(cFractionVsPt, gFractionVsPt, optionF);
      lpf->AddEntry(gFractionVsPt, name.Data(), "LP");
    }
    
    
    // Rapidity
    TH1D* hJpsiVsY = hJpsi->ProjectionX();
    
    if(!cutBin) hJpsiVsYTotal = (TH1D*) hJpsiVsY->Clone();
    
    iPoint = 0;
    for(unsigned int iBin=0; iBin<nRapidityBins; ++iBin){
      double maxYieldVsY = hJpsiVsY->GetMaximum();
      double yield       = maxYieldVsY ? hJpsiVsY->GetBinContent(iBin+1) / maxYieldVsY : 0.;
      double yieldError  = maxYieldVsY ? hJpsiVsY->GetBinError(iBin+1)   / maxYieldVsY : 0.;
      double y = hJpsiVsY->GetBinCenter(iBin+1);
      
      if(yield) { 
        gYieldVsY->SetPoint(      iPoint, y,  yield);
        gYieldVsY->SetPointError( iPoint, 0., yieldError);
        
        if(cutBin) {
          double total = hJpsiVsYTotal->GetBinContent(iBin+1);
          double fraction      = total ? hJpsiVsY->GetBinContent(iBin+1) / total : 0;
          double fractionError = total ? hJpsiVsY->GetBinError(iBin+1)   / total : 1;
        
          gFractionVsY->SetPoint(      iPoint, y,  fraction);
          gFractionVsY->SetPointError( iPoint, 0., fractionError);
        }
        ++iPoint;
      }
    }  // end loop over rapidity bins
    
    stylePlot(gYieldVsY, color, style, lineStyle, title, "y", "yield (a.u.)", -6., 6.);
    drawInCanvas(cYieldVsY, gYieldVsY, optionY);
    lyy->AddEntry(gYieldVsY, name.Data(), "LP");
    
    if(cutBin) {
      stylePlot(gFractionVsY, color, style, lineStyle, title, "y", "Fraction of total yield", -6., 6., 0.,1.);
      drawInCanvas(cFractionVsY, gFractionVsY, optionF);
      lyf->AddEntry(gFractionVsY, name.Data(), "LP");
    }
    
    gFractionVsMult->Write(Form("fractonVsMult_%d", iCut));
    gFractionVsPt  ->Write(Form("fractionVsPt_%d",  iCut));
    gFractionVsY   ->Write(Form("fractionVsY_%d",   iCut));
      
    gYieldVsMult   ->Write(Form("yieldVsMult_%d", iCut));
    gRatioToLin    ->Write(Form("ratioToLin_%d",  iCut));
    gYieldVsPt     ->Write(Form("yieldVsPt_%d",   iCut));
    gYieldVsY      ->Write(Form("yieldVsY_%d",    iCut));
  }
  
  
  lmy->AddEntry(linear, "linear", "L");
  linear->Write("linear");
  nullPt->Write("nullPt");
  onePt ->Write("onePt");

  nullY->Write("nullY");
  oneY ->Write("oneY");

  lmy->Write("legend");
  lpy->Write("legend2");
  out.Close();


  drawInCanvas(cYieldVsMult,    lmy);
  drawInCanvas(cRatioToLin,     lmy);
  drawInCanvas(cFractionVsMult, lmf);

  save(cYieldVsMult,    folder, Form("yieldVsMult_%s",    fn_add.Data()));
  save(cFractionVsMult, folder, Form("fractionVsMult_%s", fn_add.Data()));
  save(cRatioToLin,     folder, Form("ratioToLin_%s",     fn_add.Data()));

  drawInCanvas(cYieldVsPt,    lpy);
  drawInCanvas(cFractionVsPt, lpf);
  save(cYieldVsPt,    folder, Form("yieldVsPt_%s",    fn_add.Data()));
  save(cFractionVsPt, folder, Form("fractionVsPt_%s", fn_add.Data()));

  drawInCanvas(cYieldVsY,    lyy);
  drawInCanvas(cFractionVsY, lyf);
  save(cYieldVsY,    folder, Form("yieldVsY_%s",    fn_add.Data()));
  save(cFractionVsY, folder, Form("fractionVsY_%s", fn_add.Data()));
}


void crash(TString message) {
  cout << message << endl;
  exit(1);
}


void drawInCanvas(TCanvas* c, TObject* o, TString option) {
  c->cd();
  o->Draw(option);
}


TString normalized(TString str) {
   return "#frac{"  + str + "}{#LT" + str + "#GT}";
}


TCanvas* makeCanvas(TString name, bool logx, bool logy) {

  TCanvas* c = new TCanvas(name.Data(), name.Data(), 800, 600);

  c->SetTickx();
  c->SetTicky();

  c->SetTopMargin(0.08);
  c->SetRightMargin(0.02);
  c->SetLeftMargin(0.18);
  c->SetBottomMargin(0.2);

  if(logx) c->SetLogx();
  if(logy) c->SetLogy();
  return c;
}


TLegend* makeLegend(TString header, int xPos, int yPos) {
  double x1, x2, y1, y2;

  switch(xPos) {
    case  1: x1 = 0.55; x2 = 0.8;  break;
    case  2: x1 = 0.4;  x2 = 0.8;  break;
    default: x1 = 0.2;  x2 = 0.65; break;
  }
  switch(yPos) {
    case  1: y1 = 0.21; y2 = 0.55; break;
    default: y1 = 0.55; y2 = 0.89; break;
  }

  TLegend* l = new TLegend(x1, y1, x2, y2);
  l->SetHeader(header.Data());
  l->SetTextSize(0.05);
  l->SetBorderSize(0);
  l->SetFillColor(0);
  l->SetFillStyle(0);
  return l;
}


void save(TCanvas* c, TString folder, TString filename) {
  c->SaveAs(Form("%s/%s.pdf", folder.Data(), filename.Data()));
  c->SaveAs(Form("%s/%s.png", folder.Data(), filename.Data()));
}


void styleFunction(TF1* g, Color_t color, int lineStyle, TString title, TString xtitle, TString ytitle,
                   double xMin, double xMax, double yMin, double yMax) {

  g->SetLineColor(color);
  g->SetLineStyle(lineStyle);
  g->SetLineWidth(2);

  g->SetTitle(Form("PYTHIA 8.3, p--Pb #sqrt{#it{s}} = 5 TeV, %s", title.Data()));
  g->GetXaxis()->SetTitle(xtitle.Data());
  g->GetYaxis()->SetTitle(ytitle.Data());
  g->GetXaxis()->SetTitleOffset(1.5);
  g->GetYaxis()->SetTitleOffset(1.2);
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetLabelSize(0.05);
  if(xMax) g->GetXaxis()->SetRangeUser(xMin, xMax);
  if(yMax) g->GetYaxis()->SetRangeUser(yMin, yMax);
}


void stylePlot(TGraphErrors* g, Color_t color, int style, int lineStyle, TString title, TString xtitle,
               TString ytitle, double xMin, double xMax, double yMin, double yMax) {

  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(style);
  g->SetMarkerSize(1.4);
//   g->SetMarkerSize(2);
  g->SetLineStyle(lineStyle);
  g->SetLineWidth(2);

//   g->SetTitle(title.Data());
  g->SetTitle(Form("PYTHIA 8.3, p--Pb #sqrt{#it{s}} = 5 TeV, %s", title.Data()));
  g->GetXaxis()->SetTitle(xtitle.Data());
  g->GetYaxis()->SetTitle(ytitle.Data());
  g->GetXaxis()->SetTitleOffset(1.5);
  g->GetYaxis()->SetTitleOffset(1.2);
  g->GetXaxis()->SetTitleSize(0.055);
  g->GetYaxis()->SetTitleSize(0.055);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetLabelSize(0.05);

  if(xMax) g->GetXaxis()->SetRangeUser(xMin, xMax);
  if(yMax) g->GetYaxis()->SetRangeUser(yMin, yMax);
}
