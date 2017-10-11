///////////////////////////////////
//
// Usage: 
//   > root -l makeHistos.C 
//
// Comments:
//   - Systematic histograms contain varied yields (in absolute values, not relative) 
///////////////////////////////////

void writeHistoToFile(TH1F* h, TString fileName)
{
  TFile* f = new TFile(fileName,"recreate");
  h->Write();
  f->Write();
  f->Close();
}

void setHistoStyle(TH1F* h, int color, int lineStyle=1)
{
  h->SetLineStyle(lineStyle);
  h->SetLineColor(color);
}

void makeHistos()
{
  gROOT->SetStyle("Plain");

  int Nbins=2;
  float binMin=0;
  float binMax=2;

  // signal
  TH1F* hSig = new TH1F("histo","hSig",Nbins,binMin,binMax);
  TH1F* hSigSyst1up = new TH1F("histo","hSigSyst1up",Nbins,binMin,binMax);
  TH1F* hSigSyst1down = new TH1F("histo","hSigSyst1down",Nbins,binMin,binMax);
  TH1F* hSigSyst2up = new TH1F("histo","hSigSyst2up",Nbins,binMin,binMax);
  TH1F* hSigSyst2down = new TH1F("histo","hSigSyst2down",Nbins,binMin,binMax);
  TH1F* hSigSyst5up = new TH1F("histo","hSigSyst5up",Nbins,binMin,binMax);
  TH1F* hSigSyst5down = new TH1F("histo","hSigSyst5down",Nbins,binMin,binMax);

  // ttbar
  TH1F* hBg1 = new TH1F("histo","hBg1",Nbins,binMin,binMax);
  TH1F* hBg1Syst1up = new TH1F("histo","hBg1Syst1up",Nbins,binMin,binMax);
  TH1F* hBg1Syst1down = new TH1F("histo","hBg1Syst1down",Nbins,binMin,binMax);
  TH1F* hBg1Syst2up = new TH1F("histo","hBg1Syst2up",Nbins,binMin,binMax);
  TH1F* hBg1Syst2down = new TH1F("histo","hBg1Syst2down",Nbins,binMin,binMax);

  // ttWZ
  TH1F* hBg2 = new TH1F("histo","hBg2",Nbins,binMin,binMax);
  TH1F* hBg2Syst1up = new TH1F("histo","hBg2Syst1up",Nbins,binMin,binMax);
  TH1F* hBg2Syst1down = new TH1F("histo","hBg2Syst1down",Nbins,binMin,binMax);
  TH1F* hBg2Syst3up = new TH1F("histo","hBg2Syst3up",Nbins,binMin,binMax);
  TH1F* hBg2Syst3down = new TH1F("histo","hBg2Syst3down",Nbins,binMin,binMax);

  // QCD
  TH1F* hBg3 = new TH1F("histo","hBg3",Nbins,binMin,binMax);
  TH1F* hBg3Syst2up = new TH1F("histo","hBg3Syst2up",Nbins,binMin,binMax);
  TH1F* hBg3Syst2down = new TH1F("histo","hBg3Syst2down",Nbins,binMin,binMax);
  TH1F* hBg3Syst4up = new TH1F("histo","hBg3Syst4up",Nbins,binMin,binMax);
  TH1F* hBg3Syst4down = new TH1F("histo","hBg3Syst4down",Nbins,binMin,binMax);

  // WW
  TH1F* hBg4 = new TH1F("histo","hBg4",Nbins,binMin,binMax);
  TH1F* hBg4Syst1up = new TH1F("histo","hBg4Syst1up",Nbins,binMin,binMax);
  TH1F* hBg4Syst1down = new TH1F("histo","hBg4Syst1down",Nbins,binMin,binMax);

  TH1F* hData = new TH1F("histo","hData",Nbins,binMin,binMax);
  hData->SetMarkerStyle(8);
  hData->SetMarkerSize(2);


  // Fill signal
  hSig->SetBinContent(1,5);
  hSig->SetBinError(1,1);
  hSig->SetBinContent(2,15.5);
  hSig->SetBinError(2,2.1);

  hSigSyst1up   ->SetBinContent(1,5.3);
  hSigSyst1down ->SetBinContent(1,4.7);
  hSigSyst2up   ->SetBinContent(1,5.25);
  hSigSyst2down ->SetBinContent(1,3.5);
  hSigSyst5up   ->SetBinContent(1,5.15);
  hSigSyst5down ->SetBinContent(1,4.95);
  hSigSyst1up   ->SetBinContent(2,14.57);
  hSigSyst1down ->SetBinContent(2,16.43);
  hSigSyst2up   ->SetBinContent(2,15.5);
  hSigSyst2down ->SetBinContent(2,15.5);
  hSigSyst5up   ->SetBinContent(2,16.74);
  hSigSyst5down ->SetBinContent(2,13.95);

  // Fill ttbar
  hBg1->SetBinContent(1,25);
  hBg1->SetBinError(1,7);
  hBg1->SetBinContent(2,60);
  hBg1->SetBinError(2,3);
  
  hBg1Syst1up   -> SetBinContent(1,27.5);
  hBg1Syst1down -> SetBinContent(1,17.5);
  hBg1Syst2up   -> SetBinContent(1,32.5);
  hBg1Syst2down -> SetBinContent(1,20);
  hBg1Syst1up   -> SetBinContent(2,65.4);
  hBg1Syst1down -> SetBinContent(2,57);
  hBg1Syst2up   -> SetBinContent(2,60);
  hBg1Syst2down -> SetBinContent(2,60);

  // Fill ttWZ
  hBg2->SetBinContent(1,25);
  hBg2->SetBinError(1,12);
  hBg2->SetBinContent(2,55);
  hBg2->SetBinError(2,9.4);

  hBg2Syst1up   -> SetBinContent(1,30);
  hBg2Syst1down -> SetBinContent(1,23.75);
  hBg2Syst3up   -> SetBinContent(1,23.5);
  hBg2Syst3down -> SetBinContent(1,28.75);
  hBg2Syst1up   -> SetBinContent(2,56.1);
  hBg2Syst1down -> SetBinContent(2,52.25);
  hBg2Syst3up   -> SetBinContent(2,60.5);
  hBg2Syst3down -> SetBinContent(2,52.8);

  // Fill QCD
  hBg3->SetBinContent(1,33);
  hBg3->SetBinError(1,3.5);
  hBg3->SetBinContent(2,42);
  hBg3->SetBinError(2,7.6);

  hBg3Syst2up   -> SetBinContent(1,29.7);
  hBg3Syst2down -> SetBinContent(1,42.9);
  hBg3Syst4up   -> SetBinContent(1,13.2);
  hBg3Syst4down -> SetBinContent(1,52.8);
  hBg3Syst2up   -> SetBinContent(2,36.12);
  hBg3Syst2down -> SetBinContent(2,52.5);
  hBg3Syst4up   -> SetBinContent(2,25.2);
  hBg3Syst4down -> SetBinContent(2,63.42);

  // Fill WW
  hBg4 -> SetBinContent(2,3);
  hBg4 -> SetBinError(2,0);

  hBg4Syst1up -> SetBinContent(2,3.24);
  hBg4Syst1down -> SetBinContent(2,2.79);

  // Fill data
  hData -> SetBinContent(1,100);
  hData -> SetBinContent(2,150);

  // Write histos
  writeHistoToFile(hSig,"histos/fSig.root");
  writeHistoToFile(hSigSyst1up,"histos/fSigSyst1up.root");
  writeHistoToFile(hSigSyst1down,"histos/fSigSyst1down.root");
  writeHistoToFile(hSigSyst2up,"histos/fSigSyst2up.root");
  writeHistoToFile(hSigSyst2down,"histos/fSigSyst2down.root");
  writeHistoToFile(hSigSyst5up,"histos/fSigSyst5up.root");
  writeHistoToFile(hSigSyst5down,"histos/fSigSyst5down.root");

  writeHistoToFile(hBg1,"histos/fBg1.root");
  writeHistoToFile(hBg1Syst1up,"histos/fBg1Syst1up.root");
  writeHistoToFile(hBg1Syst1down,"histos/fBg1Syst1down.root");
  writeHistoToFile(hBg1Syst2up,"histos/fBg1Syst2up.root");
  writeHistoToFile(hBg1Syst2down,"histos/fBg1Syst2down.root");

  writeHistoToFile(hBg2,"histos/fBg2.root");
  writeHistoToFile(hBg2Syst1up,"histos/fBg2Syst1up.root");
  writeHistoToFile(hBg2Syst1down,"histos/fBg2Syst1down.root");
  writeHistoToFile(hBg2Syst3up,"histos/fBg2Syst3up.root");
  writeHistoToFile(hBg2Syst3down,"histos/fBg2Syst3down.root");

  writeHistoToFile(hBg3,"histos/fBg3.root");
  writeHistoToFile(hBg3Syst2up,"histos/fBg3Syst2up.root");
  writeHistoToFile(hBg3Syst2down,"histos/fBg3Syst2down.root");
  writeHistoToFile(hBg3Syst4up,"histos/fBg3Syst4up.root");
  writeHistoToFile(hBg3Syst4down,"histos/fBg3Syst4down.root");

  writeHistoToFile(hBg4,"histos/fBg4.root");
  writeHistoToFile(hBg4Syst1up,"histos/fBg4Syst1up.root");
  writeHistoToFile(hBg4Syst1down,"histos/fBg4Syst1down.root");

  writeHistoToFile(hData,"histos/fData.root");
}
