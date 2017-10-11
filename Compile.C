void Compile() {
  gROOT->LoadMacro("OTHRdmGenerator.C+");
  gROOT->LoadMacro("OTHSystematics.C+");
  gROOT->LoadMacro("OTHAlgorithms.C+");
  gROOT->LoadMacro("OTHBase.C+");
  gROOT->LoadMacro("OTHPdfGenerator.C+");
  gROOT->LoadMacro("OTHSingleSyst.C+");
  gROOT->LoadMacro("OTHYieldWithUncert.C+");
  gROOT->LoadMacro("OTHSample.C+");
  gROOT->LoadMacro("OTHObserved.C+");
  gROOT->LoadMacro("OTHMuVsObs.C+");
  gROOT->LoadMacro("OTHChannel.C+");
  gROOT->LoadMacro("OTHShapeSyst.C+");
  gROOT->LoadMacro("OTHShape.C+");
  gROOT->LoadMacro("OpTHyLiC.C+");
}
