///////////////////////////////////////////////////////////////////////////////////
// This file is part the software OpTHyLiC
// Copyright © 2015 Laboratoire de Physique Corpusculaire de Clermont-Ferrand
// Developpers: David Calvet, Emmanuel Busato, Timothée Theveneaux-Pelzer
// Contact: opthylic@in2p3.fr
//
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
using namespace std;

#include "TH1.h"
#include "TFile.h"

#include "OTHRdmGenerator.h"
#include "OTHSystematics.h"
#include "OTHPdfGenerator.h"
#include "OTHShape.h"
#include "OTHShapeSyst.h"

#include "OpTHyLiC.h"
using namespace OTH;

OpTHyLiC::OpTHyLiC(const SystType systInterpExtrapStyle, const StatType statSampling, const int RandomEngineType, const int seed,const CombType systCombinationType) :
  Base(),
  m_pRdmGen(0),
  m_pSyste(0),
  m_pStatSampling(0),
  m_pChannels(),
  m_sigStrength(1),
  m_sumMu(0),
  m_nbMu(0),
  m_pHs(nbHistos,0),
  m_muObs(),
  m_muObsInterpol()
{
  // using pseudo-random number generator provided by TRandom3 class (default)
  if (RandomEngineType==TR3) {
    m_pRdmGen = new RdmGenerator_TR3(seed);
    cout << "OpTHyLiC Info: using pseudo-random number generator implemented in TRandom3 class" <<endl;
  }
#if defined CPP11
  // using pseudo-random number generators provided by C++11 standard library
  else if (RandomEngineType==STD_minstd_rand) {
    m_pRdmGen = new RdmGenerator_STD<std::minstd_rand>(seed);
    cout << "OpTHyLiC Info: using pseudo-random number generator implemented in std::minstd_rand class" <<endl;
  }
  else if (RandomEngineType==STD_minstd_rand0) {
    m_pRdmGen = new RdmGenerator_STD<std::minstd_rand0>(seed);
    cout << "OpTHyLiC Info: using pseudo-random number generator implemented in std::minstd_rand0 class" <<endl;
  }
  else if (RandomEngineType==STD_mt19937) {
    m_pRdmGen = new RdmGenerator_STD<std::mt19937>(seed);
    cout << "OpTHyLiC Info: using pseudo-random number generator implemented in std::mt19937 class" <<endl;
  }
  else if (RandomEngineType==STD_mt19937_64) {
    m_pRdmGen = new RdmGenerator_STD<std::mt19937_64>(seed);
    cout << "OpTHyLiC Info: using pseudo-random number generator implemented in std::mt19937_64 class" <<endl;
  }
  else if (RandomEngineType==STD_ranlux24_base) {
    m_pRdmGen = new RdmGenerator_STD<std::ranlux24_base>(seed);
    cout << "OpTHyLiC Info: using pseudo-random number generator implemented in std::ranlux24_base class" <<endl;
  }
  else if (RandomEngineType==STD_ranlux48_base) {
    m_pRdmGen = new RdmGenerator_STD<std::ranlux48_base>(seed);
    cout << "OpTHyLiC Info: using pseudo-random number generator implemented in std::ranlux48_base class" <<endl;
  }
  else if (RandomEngineType==STD_ranlux24) {
    m_pRdmGen = new RdmGenerator_STD<std::ranlux24>(seed);
    cout << "OpTHyLiC Info: using pseudo-random number generator implemented in std::ranlux24 class" <<endl;
  }
  else if (RandomEngineType==STD_ranlux48) {
    m_pRdmGen = new RdmGenerator_STD<std::ranlux48>(seed);
    cout << "OpTHyLiC Info: using pseudo-random number generator implemented in std::ranlux48 class" <<endl;
  }
  else if (RandomEngineType==STD_knuth_b) {
    m_pRdmGen = new RdmGenerator_STD<std::knuth_b>(seed);
    cout << "OpTHyLiC Info: using pseudo-random number generator implemented in std::knuth_b class" <<endl;
  }
#endif
  else {
    cerr << "OpTHyLiC Error ! Unknown random generator engine type" << endl;
    throw runtime_error("Unknown random generator engine type !");
  }
  
  // treatment of systematic uncertainties
  m_pSyste = new Systematics(m_pRdmGen, systInterpExtrapStyle);

  // treatment of statistical uncertainties
  m_pStatSampling = new PdfGenerator(m_pRdmGen, statSampling);

  // treatment of combination of systematics uncertainties
  if(systCombinationType!=CombAutomatic && systCombinationType!=CombAdditive && systCombinationType!=CombMultiplicative) {
    cerr << "OpTHyLiC Error ! Unknown combination type for systematics " << endl;
    throw runtime_error("Unknown combination type !");
  }
  else {
    if(systCombinationType==CombAutomatic) { // default behaviour for combination type
      m_additiveSystComb=(systInterpExtrapStyle==SystMclimit || systInterpExtrapStyle==SystLinear);
    }
    else m_additiveSystComb=(systCombinationType==CombAdditive);
  }
}

OpTHyLiC::~OpTHyLiC()
{
  delete m_pRdmGen;
  delete m_pSyste;
  delete m_pStatSampling;
  for(unsigned int i=0 ; i<m_pChannels.size() ; ++i) {
    delete m_pChannels[i];
  }
  // do not delete histos, belong to ROOT
}

bool OpTHyLiC::isShape(const string &fileName) const
{
  ifstream in(fileName.c_str());
  if (!in) {
    cerr << "ERROR ! Unable to open file '" << fileName << "' !" << endl;
  }
  for(; !in.eof() ;) {
    string line;
    if (!getline(in,line)) break;
    if (!line.empty() && line[0]!='#') {
      istringstream istr(line);
      if(line.find(".root") != string::npos)
	return true;
    }
  }
  return false;
}

void OpTHyLiC::makeInputsFromShapes(const string &channelName, const string &fileName,const bool removeFiles) 
{
  ifstream in(fileName.c_str());
  if (!in) {
    cerr << "ERROR ! Unable to open file '" << fileName << "' !" << endl;
    return;
  }

  enum {None,Background,Signal,Setup};
  int type=None;
  bool signalSet=false;

  string directory("");
  string histoName("");

  // Objects that contain shapes and that are used afterward to make inputs
  Shape* dataShape = NULL;
  Shape* sigShape = NULL;
  deque<Shape*> bgShapes;
  string channelNameLaTeX("");

  // Parse input file to fill dataShape, sigShape and bgShapes
  for(; !in.eof() ;) {
    string line;
    if (!getline(in,line)) break;
    if (!line.empty() && line[0]!='#') {
      istringstream istr(line);
      string keyword;
      string name="";
      string file1="";
      string file2="";
      istr >> keyword >> name >> file1 >> file2;
      //cout << "printing in shapes: " << keyword << "  " << name << "  " << file1 << "  " << file2 << endl;
      if ("+setup"==keyword) {
	if(type==Background || type==Signal) {
	  cerr << "ERROR ! +setup should be before any signal or background" << endl;
	  return;
	}
	type=Setup; 
      } else if(type==Setup && ".directory"==keyword) {
	directory=name;
	if(directory!="")
	  directory+="/";
      } else if(type==Setup && ".histoName"==keyword) {
	histoName=name;
      } else if("+bg"==keyword || "+sig"==keyword || "+data"==keyword) {
	string histoNameLocal(histoName);
	if("+data"==keyword) 
	  file1 = name;
	if(histoNameLocal=="") {
	  size_t npos1=file1.find_last_of("(");
	  size_t npos2=file1.find_last_of(")");
	  if(npos1!=string::npos && npos2!=string::npos) {
	    histoNameLocal=file1.substr(npos1+1,npos2-npos1-1);
	  } else {
	    cerr << "ERROR ! histoName not set -> set it !" << endl;
	    return;
	  }
	  file1=file1.substr(0,npos1);
	}
	TFile f(Form("%s%s",directory.c_str(),file1.c_str()),"read");
	TH1* h = (TH1*) f.Get(histoNameLocal.c_str());	
	h->SetDirectory(0);

	if("+bg"==keyword) {
	  type=Background;
	  bgShapes.push_back(new Shape(name,h));
	} else if ("+sig"==keyword) {
	  type=Signal;
	  sigShape = new Shape(name,h);
	  if (signalSet) {
	    cerr << "WARNING ! Signal already set while setting new signal '" << name << "' !" << endl;
	  }
	  signalSet=true;
	} else if ("+data"==keyword) {
	  dataShape = new Shape("data",h);
	} 
      } else if (".syst"==keyword) {
	if (!(type==Signal || type==Background)) {
	  cerr << "SYNTAX ERROR ! Trying to define a systematic uncertainty '" << name 
	       << "' before defining any sample !" << endl;
	} else if (line.find(".root") != string::npos){ // shape syst
	  string histoNameLocal1(histoName);
	  string histoNameLocal2(histoName);
	  if(histoNameLocal1=="") {
	    size_t npos1=file1.find_last_of("(");
	    size_t npos2=file1.find_last_of(")");
	    if(npos1!=string::npos && npos2!=string::npos) {
	      histoNameLocal1=file1.substr(npos1+1,npos2-npos1-1);
	    } else {
	      cerr << "ERROR ! histoName1 not set -> set it !" << endl;
	      return;
	    }
	    file1=file1.substr(0,npos1);
	  }
	  if(histoNameLocal2=="") {
	    size_t npos1=file2.find_last_of("(");
	    size_t npos2=file2.find_last_of(")");
	    if(npos1!=string::npos && npos2!=string::npos) {
	      histoNameLocal2=file2.substr(npos1+1,npos2-npos1-1);
	    } else {
	      cerr << "ERROR ! histoName1 not set -> set it !" << endl;
	      return;
	    }
	    file2=file2.substr(0,npos1);
	  }

	  TFile fUp(Form("%s%s",directory.c_str(),file1.c_str()),"read");
	  TFile fDown(Form("%s%s",directory.c_str(),file2.c_str()),"read");
	  TH1* hUp = (TH1*) fUp.Get(histoNameLocal1.c_str());
	  TH1* hDown = (TH1*) fDown.Get(histoNameLocal2.c_str());
	  hUp->SetDirectory(0);
	  hDown->SetDirectory(0);
	  ShapeSyst* shapeSyst = new ShapeSyst(name,hUp,hDown);
	  if (Background==type) {
	    bgShapes.back()->addShapeSystematic(shapeSyst);
	  } else if (Signal==type) {
	    sigShape->addShapeSystematic(shapeSyst);
	  } else {
	    cerr << "WHAT THE HELL !" << endl;
	  }
	} else { // global syst
	  istringstream s1(file1);
	  float systUp;
	  s1 >> systUp;
	  istringstream s2(file2);
	  float systDown;
	  s2 >> systDown;
	  SingleSyst* globalSyst = new SingleSyst(name,0,systDown,systUp);
	  if (Background==type) {
	    bgShapes.back()->addGlobalSystematic(globalSyst);
	  } else if (Signal==type) {
	    sigShape->addGlobalSystematic(globalSyst);
	  } else {
	    cerr << "WHAT THE HELL !" << endl;
	  }
	}
      } else if (".nameLaTeX"==keyword) {
	const string nameLaTeX=line.substr(line.find_first_of(' ')+1);
	if (None==type || Setup==type) {
	  cerr << "SYNTAX ERROR ! Trying to define a LaTeX name before defining any sample !" << endl;
	} else if (Background==type) bgShapes.back()->setNameLaTeX(nameLaTeX);
	else if (Signal==type) sigShape->setNameLaTeX(nameLaTeX);
	else cerr << "WHAT THE HELL !" << endl;
      } else if ("+nameLaTeX"==keyword) {
	channelNameLaTeX=line.substr(line.find_first_of(' ')+1);
      } else {
	cerr << "Unknown line '" << line << "' !" << endl;
      }
    }
  }
  /*
  dataShape->print(false);
  sigShape->print();
  for(unsigned int i=0; i<bgShapes.size(); ++i) bgShapes[i]->print();
  */

  // Sanity check
  int nBinsSig = 0;
  if(sigShape) {
    nBinsSig = sigShape->getNbins();
    if(dataShape) {
      if(dataShape->getNbins() != nBinsSig) {
	cerr << "ERROR ! Number of bins not the same for data and signal" << endl;
	return;
      }
    }
    for(unsigned int i=0; i<bgShapes.size(); ++i) {
      if(bgShapes[i]->getNbins() != nBinsSig) {
	cerr << "ERROR ! Number of bins not the same for data and background " << bgShapes[i]->getName() << endl;
	return;
      }
    }
  }
  else {
    cerr << "OpTHyLiC Error ! No signal has been defined in the input file !" << endl;
    throw runtime_error("sigShape not set !");
  }

  // Dump one input per bin from dataShape, sigShape and bgShapes objects
  for(int i=1; i<=nBinsSig; ++i) {
    string fileNameThisBin = fileName.substr(0,fileName.find_last_of('.'));
    string fileNameSuffix = fileName.substr(fileName.find_last_of('.')+1);
    const char* fileNameBin = Form("%s_bin%i.%s",fileNameThisBin.c_str(),i,fileNameSuffix.c_str());
    ofstream of(fileNameBin);

    // Channel name
    of << "+nameLaTeX " << channelNameLaTeX << " (bin " << i << ")" << endl ;
    
    // Backgrounds
    bool bgYieldIsNull=true;
    for(unsigned int j=0; j<bgShapes.size(); ++j) {
      if(bgShapes[j]->getBinContent(i) != 0 || bgShapes[j]->getBinError(i) != 0) {
	bgYieldIsNull=false;
	of << endl << "+bg ";
	bgShapes[j]->writeInputFile(of,i);
      }
      else {
	cout << "Warning: bin " << i << " (value=" << bgShapes[j]->getBinCenter(i) << ") of background " << bgShapes[j]->getName() << " has 0 events and no statistical uncertainty" << endl;
      }
    }
    
    // Signal
    bool sigYieldIsNull=true;
    if(sigShape->getBinContent(i) != 0 || sigShape->getBinError(i) != 0) {
      sigYieldIsNull=false;
      of << endl << "+sig ";
      sigShape->writeInputFile(of,i);
    }
    else {
      cout << "Warning: bin " << i << " (value=" << sigShape->getBinCenter(i) << ") of signal has 0 events and no statistical uncertainty" << endl;
    }

    // Data
    if(dataShape) {
      of << endl 
	 << "+data " << dataShape->getBinContent(i) << endl;
    }
    
    if(bgYieldIsNull==false && sigYieldIsNull==false) {
      const char* channelNameBin = Form("%s_bin%i",channelName.c_str(),i);
      m_pChannels.push_back(new Channel(channelNameBin,*m_pSyste,*m_pStatSampling));
      m_pChannels.back()->addSamples(fileNameBin);
      m_pChannels.back()->setCombinationType(m_additiveSystComb);
      m_pChannels.back()->setConfLevel(m_confLevel);
      if(removeFiles)
	system(Form("rm -f %s",fileNameBin));
    }
    else {
      cout << "Warning: yield and statistical uncertainty of signal and/or total background in bin " << i << " are equal to 0" << endl;
      cout << "            -> this bin is ignored" << endl;
      system(Form("rm -f %s",fileNameBin));
    }
  }

  if(dataShape) delete dataShape;
  if(sigShape) delete sigShape;
  for(unsigned int i=0; i<bgShapes.size(); ++i) {
    if(bgShapes[i]) delete bgShapes[i];
  }
}

void OpTHyLiC::setConfLevel(const double cl)
{
  Base::setConfLevel(cl);
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    m_pChannels[c]->setConfLevel(cl);
  }  
}

unsigned int OpTHyLiC::addChannel(const string &name)
{
  // add new channel
  m_pChannels.push_back(new Channel(name,*m_pSyste,*m_pStatSampling));
  m_pChannels.back()->setCombinationType(m_additiveSystComb);
  m_pChannels.back()->setConfLevel(m_confLevel);
  return m_pChannels.size()-1;
}

unsigned int OpTHyLiC::addChannel(const string &name,const string &fileName,const bool removeFiles)
{
  if(isShape(fileName)) {
    // add as many channels as there are bins in the input histogram
    makeInputsFromShapes(name,fileName,removeFiles);
  } else {
    // add one new channel
    m_pChannels.push_back(new Channel(name,*m_pSyste,*m_pStatSampling));
    m_pChannels.back()->addSamples(fileName);
    m_pChannels.back()->setCombinationType(m_additiveSystComb);
    m_pChannels.back()->setConfLevel(m_confLevel);
  }
  return m_pChannels.size()-1;
}

Channel* OpTHyLiC::getChannel(const unsigned int iChannel)
{
  if (iChannel<m_pChannels.size()) return m_pChannels[iChannel];
  throw runtime_error("Unknown channel index !");
}

Channel* OpTHyLiC::getChannel(const string name)
{
  for(deque<Channel*>::iterator it = m_pChannels.begin() ; it!=m_pChannels.end() ; ++it) {
    if ((*it)->getName()==name) return *it;
  }
  throw runtime_error("Unknown channel name !");
}

void OpTHyLiC::setSigStrength(const double mu)
{
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    m_pChannels[c]->setSigStrength(mu);
  }
}

double OpTHyLiC::computeLLRdata() const
{
  double llr=0;
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    llr+=m_pChannels[c]->computeLLRdata();
  }
  return llr;
}

void OpTHyLiC::generateDistrLLR(const int nbExp)
{
  // resetting
  for(unsigned int h=0 ; h<m_pHs.size() ; ++h) {
    if (m_pHs[h]) {
      delete m_pHs[h];
      m_pHs[h]=0;
    }
  }
  double llrMini=0,llrMaxi=0;
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    double llrMin,llrMax;
    m_pChannels[c]->initDistrLLR(llrMin,llrMax);
    llrMini+=llrMin;
    llrMaxi+=llrMax;
  }

  if (nbExp<1) return;

  // creation of histograms to store distributions
  const double llrMin=llrMini;
  const double llrMax=llrMaxi;
  m_pHs[hLLRb]=new TH1F("LLRb",";LLR;Probability",10000,llrMin,llrMax);
  m_pHs[hLLRsb]=new TH1F("LLRsb",";LLR;Probability",10000,llrMin,llrMax);

  // loop on all pseudo-experiments
  for(int i=0 ; i<nbExp ; ++i) {
    // systematic uncertainties variations
    m_pSyste->variate();

    // compute test-statistic
    double sumLLRb=0,sumLLRsb=0;
    for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
      double llrB,llrSB;
      m_pChannels[c]->generateSinglePseudoExp(llrB,llrSB);
      sumLLRb+=llrB;
      sumLLRsb+=llrSB;
    }

    m_pHs[hLLRb]->Fill(sumLLRb);
    m_pHs[hLLRsb]->Fill(sumLLRsb);
  }
  
  // normalization of distributions
  m_pHs[hLLRb]->Scale(1/static_cast<float>(nbExp));
  m_pHs[hLLRsb]->Scale(1/static_cast<float>(nbExp));
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    m_pChannels[c]->endDistrLLR(nbExp);
  }
}

double OpTHyLiC::computeCLsData() const
{
  if (!m_pHs[hLLRsb] || !m_pHs[hLLRb]) {
    cout << "ERROR ! No LLR distribution found, use generateDistrLLR first..." << endl;
    return -1;
  }

  return Algorithms::computeCLs(m_pHs[hLLRsb],m_pHs[hLLRb],computeLLRdata());
}

double OpTHyLiC::pValueData() const
{
  if (!m_pHs[hLLRb]) {
    cout << "ERROR ! No background distribution found, use generateDistrLLR first..." << endl;
    return 0;
  }
  return m_pHs[hLLRb]->Integral(0,m_pHs[hLLRb]->FindBin(computeLLRdata()));
}

double OpTHyLiC::generateForCLs(const double mu,const int nbExp,const int type)
{
  setSigStrength(mu);
  generateDistrLLR(nbExp);
  if(LimObserved==type) {
    return computeCLsData();
  }
  else if(type>=LimExpectedP2sig && type<=LimExpectedM2sig) {
    return Algorithms::getCLsFromLLR(type,m_pHs[hLLRsb],m_pHs[hLLRb]);
  }
  else {
    throw runtime_error("Unknown limit type !");
  }
}

double OpTHyLiC::sigStrengthExclusion(const LimitType type,const int nbExp,double &cls,
				      const double muHint,const MethType method)
{
  double mu=0.5,muStep=3;
  if (muHint!=1) mu=muHint/2;
  else if(LimObserved==type) {
    Observed obs(m_pChannels.size());
    for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
      obs.set(c,m_pChannels[c]->getYieldData());
    }

    cout << "--------- Searching mu for obs = ( ";
    obs.print();
    cout << ") -----------" << endl;

    // if already two mus, interpolate to find first value of mu
    for(unsigned int c=0 ; c<m_muObsInterpol.size() && c<m_pChannels.size() ; ++c) {
      Observed obs1=obs;
      obs1.set(c,-1);
      const map<Observed,MuVsObs> &muObs1=m_muObsInterpol[c];
      map<Observed,MuVsObs>::const_iterator it=muObs1.find(obs1);
      if (it!=muObs1.end()) {
	double mu1=it->second.interpolateMu(m_pChannels[c]->getYieldData());
	if (mu1>0) {
	  mu=mu1;
	  muStep=1.2;
	  break;
	}
      }
    }

    // if no good mu, try average
    if (1==mu && m_nbMu>0) {
      mu=m_sumMu/m_nbMu;
      if (mu<1e-5) mu=1;
    }
  }

  return Algorithms::sigStrengthExclusion(*this,mu,muStep,nbExp,type,cls,m_confLevel,method!=MethDichotomy);
}

double OpTHyLiC::expectedSigStrengthExclusion(const int nbMu,const int nbExp)
{
  // resetting average mu
  m_sumMu=0;
  m_nbMu=0;

  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    m_pChannels[c]->saveYieldData();
  }

  // resetting list of mu values
  m_muObs.clear();
  m_muObsInterpol.resize(m_pChannels.size());
  Observed obs(m_pChannels.size());
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    m_pChannels[c]->setYieldDataToBkg();
    obs.set(c,m_pChannels[c]->getYieldData());
  }
  double cls=0;
  const double mu0=sigStrengthExclusion(LimObserved,nbExp,cls);  
  m_muObs[obs]=mu0;
  setMuVsObs(obs,mu0);

  // resetting histograms
  if (m_pExpMu) {
    delete m_pExpMu;
    m_pExpMu=0;
  }
  m_pExpMu=new TH1F("hmu",";Expected signal strength;Probability",10000,0,mu0*10);
  if (m_pCLs) {
    delete m_pCLs;
    m_pCLs=0;
  }
  m_pCLs=new TH1F("hCLs",";CL_{s};Entries",1000,(1-m_confLevel)*0.8,(1-m_confLevel)*1.2);
  m_pCLs->Fill(cls);

  // loop on background only pseudo-experiments
  for(int i=0 ; i<nbMu ; ++i) {
    // systematic uncertainties variations
    m_pSyste->variate();

    // generate pseudo-experiments
    for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
      obs.set(c,m_pChannels[c]->generateSinglePseudoData());
    }

    // search for mu_95 in map, otherwise compute it
    map<Observed,double>::const_iterator it=m_muObs.find(obs);
    double mu=0;
    if (it!=m_muObs.end()) mu=it->second;
    else {
      mu=sigStrengthExclusion(LimObserved,nbExp,cls);
      m_muObs[obs]=mu;
      setMuVsObs(obs,mu);
      m_pCLs->Fill(cls);
    } 
    m_pExpMu->Fill(mu);
    if ((i+1)%10000==0) cout << "---- Already " << i+1 << " mus computed" << endl;
  }
  m_pExpMu->Scale(1/static_cast<float>(nbMu));

  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    m_pChannels[c]->restoreYieldData();
  }

  vector<double> muQ=Algorithms::getQuantiles(m_pExpMu);
  return muQ[2];
}

void OpTHyLiC::printSamples() const
{
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    m_pChannels[c]->printSamples();
  }
}

void OpTHyLiC::printSyst() const
{
  m_pSyste->print();
}

void OpTHyLiC::printMuVsObs() const 
{
  /*
  cout << "-> printing mu Vs observation" << endl;
  for (auto it=m_muObs.begin(); it!=m_muObs.end(); ++it) {
    it->first.print();
    cout << ": " << it->second << endl;
  }
  */
}

void OpTHyLiC::createInputYieldTable(ostream &latex,const int precision) const
{
  createYieldTable(0,latex,precision);
}

void OpTHyLiC::createGeneratedYieldTable(ostream &latex,const int precision,const int nbExp) const
{
  createYieldTable(nbExp,latex,precision);
}

void OpTHyLiC::createSysteTables(ostream &latex,const string fileName,const int precision) const
{
  map<string,string> systeNames;
  ifstream in(fileName.c_str());
  if (!in) {
    cerr << "ERROR ! Unable to open file '" << fileName << "' !" << endl;
  } else {
    for(; !in.eof() ;) {
      string line;
      if (!getline(in,line)) break;
      if (!line.empty() && line[0]!='#') {
	istringstream istr(line);
	string name;
	istr >> name;
	systeNames[name]=line.substr(line.find_first_of(' ')+1);
      }
    }
    in.close();
  }

  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    const Sample &signal=m_pChannels[c]->getSigSample();
    const deque<Sample> &samples=m_pChannels[c]->getBkgSamples();
    latex << "\\begin{table}\\begin{center}" << endl
	  << "\\caption{List of relative systematic uncertainties (in \\%) for channel " 
	  << m_pChannels[c]->getNameLaTeX() << ".}" << endl
	  << "\\renewcommand{\\arraystretch}{1.3}" << endl
	  << "\\begin{tabular}{l *{" << samples.size()+1 << "}{c}}" << endl 
	  << "\\hline\\hline" << endl << "Uncertainty";
    for(unsigned int s=0 ; s<samples.size() ; ++s) {
      latex << " & " << samples[s].getNameLaTeX();
    }
    latex << " & Signal \\\\" << endl << "\\hline\\hline" << endl;

    for(unsigned int sy=0 ; sy<m_pSyste->getSize() ; ++sy) {
      const string systeName=m_pSyste->getName(sy);
      map<string,string>::const_iterator it=systeNames.find(systeName);
      if (it!=systeNames.end()) latex << it->second;
      else latex << systeName;
      for(unsigned int s=0 ; s<samples.size() ; ++s) {
	latex << " & " << samples[s].getLaTeXSystFromId(sy,precision);
      }
      latex << " & " << signal.getLaTeXSystFromId(sy,precision)
	    << " \\\\" << endl << "\\hline" << endl;
    }

    latex << "\\hline" << endl << "Total";
    for(unsigned int s=0 ; s<samples.size() ; ++s) {
      latex << " & " << samples[s].getLaTeXTotalSyst(precision);
    }
    latex << " & " << signal.getLaTeXTotalSyst(precision)
	  << " \\\\" << endl;

    latex << "\\hline" << endl << "\\end{tabular}" << endl 
	  << "\\end{center}\\end{table} " << endl << endl;
  }
}

TH1 *OpTHyLiC::getSystGaussDistr() const
{
  return m_pSyste->getDistr();
}

TH1 *OpTHyLiC::getHisto(const int i) const
{
  if (i<nbHistos) {
    TH1 *pH=m_pHs[i];
    pH->SetLineWidth(2);
    if (i==hLLRb) pH->SetLineColor(kBlue);
    else if (i==hLLRsb) pH->SetLineColor(kRed);
    return pH;
  }
  return 0;
}


TH1 *OpTHyLiC::getHistoLLRsb() const
{
  return getHisto(hLLRsb);
}

TH1 *OpTHyLiC::getHistoLLRb() const
{
  return getHisto(hLLRb);
}

void OpTHyLiC::setMuVsObs(const Observed &obs,const double mu)
{
  ++m_nbMu;
  m_sumMu+=mu;

  for(unsigned int c=0 ; c<m_muObsInterpol.size() && c<m_pChannels.size() ; ++c) {
    Observed obs1=obs;
    obs1.set(c,-1);
    map<Observed,MuVsObs> &muObs1=m_muObsInterpol[c];
    map<Observed,MuVsObs>::iterator it=muObs1.find(obs1);
    if (it!=muObs1.end()) it->second.add(m_pChannels[c]->getYieldData(),mu);
    else {
      MuVsObs muObs;
      muObs.add(m_pChannels[c]->getYieldData(),mu);
      muObs1[obs1]=muObs;
    }
  }  
}

void OpTHyLiC::createYieldTable(const int nbExp,ostream &latex,const int precision) const
{
  latex << "\\begin{table}\\begin{center}" << endl;
  if (0==nbExp) {
    latex << "\\caption{Observed yields and nominal expected yields. For each nominal expected yield, the first quoted uncertainty represent the statistical uncertainty, while the second is an approximation of the total systematic uncertainty, without taking into account the correlations between them.}" << endl;
  } else {
    latex << "\\caption{Observed yields and median expected yields. For each median expected yield, derived with pseudo-experiments, the quoted uncertainty is a combination of the statistical and systematic uncertainties, taking into account the correlations between systematics.}" << endl;
  }
  latex << "\\renewcommand{\\arraystretch}{1.3}" << endl
	<< "\\begin{tabular}{l *{" << m_pChannels.size() << "}{c}}" << endl 
	<< "\\hline\\hline" << endl << "Sample";
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    latex << " & " << m_pChannels[c]->getNameLaTeX();
    if (nbExp>0) m_pChannels[c]->generateDistrYield(nbExp);
  }
  latex << " \\\\" << endl << "\\hline\\hline" << endl;

  deque<string> bgList;
  set<string> bgSet;
  if (!m_pChannels.empty()) {
    const deque<Sample> &samples0=m_pChannels[0]->getBkgSamples();
    for(unsigned int s=0 ; s<samples0.size() ; ++s) {
      const string &bgName=samples0[s].getNameLaTeX();
      bgList.push_back(bgName);
      bgSet.insert(bgName);
    }
    for(unsigned int c=1 ; c<m_pChannels.size() ; ++c) {
      const deque<Sample> &samples=m_pChannels[c]->getBkgSamples();
      for(unsigned int s=0 ; s<samples.size() ; ++s) {
	const string &bgName=samples[s].getNameLaTeX();
	if (bgSet.find(bgName)==bgSet.end()) {
	  bgList.push_back(bgName);
	  bgSet.insert(bgName);
	}
      }  
    }
  }

  for(unsigned int s=0 ; s<bgList.size() ; ++s) {
    const string name=bgList[s];
    latex << name;

    for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
      latex << " & ";
      bool found=false;
      const deque<Sample> &samples=m_pChannels[c]->getBkgSamples();
      for(unsigned int s2=0 ; s2<samples.size() ; ++s2) {
	if (samples[s2].getNameLaTeX()==name) {
	  if (0==nbExp) latex << samples[s2].getYield().getLaTeX(precision);
	  else latex << samples[s2].getGeneratedYield().getLaTeX(precision,false);
	  found=true;
	  break;
	}
      }
      if (!found) latex << " --- ";
    }
    latex << " \\\\" << endl << "\\hline" << endl;
  }

  latex << "\\hline" << endl << "Total bkg.";
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    latex << " & ";
    if (0==nbExp) latex << m_pChannels[c]->getYieldBkg().getLaTeX(precision);
    else latex << m_pChannels[c]->getGeneratedYieldBkg().getLaTeX(precision,false);
  }
  latex << " \\\\" << endl;

  latex << "\\hline\\hline" << endl << "Data";
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    latex << " & " << m_pChannels[c]->getYieldData();
  }
  latex << " \\\\" << endl;

  latex << "\\hline\\hline" << endl << "Signal";
  for(unsigned int c=0 ; c<m_pChannels.size() ; ++c) {
    latex << " & ";
    if (0==nbExp) latex << m_pChannels[c]->getSigSample().getYield().getLaTeX(precision);
    else latex << m_pChannels[c]->getSigSample().getGeneratedYield().getLaTeX(precision,false);
  }
  latex << " \\\\" << endl;

  latex << "\\hline" << endl << "\\end{tabular}" << endl 
	<< "\\end{center}\\end{table} " << endl << endl;
}

