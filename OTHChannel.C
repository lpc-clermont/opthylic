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

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
using namespace std;

#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"

#include "OTHSystematics.h"
#include "OTHPdfGenerator.h"

#include "OTHChannel.h"
using namespace OTH;

Channel::Channel(const string &name,Systematics &syste,PdfGenerator &statSampling) :
  Base(),
  m_name(name),
  m_nameLaTeX(name),
  m_syste(syste),
  m_statSampling(statSampling),
  m_bgSamples(),
  m_sigSample(),
  m_sigStrength(1),
  m_yieldData(0),
  m_yieldSaved(0),
  m_yieldBg(0),
  m_yieldSB(0),
  m_cacheLog(0),
  m_cacheYieldS(0),
  m_cached(false),
  m_pHs(nbHistos,0),
  m_hNames(nbHistos,name),
  m_muObs(),
  m_muVsObs()
{
  m_hNames[hDistrBg]+="_Nb";
  m_hNames[hDistrSB]+="_Nsb";
  m_hNames[hLLRb]+="_LLRb";
  m_hNames[hLLRsb]+="_LLRsb";
  m_hNames[hYieldBg]+="_YieldBg";
}

Channel::~Channel()
{
  // do not delete histograms, belong to ROOT
}

YieldWithUncert Channel::getYieldBkg() const
{
  YieldWithUncert total;
  for(unsigned int s=0 ; s<m_bgSamples.size() ; ++s) {
    total+=m_bgSamples[s].getYield();
  }
  return total;
}

void Channel::addSamples(const std::string &fileName)
{
  ifstream in(fileName.c_str());
  if (!in) {
    cerr << "ERROR ! Unable to open file '" << fileName << "' !" << endl;
    return;
  }

  enum {None,Background,Signal};
  int type=None;
  unsigned int index=0;
  bool signalSet=false;
  for(; !in.eof() ;) {
    string line;
    if (!getline(in,line)) break;
    if (!line.empty() && line[0]!='#') {

      istringstream istr(line);
      string keyword,name;
      double val1,val2=0;
      istr >> keyword >> name >> val1 >> val2;

      if ("+bg"==keyword) {
	type=Background;
	index=addBkgSample(name,val1,val2);

      } else if ("+sig"==keyword) {
	type=Signal;
	setSigSample(name,val1,val2);
	if (signalSet) {
	  cerr << "WARNING ! Signal already set while setting new signal '" << name << "' !" << endl;
	}
	signalSet=true;

      } else if (".syst"==keyword) {
	if (None==type) {
	  cerr << "SYNTAX ERROR ! Trying to define a systematic uncertainty '" << name 
	       << "' before defining any sample !" << endl;
	} else if (Background==type) addBkgSystematics(index,name,val1,val2);
	else if (Signal==type) addSigSystematics(name,val1,val2);
	else cerr << "WHAT THE HELL !" << endl;

      } else if ("+data"==keyword) {
	istringstream istr2(name);
	double value=0;
	istr2 >> value;
	m_yieldData=static_cast<int>(value+0.5);

      } else if (".nameLaTeX"==keyword) {
	const string nameLaTeX=line.substr(line.find_first_of(' ')+1);
	if (None==type) {
	  cerr << "SYNTAX ERROR ! Trying to define a LaTeX name before defining any sample !" << endl;
	} else if (Background==type) m_bgSamples[index].setNameLaTeX(nameLaTeX);
	else if (Signal==type) m_sigSample.setNameLaTeX(nameLaTeX);
	else cerr << "WHAT THE HELL !" << endl;

      } else if ("+nameLaTeX"==keyword) {
	m_nameLaTeX=line.substr(line.find_first_of(' ')+1);

      } else {
	cerr << "Unknown line '" << line << "' !" << endl;
      }

    }
  }

  in.close();
}

unsigned int Channel::addBkgSample(const string &name,const double nominal,const double stat)
{
  string sName=m_name+"_"+name;
  // add new background sample
  m_bgSamples.push_back(Sample(sName,name,nominal,stat));

  // increment expected yields
  m_yieldBg+=nominal;
  m_yieldSB+=nominal;
  // reset cache for LLR
  m_cached=false;

  return m_bgSamples.size()-1;
}

void Channel::addBkgSystematics(const unsigned int iSample,
				const string &systName,
				const double up,const double down)
{
  // add systematic uncertainty for given background
  if (iSample<m_bgSamples.size()) {
    const unsigned int id=m_syste.add(systName);
    m_bgSamples[iSample].addSyst(systName,id,down,up);
  }
}

void Channel::setSigSample(const string &name,const double nominal,const double stat)
{
  string sName=m_name+"_"+name;
  // set signal sample
  m_sigSample=Sample(sName,name,nominal,stat);

  // recompute expected yield
  m_yieldSB=m_yieldBg+nominal*m_sigStrength;
  // reset cache for LLR
  m_cached=false;
}

void Channel::addSigSystematics(const string &systName,
				const double up,const double down)
{
  // add systematic uncertainty for signal
  const unsigned int id=m_syste.add(systName);
  m_sigSample.addSyst(systName,id,down,up);
}

void Channel::setSigStrength(const double mu)
{
  // set signal strength
  m_sigStrength=mu;

  // recompute expected yield
  m_yieldSB=m_yieldBg+m_sigSample.getNominal()*m_sigStrength;
  // reset cache for LLR
  m_cached=false;
}

double Channel::computeLLR(const int obs) const
{
  if (obs<0 || m_yieldSB<=0 || m_yieldBg<=0) {
    cerr << "OTHChannel Error ! (channel name=" << m_name << ") ";
    if (obs<0) cerr << "obs=" << obs << " ! ";
    if (m_yieldSB<=0) cerr << "m_yieldSB=" << m_yieldSB << " ! ";
    if (m_yieldBg<=0) cerr << "m_yieldBg=" << m_yieldBg << " ! ";
    throw runtime_error("Impossible to compute LLR !");
  }
  const double dobs=static_cast<double>(obs);
  if (!m_cached) {
    m_cacheLog=TMath::Log(m_yieldSB/m_yieldBg);
    m_cacheYieldS=m_yieldSB-m_yieldBg;
    m_cached=true;
  }
  return 2*(m_cacheYieldS-dobs*m_cacheLog);
}

void Channel::initDistrLLR(double &llrMin,double &llrMax)
{
  // resetting
  for(unsigned int h=0 ; h<m_pHs.size() ; ++h) {
    if (m_pHs[h]) {
      delete m_pHs[h];
      m_pHs[h]=0;
    }
  }
  m_cached=false;

  // creation of histograms to store distributions
  const int maxEvt=static_cast<int>(5*m_yieldSB)+1;
  m_pHs[hDistrBg]=new TH1F(m_hNames[hDistrBg].c_str(),";Events;Probability",maxEvt,-0.5,static_cast<float>(maxEvt)-0.5);
  m_pHs[hDistrSB]=new TH1F(m_hNames[hDistrSB].c_str(),";Events;Probability",maxEvt,-0.5,static_cast<float>(maxEvt)-0.5);

  llrMin=computeLLR(maxEvt);
  llrMax=computeLLR(0)*2;
  m_pHs[hLLRb]=new TH1F(m_hNames[hLLRb].c_str(),";LLR;Probability",1000,llrMin,llrMax);
  m_pHs[hLLRsb]=new TH1F(m_hNames[hLLRsb].c_str(),";LLR;Probability",1000,llrMin,llrMax);
}

void Channel::generateSinglePseudoExp(double &llrB,double &llrSB)
{
  double expected=0;

  // compute expected background contribution and Poisson varied
  const int expBg=generateSinglePseudoExpBg(expected);
  m_pHs[hDistrBg]->Fill(expBg);

  // compute test-statistic for b
  llrB=computeLLR(expBg);
  m_pHs[hLLRb]->Fill(llrB);

  // add signal
  expected+=generateSingleSample(m_sigSample,m_sigStrength);
  // vary expectation with Poisson statistics
  const int expSB=m_statSampling.poisson(expected);
  m_pHs[hDistrSB]->Fill(expSB);

  // compute test-statistic for s+b
  llrSB=computeLLR(expSB);
  m_pHs[hLLRsb]->Fill(llrSB);
}

void Channel::generateDistrLLR(const int nbExp)
{
  double dummy1,dummy2;
  initDistrLLR(dummy1,dummy2);

  if (nbExp<1) return;

  // loop on all pseudo-experiments
  for(int i=0 ; i<nbExp ; ++i) {
    // systematic uncertainties variations
    m_syste.variate();

    // single exp generation
    generateSinglePseudoExp(dummy1,dummy2);
  }

  endDistrLLR(nbExp);
}

void Channel::endDistrLLR(const int nbExp)
{
  // normalization of distributions
  m_pHs[hDistrBg]->Scale(1/static_cast<float>(nbExp));
  m_pHs[hDistrSB]->Scale(1/static_cast<float>(nbExp));
  m_pHs[hLLRb]->Scale(1/static_cast<float>(nbExp));
  m_pHs[hLLRsb]->Scale(1/static_cast<float>(nbExp));
}

double Channel::pValue(const int obs) const
{
  if (!m_pHs[hDistrBg]) {
    cout << "ERROR ! No background event distribution found, use generateDistrLLR first..." << endl;
    return 0;
  }
  return m_pHs[hDistrBg]->Integral(m_pHs[hDistrBg]->FindBin(obs),-1);
}

double Channel::pValueData() const
{
  return pValue(m_yieldData);
}

double Channel::computeCLs(const int obs) const
{
  if (!m_pHs[hLLRsb] || !m_pHs[hLLRb]) {
    cout << "ERROR ! No LLR distribution found, use generateDistrLLR first..." << endl;
    return -1;
  }

  return Algorithms::computeCLs(m_pHs[hLLRsb],m_pHs[hLLRb],computeLLR(obs));
}

int Channel::findObsExclusion() const
{
  for(int obs=0 ; obs<static_cast<int>(1+100*m_yieldSB) ; ++obs) {
    const double cls=computeCLs(obs);
    if (cls>(1-m_confLevel)) return obs-1;
  }
  return -1;
}

double Channel::generateForCLs(const double mu,const int nbExp,const int type)
{
  setSigStrength(mu);
  generateDistrLLR(nbExp);
  if(LimObserved==type) {
    return computeCLs(m_yieldData);
  }
  else if(type>=LimExpectedP2sig && type<=LimExpectedM2sig) {
    return Algorithms::getCLsFromLLR(type,m_pHs[hLLRsb],m_pHs[hLLRb]);
  }
  else {
    throw runtime_error("Unknown limit type !");
  }
}

void Channel::generateDistrYield(const int nbExp)
{
  for(unsigned int s=0 ; s<m_bgSamples.size() ; ++s) {
    m_bgSamples[s].createYieldHisto();
  }
  m_sigSample.createYieldHisto();
  if (m_pHs[hYieldBg]) delete m_pHs[hYieldBg];
  m_pHs[hYieldBg]=new TH1F(m_hNames[hYieldBg].c_str(),";Events;Probability",1000,0,5*m_yieldBg);

  if (nbExp<1) return;

  // loop on all pseudo-experiments
  for(int i=0 ; i<nbExp ; ++i) {
    // systematic uncertainties variations
    m_syste.variate();

    double expected=0;
    // sum up all background samples
    for(unsigned int s=0 ; s<m_bgSamples.size() ; ++s) {
      const double yield=generateSingleSample(m_bgSamples[s]);
      m_bgSamples[s].fillYieldHisto(yield);
      expected+=yield;
    }
    m_pHs[hYieldBg]->Fill(expected);

    m_sigSample.fillYieldHisto(generateSingleSample(m_sigSample));
  }

  for(unsigned int s=0 ; s<m_bgSamples.size() ; ++s) {
    m_bgSamples[s].getYieldHisto()->Scale(1/static_cast<float>(nbExp));
  }
  m_sigSample.getYieldHisto()->Scale(1/static_cast<float>(nbExp));
  m_pHs[hYieldBg]->Scale(1/static_cast<float>(nbExp));
}

YieldWithUncert Channel::getGeneratedYieldBkg() const
{
  YieldWithUncert yield;
  if (m_pHs[hYieldBg]) {
    vector<double> quant=Algorithms::getQuantiles(m_pHs[hYieldBg],false);
    yield.setYield(quant[2]);
    yield.setSystHigh(quant[3]-quant[2]);
    yield.setSystLow(quant[1]-quant[2]);
  }
  return yield;
}

double Channel::sigStrengthExclusion(const LimitType type,const int nbExp,double &cls,
				     const double muHint,const MethType method)
{
  double mu=0.5,muStep=3;
  if (muHint!=1) mu=muHint/2;
  else if(LimObserved==type) {
    cout << "--------- Searching mu for obs=" << m_yieldData << " -----------" << endl;

    // if already two mus, interpolate to find first value of mu
    mu=m_muVsObs.interpolateMu(m_yieldData);
    if (mu<0) mu=1;
    else muStep=1.2;
  }

  return Algorithms::sigStrengthExclusion(*this,mu,muStep,nbExp,type,cls,m_confLevel,method!=MethDichotomy);
}

double Channel::expectedSigStrengthExclusion(const int nbMu,const int nbExp)
{
  saveYieldData();

  // resetting list of CLs values
  m_muObs.clear();
  m_muVsObs.reset();
  double cls=0;
  m_yieldData=static_cast<int>(m_yieldBg);
  m_muObs[m_yieldData]=sigStrengthExclusion(LimObserved,nbExp,cls);
  m_muVsObs.add(m_yieldData,m_muObs[m_yieldData]);

  // resetting histograms
  if (m_pExpMu) {
    delete m_pExpMu;
    m_pExpMu=0;
  }
  string hName=m_name+"_mu";
  m_pExpMu=new TH1F(hName.c_str(),";Expected signal strength;Probability",1000,0,m_muObs[m_yieldData]*10);
  if (m_pCLs) {
    delete m_pCLs;
    m_pCLs=0;
  }
  hName=m_name+"_CLs";
  m_pCLs=new TH1F(hName.c_str(),";CL_{s};Entries",1000,(1-m_confLevel)*0.8,(1-m_confLevel)*1.2);
  m_pCLs->Fill(cls);

  // loop on background only pseudo-experiments
  for(int i=0 ; i<nbMu ; ++i) {
    // systematic uncertainties variations
    m_syste.variate();

    // generate pseudo-experiments
    double initial=0;
    const int obs=generateSinglePseudoExpBg(initial);

    // search for mu_95 in map, otherwise compute it
    map<int,double>::const_iterator it=m_muObs.find(obs);
    double mu=0;
    if (it!=m_muObs.end()) mu=it->second;
    else {
      m_yieldData=obs;
      mu=sigStrengthExclusion(LimObserved,nbExp,cls);
      m_muObs[obs]=mu;
      m_muVsObs.add(obs,mu);
      m_pCLs->Fill(cls);
    } 
    m_pExpMu->Fill(mu);
  }
  m_pExpMu->Scale(1/static_cast<float>(nbMu));

  // filling mu_95 vs obs graph
  if (m_pMuObs) {
    delete m_pMuObs;
    m_pMuObs=0;
  }
  m_pMuObs=new TGraph(m_muObs.size());
  m_pMuObs->SetMarkerStyle(kCircle);
  int i=0;
  for(map<int,double>::const_iterator it=m_muObs.begin() ; it!=m_muObs.end() ; ++it,++i) {
    m_pMuObs->SetPoint(i,it->first,it->second);
  }

  restoreYieldData();
  vector<double> muQ=Algorithms::getQuantiles(m_pExpMu);
  return muQ[2];
}

int Channel::generateSinglePseudoData(const double mu)
{
  double expected=0;
  generateSinglePseudoExpBg(expected);
  if(mu!=0) {
    expected+=generateSingleSample(m_sigSample,mu);
  }
  m_yieldData=m_statSampling.poisson(expected);
  return m_yieldData;
}

void Channel::printSamples() const
{
  cout << "======= List of samples =============" << endl
       << "-> channel '" << m_name << "'" << endl
       << "--------- background ----------------" << endl;
  for(unsigned int b=0 ; b<m_bgSamples.size() ; ++b) m_bgSamples[b].print();
  cout << "---> Total expected background = " << m_yieldBg << endl;
  cout << "----------- signal ------------------" << endl;
  m_sigSample.print();
  cout << "-> signal strength: mu = " << m_sigStrength << endl
       << "---> Total expected background+mu*signal = " << m_yieldSB << endl
       << "---> Observed yield in data = " << m_yieldData << endl;
}

TH1 *Channel::getSigSystDistr(const std::string &systName) const
{
  return m_sigSample.getSystDistr(systName);
}

TH1 *Channel::getBkgSystDistr(const std::string &bkgName,const std::string &systName) const
{
  const string name=m_name+"_"+bkgName;
  for(unsigned int i=0 ; i<m_bgSamples.size() ; ++i) {
    const Sample &smpl=m_bgSamples[i];
    if (smpl.getName()==name) return smpl.getSystDistr(systName);
  }
  throw runtime_error("Unknown background sample name !");
}

TH1 *Channel::getBkgYieldDistr(const std::string &bkgName) const
{
  const string name=m_name+"_"+bkgName;
  for(unsigned int i=0 ; i<m_bgSamples.size() ; ++i) {
    const Sample &smpl=m_bgSamples[i];
    if (smpl.getName()==name) return smpl.getYieldHisto();
  }
  throw runtime_error("Unknown background sample name !");
}

TH1* Channel::getHisto(const int i) const
{
  if (i<nbHistos) {
    TH1 *pH=m_pHs[i];
    pH->SetLineWidth(2);
    if (i==hDistrBg || i==hLLRb || i==hYieldBg) pH->SetLineColor(kBlue);
    else if (i==hDistrSB || i==hLLRsb) pH->SetLineColor(kRed);
    return pH;
  }
  return 0;
}

TH1 *Channel::getHistoLLRsb() const
{
  return getHisto(hLLRsb);
}

TH1 *Channel::getHistoLLRb() const
{
  return getHisto(hLLRb);
}

//////////////////////////////////////////
// Private methods

int Channel::generateSinglePseudoExpBg(double &expected) const
{
  // sum up all background samples
  for(unsigned int s=0 ; s<m_bgSamples.size() ; ++s) {
    expected+=generateSingleSample(m_bgSamples[s]);
  }
  // vary expectation with Poisson statistics
  return m_statSampling.poisson(expected);
}

double Channel::generateSingleSample(const Sample &sample,const double mu) const
{
  // apply statistical uncertainty to sample
  double expSamp=0;
  if(sample.getStat()==0) expSamp=sample.getNominal()*mu;
  else expSamp=m_statSampling.draw(sample.getNominal()*mu,sample.getStat()*mu);
  
  // apply systematics
  double systScale=m_additiveSystComb?0:1;
  for(unsigned int i=0 ; i<sample.getSystSize() ; ++i) {
    const double var=m_syste.getScaleFactor(sample.getSystId(i),
					    sample.getSystLow(i),
					    sample.getSystHigh(i));
    if(m_additiveSystComb) systScale+=var-1;
    else systScale*=var;
    sample.fillSystDistr(i,var);
  }
  if(m_additiveSystComb) systScale=1+systScale;
  expSamp*=systScale;

  return expSamp;
}

