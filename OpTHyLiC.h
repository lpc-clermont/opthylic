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

#ifndef OPTHYLIC_H
#define OPTHYLIC_H

#include "OTHObserved.h"
#include "OTHChannel.h"

class OpTHyLiC: public OTH::Base {

public:

  OpTHyLiC(const OTH::SystType systInterpExtrapStyle,
	   const OTH::StatType statSampling,
	   const int RandomEngineType=OTH::TR3,
	   const int seed=0,
	   const OTH::CombType systCombinationType=OTH::CombAutomatic);

  virtual ~OpTHyLiC();

  // make multiple input files from single input with shapes
  void makeInputsFromShapes(const std::string &channelName,const std::string &fileName,const bool removeFiles=true);

  // set confidence level of computed limits
  virtual void setConfLevel(const double cl);

  // setting of samples yields and uncertainties
  unsigned int addChannel(const std::string &name);
  unsigned int addChannel(const std::string &name,const std::string &fileName,const bool removeFiles=true);
  
  // get pointer to specified channel, using its index
  OTH::Channel* getChannel(const unsigned int iChannel);
  // get pointer to specified channel, using its name
  OTH::Channel* getChannel(const std::string name);

  // setting of signal strength to be used for all computations
  virtual void setSigStrength(const double mu);

  // computation of the LLR value
  // combining all channels for the observed events in data
  double computeLLRdata() const;

  // generation of nbExp pseudo-experiments to compute the LLR distributions
  // must be called before trying to compute any CLs or p-value
  virtual void generateDistrLLR(const int nbExp);

  // computation of the p-value
  // the LLR distributions must have been generated before
  virtual double pValueData() const;

  // computation of the CLs value
  // combining all channels for the observed events in data
  // the LLR distributions must have been generated before
  double computeCLsData() const;

  // generation of nbExp pseudo-experiments with the signal strength mu
  // computation of the CLs for the observed events in data
  // (calls setSigStrength, generateDistrLLR and computeCLsData)
  virtual double generateForCLs(const double mu,const int nbExp,const int type);

  // methods called for observed and expected (median, -+1 sigma, +-2 sigma) limit computation
  virtual double sigStrengthExclusion(const OTH::LimitType type,const int nbExp,double &cls,
				      const double muHint=1,const OTH::MethType method=OTH::MethDichotomy);

  // computation of the distribution of the observed signal strengths if no signal exists
  // computation of the quantiles of this distribution (returns the median)
  // combining all channels
  double expectedSigStrengthExclusion(const int nbMu,const int nbExp);

  // print samples
  void printSamples() const;

  // print last systematics
  void printSyst() const;

  // print mu vs observation
  void printMuVsObs() const;

  // create LaTeX yield table
  void createInputYieldTable(std::ostream &latex,const int precision=-1) const;
  void createGeneratedYieldTable(std::ostream &latex,const int precision=-1,const int nbExp=1000000) const;

  // create LaTeX systematics tables
  // fileName contains the translation of systematics to LaTeX names
  void createSysteTables(std::ostream &latex,const std::string fileName,const int precision=2) const;

  // get systematic uncertainties base distribution
  TH1 *getSystGaussDistr() const;
  
  // histograms
  enum {hLLRb,hLLRsb,nbHistos};
  TH1 *getHisto(const int i) const;
  virtual TH1 *getHistoLLRsb() const;
  virtual TH1 *getHistoLLRb() const;

private:
  OpTHyLiC();
  OpTHyLiC(const OpTHyLiC&);
  OpTHyLiC &operator=(const OpTHyLiC&);

  bool isShape(const std::string &fileName) const;
  void setMuVsObs(const OTH::Observed &obs,const double mu);
  void createYieldTable(const int nbExp,std::ostream &latex,const int precision) const;

  OTH::RdmGenerator *m_pRdmGen; // random number generator
  OTH::Systematics *m_pSyste; // list of systematic uncertainties
  OTH::PdfGenerator *m_pStatSampling; // sampling method for stat uncertainty
  std::deque<OTH::Channel*> m_pChannels; // channels
  double m_sigStrength; // signal strength (scale factor of signal)
  double m_sumMu; // to compute average mu
  int m_nbMu; // to compute average mu

  // distributions
  std::vector<TH1*> m_pHs; // main histos
  std::map<OTH::Observed,double> m_muObs; // values of mu_95 for given observed events
  std::vector< std::map<OTH::Observed,OTH::MuVsObs> > m_muObsInterpol; // for interpolation
};
#endif // OPTHYLIC_H
