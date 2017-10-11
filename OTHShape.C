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

#include "TH1.h"

#include "OTHShape.h"
using namespace OTH;
using namespace std;

Shape::Shape(const string &name,TH1* histoNominal) :
  m_name(name),
  m_nameLaTeX(""),
  m_histoNominal(histoNominal),
  m_shapeSyst(),
  m_globalSyst()
{}

Shape::~Shape()
{
  delete m_histoNominal;
  for(unsigned int i=0; i<m_shapeSyst.size(); ++i) {
    delete m_shapeSyst[i];
  }
  for(unsigned int i=0; i<m_globalSyst.size(); ++i) {
    delete m_globalSyst[i];
  }
}

void Shape::addShapeSystematic(ShapeSyst* shapeSyst) 
{
  if(m_histoNominal->GetNbinsX() != shapeSyst->getNbins()) {
    cerr << "ERROR ! number of bins in histoNominal and varied histograms is different (histoNominal:" 
	 << m_histoNominal->GetNbinsX() 
	 << " bins, varied histos: " 
	 << shapeSyst->getNbins()
	 << " bins)" << endl;
    return;
  }
  m_shapeSyst.push_back(shapeSyst);
}

void Shape::addGlobalSystematic(SingleSyst* globalSyst) 
{
  m_globalSyst.push_back(globalSyst);
}

float Shape::getBinContent(const int i) const 
{
  return m_histoNominal->GetBinContent(i);
}

float Shape::getBinError(const int i) const 
{
  return m_histoNominal->GetBinError(i);
}

float Shape::getBinCenter(const int i) const 
{
  return m_histoNominal->GetBinCenter(i);
}

void Shape::print(const bool printStatUncert) const
{
  cout << "======= Shape =========" << endl;
  cout << " -> name: " << m_name << endl;
  cout << " -> Nominal histo content: " << endl; 
    for(int i=0; i<=m_histoNominal->GetNbinsX()+1; ++i) {
      cout << "     o bin " << i << ": " 
	   << m_histoNominal->GetBinContent(i) ;
      if(printStatUncert) 
	cout << " +- " << m_histoNominal->GetBinError(i) << endl;
      else
	cout << endl;
    }
    cout << " -> Varied histos (total number of systematics: " << m_shapeSyst.size() << "): " << endl;
    for(unsigned int i=0; i<m_shapeSyst.size(); ++i) {
       m_shapeSyst[i]->print();
    }  
}

void Shape::writeInputFile(ostream &out,const int iBin) const
{
  const float binContent=getBinContent(iBin);
  out << m_name << " " << binContent << " " << getBinError(iBin) << endl;
  if (m_nameLaTeX!="") {
    out << ".nameLaTeX " << m_nameLaTeX << endl;
  }
  if (binContent!=0) {
    for(unsigned int k=0; k<m_shapeSyst.size(); ++k) {
      const ShapeSyst* pSyst=m_shapeSyst[k];
      if (pSyst->getBinContentUp(iBin)!=binContent || pSyst->getBinContentDown(iBin)!=binContent) {
	out << ".syst " << pSyst->getName() 
	    << " " << (pSyst->getBinContentUp(iBin)-binContent)/binContent 
	    << " " << (pSyst->getBinContentDown(iBin)-binContent)/binContent 
	    << endl;
      }	
    }
  }
  for(unsigned int k=0; k<m_globalSyst.size(); ++k) {
    const SingleSyst* pSyst=m_globalSyst[k];
    if (pSyst->getLow()!=0 || pSyst->getHigh()!=0) {
      out << ".syst " << pSyst->getName() 
	  << " " << pSyst->getHigh()
	  << " " << pSyst->getLow()
	  << endl;
    }
  }
}


