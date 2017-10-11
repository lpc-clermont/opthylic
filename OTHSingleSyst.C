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
using namespace std;

#include "TH1.h"

#include "OTHSingleSyst.h"
using namespace OTH;

SingleSyst::SingleSyst() :
  m_name(),
  m_id(0),
  m_low(0),
  m_high(0),
  m_pH(0)
{}

SingleSyst::SingleSyst(const string &name,const unsigned int id,
		       const double low,const double high) :
  m_name(name),
  m_id(id),
  m_low(low),
  m_high(high),
  m_pH(0)
{}

SingleSyst::~SingleSyst()
{
  // do not delete m_pH, belongs to ROOT
}

SingleSyst::SingleSyst(const SingleSyst &syst) :
  m_name(syst.m_name),
  m_id(syst.m_id),
  m_low(syst.m_low),
  m_high(syst.m_high),
  m_pH(0)
{}

SingleSyst &SingleSyst::operator=(const SingleSyst &syst)
{
  if (this!=&syst) {
    m_name=syst.m_name;
    m_id=syst.m_id;
    m_low=syst.m_low;
    m_high=syst.m_high;
    m_pH=0;
  }
  return *this;
}

void SingleSyst::createDistr(const string &smplName)
{
  double maxi=2;
  if (m_high>m_low && m_high>0) maxi=(1+m_high)*2;
  else if (m_low>m_high && m_low>0) maxi=(1+m_low)*2;

  string hName=smplName;
  hName+="_";
  hName+=m_name;
  string hTitle=hName;

  for(size_t it=hName.find(' ') ; it!=string::npos ; it=hName.find(' ')) hName[it]='_';
  for(size_t it=hTitle.find(';') ; it!=string::npos ; it=hTitle.find(';')) hTitle[it]=',';
  hTitle+=";Variation;Entries";

  if (0!=m_pH) delete m_pH;
  m_pH=new TH1I(hName.c_str(),hTitle.c_str(),100,0,maxi);
}

void SingleSyst::fillDistr(const double value) const
{
  if (m_pH) m_pH->Fill(value);
}

void SingleSyst::print() const
{
  cout << " -- syst '" << m_name << "' (" << m_id << "): "
       << m_high*100 << "%" << " " << m_low*100 << "%" << endl;
}

