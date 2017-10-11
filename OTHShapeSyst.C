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

#include "OTHShapeSyst.h"
using namespace OTH;

ShapeSyst::ShapeSyst(const string& name, TH1* histoUp, TH1* histoDown) :
  m_name(name),
  m_histoUp(histoUp),
  m_histoDown(histoDown)
{
  if(m_histoUp->GetNbinsX() != m_histoDown->GetNbinsX()) {
    cerr << "ERROR ! number of bins in histoUp and histoDown is different (histoUp:" 
	 << m_histoUp->GetNbinsX() 
	 << " bins, histoDown: " 
	 << m_histoDown->GetNbinsX()
	 << " bins)" << endl;
    return;
  }
}

ShapeSyst::~ShapeSyst()
{
  delete m_histoUp;
  delete m_histoDown;
}

int ShapeSyst::getNbins() const 
{
  return m_histoUp->GetNbinsX();
}

float ShapeSyst::getBinContentUp(const int i) const 
{
  return m_histoUp->GetBinContent(i);
}

float ShapeSyst::getBinContentDown(const int i) const 
{
  return m_histoDown->GetBinContent(i);
}

void ShapeSyst::print() const
{
  cout << "  * name: " << m_name << endl;
  cout << "   -> histoUp content: " << endl;
  for(int i=0; i<=m_histoUp->GetNbinsX()+1; ++i) {
    cout << "     o bin " << i << ": " 
	 << m_histoUp->GetBinContent(i) 
	 << endl;
    }
  cout << "   -> histoDown content: " << endl;
  for(int i=0; i<=m_histoDown->GetNbinsX()+1; ++i) {
    cout << "     o bin " << i << ": " 
	 << m_histoDown->GetBinContent(i) 
	 << endl;
    }
}

