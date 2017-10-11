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

#ifndef OTH_TYPES_H
#define OTH_TYPES_H

namespace OTH {

  // Type of limit
  enum LimitType {LimExpectedP2sig, // expected +2 sigma
		  LimExpectedP1sig, // expected +1 sigma
		  LimExpectedMed, // expected median
		  LimExpectedM1sig, // expected -1 sigma
		  LimExpectedM2sig, // expected -2 sigma
		  LimObserved}; // observed
    
  // Significance type
  enum SignifType {SignifExpectedP2sig, // expected +2 sigma
		   SignifExpectedP1sig, // expected +1 sigma
		   SignifExpectedMed, // expected median
		   SignifExpectedM1sig, // expected -1 sigma
		   SignifExpectedM2sig, // expected -2 sigma
		   SignifObserved}; // observed

  // Type of engine for pseudo-random number generation
#if defined CPP11 || defined __CLING__
  enum {TR3, // Provided by TRandom3 class
	STD_minstd_rand, // Provided by the C++11 standard library
	STD_minstd_rand0, // Provided by the C++11 standard library
	STD_mt19937, // Provided by the C++11 standard library
	STD_mt19937_64, // Provided by the C++11 standard library
	STD_ranlux24_base, // Provided by the C++11 standard library
	STD_ranlux48_base, // Provided by the C++11 standard library
	STD_ranlux24, // Provided by the C++11 standard library
	STD_ranlux48, // Provided by the C++11 standard library
	STD_knuth_b}; // Provided by the C++11 standard library
#else
  enum {TR3}; // Provided by TRandom3 class
#endif
    
  // Type of combination of systematics
  enum CombType {CombAdditive, // additive
		 CombMultiplicative, // multiplicative
		 CombAutomatic}; // depends on the type of systematics interp/extrapol functions

  // Type of systematics interplation/extrapolation
  enum SystType {SystMclimit, // a la McLimit
		 SystLinear, // linear
		 SystExpo, // exponential
		 SystPolyexpo}; // polynomial interpolation and exponential extrapolation

  // Types of statistical PDF generator
  enum StatType {StatNormal, // normal distribution (Gauss)
		 StatLogN, // log-normal
		 StatGammaHyper, // gamma with hyperbolic prior
		 StatGammaUni, // gamma with uniform prior
		 StatGammaJeffreys}; // gamma with Jeffreys prior

  // Type of method for CLs(mu) computation
  enum MethType {MethDichotomy, // using log-dichotomy method
		 MethExtrapol}; // using simple extrapolation
}

#endif // OTH_TYPES_H
