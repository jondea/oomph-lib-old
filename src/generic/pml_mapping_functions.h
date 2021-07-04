//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#ifndef OOMPH_PML_MAPPING_FUNCTIONS_HEADER
#define OOMPH_PML_MAPPING_FUNCTIONS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

#include "complex_matrices.h"
#include "oomph_utilities.h"

namespace oomph
{

//=======================================================================
/// Class to hold the mapping function (gamma) for the Pml which defines
/// how the coordinates are transformed in the Pml. This class holds
/// the one dimensional or uniaxial case which is the most common 
//=======================================================================
class UniaxialPMLMapping
{

public:

  /// Default constructor (empty)
  UniaxialPMLMapping(){};

  /// \short Pure virtual to return Pml mapping gamma, where gamma is the
  /// \f$d\tilde x / d x\f$ as  function of \f$\nu\f$ where \f$\nu = x - h\f$ where h is
  /// the vector from the origin to the start of the Pml
  virtual std::complex<double> dtransformed_nu_dnu(const double& nu,
    const double& delta,
    const double& wavenumber_squared,
    const double& alpha_shift=0.0) = 0;
  
  //
  virtual std::complex<double> transformed_nu(const double& nu,
    const double& delta,
    const double& wavenumber_squared,
    const double& alpha_shift=0.0) = 0;
 
};

//=======================================================================
/// A mapping function propsed by Bermudez et al, appears to be the best
/// for the Helmholtz equations and so this will be the default mapping 
/// (see definition of PmlHelmholtzEquations)
//=======================================================================
class BermudezPMLMapping : public UniaxialPMLMapping
{

public:

  /// Default constructor (empty)
  BermudezPMLMapping(){};

  /// \short Overwrite the pure Pml mapping coefficient function to return the
  /// coeffcients proposed by Bermudez et al
  std::complex<double> dtransformed_nu_dnu(const double& nu,
                                           const double& delta,
                                           const double& k,
                                           const double& alpha_shift=0.0)
  {
    return 1.0 + MathematicalConstants::I / k * (1.0/std::fabs(delta - nu));
  }
  
  //
  std::complex<double> transformed_nu(const double& nu,
                                      const double& delta,
                                      const double& k,
                                      const double& alpha_shift=0.0)
  {
    return nu - MathematicalConstants::I/k * log(1.0 - std::fabs(nu/delta));
  }

};

//=======================================================================
/// A mapping function proposed by Bermudez et al, similar to the one above
/// but is continuous across the inner Pml boundary
/// appears to be the best for TimeHarmonicLinearElasticity
/// and so this will be the default mapping 
//=======================================================================
class C1BermudezPMLMapping : public UniaxialPMLMapping
{

public:

  /// Default constructor (empty)
  C1BermudezPMLMapping(){};

  /// \short Overwrite the pure Pml mapping coefficient function to return the
  /// coeffcients proposed by Bermudez et al
  std::complex<double> dtransformed_nu_dnu(const double& nu,
                                           const double& delta,
                                           const double& k,
                                           const double& alpha_shift=0.0)
  {
    /// return \f$\gamma=1 + (i/k)(1/|outer_boundary - x|-1/|pml width|)\f$ 
    return 1.0 + MathematicalConstants::I / k
                  *( 1.0/std::fabs(delta-nu) - 1.0/std::fabs(delta));
  }

  //
  std::complex<double> transformed_nu(const double& nu,
                                      const double& delta,
                                      const double& k,
                                      const double& alpha_shift=0.0)
  {
    return nu - MathematicalConstants::I/k
                  *( log(1.0-std::fabs(nu/delta)) - nu/std::fabs(delta) );
  }

};

//=======================================================================
/// A slight modification to the mapping function propsed by Bermudez et al.
/// The transformation is independent of the scale (thickness) of the PML.
/// See Jonathan Deakin's PhD thesis (University of Manchester) for more
/// information.
//=======================================================================
class ScaleFreeBermudezPMLMapping : public UniaxialPMLMapping
{

public:

  /// Default constructor (empty)
  ScaleFreeBermudezPMLMapping(){};

  /// \short Overwrite the pure Pml mapping coefficient function to return the
  /// coeffcients proposed by Bermudez et al
  std::complex<double> dtransformed_nu_dnu(const double& nu,
                                           const double& delta,
                                           const double& k,
                                           const double& alpha_shift=0.0)
  {
    return MathematicalConstants::I / k * (1.0/std::fabs(delta - nu));
  }
  
  //
  std::complex<double> transformed_nu(const double& nu,
                                      const double& delta,
                                      const double& k,
                                      const double& alpha_shift=0.0)
  {
    return - MathematicalConstants::I/k * log(1.0 - std::fabs(nu/delta));
  }

};

//=======================================================================
/// Class to hold the mapping function (gamma) for the PML which defines
/// how the coordinates are transformed in the PML. This PML mapping
/// aims to transform the solution into a straight line, and requires a 
/// lot of information from the element. Returns a mapping Jacobian which
/// I haven't quite worked out yet
//=======================================================================
class TangentiallyVaryingConformalPMLMapping
{

public:

  /// Default constructor (empty)
  TangentiallyVaryingConformalPMLMapping(){};

  /// \short Pure virtual to return PML mapping Jacobian
  virtual void get_mapping_jacobian(
    const Vector<double>& x,
    const Vector<double>& x_inner,
    const Vector<double>& x_outer,
    const Vector<double>& dx_inner_dacross,
    const Vector<double>& dp_dacross,
    const double& k,
    std::complex<double>& tnu,
    std::complex<double>& dtnu_dnu,
    std::complex<double>& dtnu_dacross,
    const double& alpha=0.0
  ) = 0;
};

//=======================================================================
/// Class used to calculate the PML transformation and its Jacobian.
/// This specific variant is defined in a conformal geometry, and the
/// transformation can be discontinuous in the tangential direction
/// (across the PML)
//=======================================================================
class TangentiallyDiscontinuousConformalPMLMapping :
public TangentiallyVaryingConformalPMLMapping
{

public:

  /// Default constructor (empty)
  TangentiallyDiscontinuousConformalPMLMapping(){};

  /// \short Search along the line 0 to nun for pole (or closest point to)
  virtual void pole_line_search(
    const Vector<double>& x_inner,
    const Vector<double>& p,
    const double& k,
    double& nun,
    std::complex<double>& tnu,
    const double& alpha=0.0) = 0;

  /// Make a Newton step towards a pole (pole is when du/dtnu=0s)
  virtual bool newton_step_to_pole(
    const Vector<double>& x_inner,
    const Vector<double>& p,
    const Vector<double>& dx_inner_dacross,
    const Vector<double>& dp_dacross,
    double& nu,
    double& zeta,
    std::complex<double>& tnu,
    std::complex<double>& alpha,
    std::complex<double>& beta) = 0;

  /// Can this be taken out?
  /// Make a Newton step towards a pole (pole is when du/dtnu=0s)
  virtual void set_initial_guess(
    const Vector<double>& x_inner,
    const Vector<double>& p,
    const double& nun,
    const std::complex<double>& tnu) = 0;
};

}

#endif
