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

#ifndef OOMPH_TWO_LAYER_ANNULAR_MESH_TEMPLATE_H
#define OOMPH_TWO_LAYER_ANNULAR_MESH_TEMPLATE_H

#include "rectangular_quadmesh.h"

namespace oomph
{


//===================================================================
/// 2D annular mesh with a unit circle in the middle and a layer
/// of thickness h surrounding it.
//==================================================================
template<class ELEMENT>
class TwoLayerAnnularMesh : public virtual RectangularQuadMesh<ELEMENT>
{

public:


  ///Constructor
  TwoLayerAnnularMesh(const unsigned &ntheta,
                      const unsigned &nr_bulk,
                      const
                      unsigned &nr_pml,
                      const double &r_inner,
                      const double &r_middle,
                      const double &r_outer,
                      const double& mesh_perturbation = 0.0,
                      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

  double x_spacing_function(unsigned xelement,
                            unsigned xnode,
                            unsigned yelement,
                            unsigned ynode);

private:

  double Mesh_perturbation;

  /// Wrap mesh into annular shape
  void wrap_into_annular_shape(const double& a,
                               const double& h,
                               const double& azimuthal_fraction,
                               const double& phi);
  // 
  void stretch_pml_region(const double& pml_start,
                          const double& stretch_factor);

};

}

#endif
