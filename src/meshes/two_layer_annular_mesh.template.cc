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

#ifndef OOMPH_ANNULAR_MESH_WITH_PML_TEMPLATE_CC
#define OOMPH_ANNULAR_MESH_WITH_PML_TEMPLATE_CC

#include "two_layer_annular_mesh.template.h"

namespace oomph
{

///Constructor
template<class ELEMENT>
TwoLayerAnnularMesh<ELEMENT>::
  TwoLayerAnnularMesh(const unsigned &ntheta,
                     const unsigned &nr_bulk,
                     const
                     unsigned &nr_pml,
                     const double &r_inner,
                     const double &r_middle,
                     const double &r_outer,
                     const double& mesh_perturbation,
                     TimeStepper* time_stepper_pt)
  : RectangularQuadMesh<ELEMENT>(ntheta, nr_bulk + nr_pml, 0.0, 1.0, 0.0,
    (nr_bulk + nr_pml)/double(nr_bulk), true, false, time_stepper_pt),
  Mesh_perturbation(mesh_perturbation)
{
  // Mesh can only be built with 2D Qelements.
  MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);
  
  // Build is delayed so that x_spacing_function can be overwritten
  this->build_mesh();
  
  // // Build the underlying Quadmesh
  // unsigned nr = nr_bulk + nr_pml;
  // // Thickness is chosen such that 0.0 is the at the inner and 1.0 is at the middle
  // double thickness = (r_outer - r_inner)/(r_middle - r_inner);
  // RectangularQuadMesh<ELEMENT>(ntheta,nr_bulk + nr_pml,1.0,(r_outer - r_inner)/(r_middle - r_inner),true,time_stepper_pt);
  
  // The wrap into annular shape will take the middle point to r_middle, 
  double stretch_factor = ((double) nr_bulk/ (double) nr_pml)*(r_outer-r_middle)/(r_middle-r_inner);
  double pml_start = 1.0;
  stretch_pml_region(pml_start,stretch_factor);
  
  // Wrap mesh into annular shape
  double phi=0.0;
  double azimuthal_fraction = 1.0;
  double h = r_middle - r_inner;
  wrap_into_annular_shape(r_inner,h,azimuthal_fraction,phi);
}

template<class ELEMENT>
double TwoLayerAnnularMesh<ELEMENT>::x_spacing_function(unsigned xelement,
                                                       unsigned xnode,
                                                       unsigned yelement,
                                                       unsigned ynode)
{
  double x = ((double) xelement + ((double) xnode)/(this->Np-1.0))/double(this->Nx);
  return x + Mesh_perturbation*sin(MathematicalConstants::Pi*x);
}

/// Wrap mesh into annular shape
template<class ELEMENT>
void TwoLayerAnnularMesh<ELEMENT>::wrap_into_annular_shape(const double& a,
                                                          const double& h,
                                                          const double& azimuthal_fraction,
                                                          const double& phi)
{
  //Create the hole
  Ellipse ellipse(a,a);

  //Set all the initial positions of the nodes
  Vector<double> xi(1);
  Vector<double> base(2);
  Vector<double> N(2);
  const unsigned n_node = this->nnode();
  for(unsigned n=0; n<n_node; n++)
  {
    // Pointer to node
    Node* nod_pt = this->node_pt(n);

    // Get the angle of the node -- rotate such that jump in angle
    // appears at periodic boundary. Shrink domain slightly
    // to keep angle unique
    xi[0] = (1.0-1.0e-10)*(-azimuthal_fraction*nod_pt->x(0))*
            2.0*MathematicalConstants::Pi+
            MathematicalConstants::Pi-
            (1.0-azimuthal_fraction)*2.0*MathematicalConstants::Pi;

    // Rotate
    xi[0]+=phi;

    //Get the node's fraction in the radial direction
    double w = nod_pt->x(1);

    //Get the position on the ellipse base
    ellipse.position(xi,base);

    //Get the unit normal, if it were a circle , by normalising the base
    double norm = sqrt(base[0]*base[0] + base[1]*base[1]);
    N[0] = base[0]/norm;
    N[1] = base[1]/norm;

    //Set circular film from the ellipse
    nod_pt->x(0) = base[0] + w*(h+a-norm)*N[0];
    nod_pt->x(1) = base[1] + w*(h+a-norm)*N[1];

    // Set boundary coordinates
    Vector<double> xi_bound(1);

    // Polar angle for boundary coordinate on boundary 0
    if(nod_pt->is_on_boundary(0))
    {
      xi_bound[0]=atan2(nod_pt->x(1),nod_pt->x(0));
      nod_pt->set_coordinates_on_boundary(0,xi_bound);
    }

    // Radius for boundary coordinate on boundary 1
    if(nod_pt->is_on_boundary(1))
    {
      xi_bound[0]=sqrt(pow(nod_pt->x(0),2)+pow(nod_pt->x(1),2));
      nod_pt->set_coordinates_on_boundary(1,xi_bound);
    }

    // Polar angle for boundary coordinate on boundary 2
    if(nod_pt->is_on_boundary(2))
    {
      xi_bound[0]=atan2(nod_pt->x(1),nod_pt->x(0));
      nod_pt->set_coordinates_on_boundary(2,xi_bound);
    }

    // Radius for boundary coordinate on boundary 3
    if(nod_pt->is_on_boundary(3))
    {
      xi_bound[0]=sqrt(pow(nod_pt->x(0),2)+pow(nod_pt->x(1),2));
      nod_pt->set_coordinates_on_boundary(3,xi_bound);
    }
  }

  this->Boundary_coordinate_exists[0]=true;
  this->Boundary_coordinate_exists[1]=true;
  this->Boundary_coordinate_exists[2]=true;
  this->Boundary_coordinate_exists[3]=true;
}

template<class ELEMENT>
void TwoLayerAnnularMesh<ELEMENT>::stretch_pml_region(
  const double& pml_start,
  const double& stretch_factor)
{
  
  // Add new inner PML boundary
  unsigned nboundary = this->nboundary();
  this->set_nboundary(nboundary + 1);

  // Loop over all th enodes in the mesh
  const unsigned n_node = this->nnode();
  for(unsigned n=0; n<n_node; n++)
  {
    // Pointer to node
    Node* nod_pt = this->node_pt(n);
    
    // If node lies on the inner PML boundary, add it to the new boundary
    if(std::fabs(nod_pt->x(1) - pml_start) < 1e-10)
    {
      // oomph_info << "Adding node at (" << nod_pt->x(0) << "," << nod_pt->x(1) << ") to middle boundary." << std::endl;
      // oomph_info << "pml_start = " << pml_start << std::endl;
      // oomph_info << "std::fabs(nod_pt->x(1) - pml_start) = " << std::fabs(nod_pt->x(1) - pml_start) << std::endl;
      this->convert_to_boundary_node(nod_pt); 
      this->add_boundary_node(nboundary,nod_pt);
    }
    
    // If it is in PML region of the initial rectangular mesh, stretch it
    if(nod_pt->x(1) > pml_start)
    {
      nod_pt->x(1) = (nod_pt->x(1) - pml_start)*stretch_factor + pml_start;
    }
  }
  
  // We added a new boundary and a load of nodes to said boundary, so we need
  // to refresh the lookup tables for boundary elements too
  this->setup_boundary_element_info();
}

}

#endif
