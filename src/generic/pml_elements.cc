//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1282 $
//LIC//
//LIC// $LastChangedDate: 2017-01-16 08:27:53 +0000 (Mon, 16 Jan 2017) $
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

#include "pml_elements.h"
#include "shape.h"

namespace oomph
{

template<> 
void PMLElementBase<1>::
compute_laplace_matrix_and_det(const DenseComplexMatrix& J, 
                               DenseComplexMatrix& laplace_matrix, 
                               std::complex<double>& detJ)
{
  const unsigned DIM = 1;
#ifdef PARANOID
  // check matrix dimensions are compatable and the DIM has been implemented
  if ( J.nrow() != DIM && J.ncol() != DIM )
  {
     throw OomphLibError("Input matrix must be 1x1", OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
  }
#endif
  detJ = J(0,0);

  // resize and intialize result
  laplace_matrix.resize(DIM, DIM, 0.0);
  
  laplace_matrix(0,0) = 1.0/J(0,0);
}

template<> 
void PMLElementBase<1>::
compute_laplace_matrix_and_det(const DiagonalComplexMatrix& J, 
                               DiagonalComplexMatrix& laplace_matrix, 
                               std::complex<double>& detJ)
{
  const unsigned DIM = 1;
  
  detJ = J(0,0);  
  
  // resize and intialize result
  laplace_matrix.resize(DIM, DIM, std::complex<double>(0.0,0.0));

  laplace_matrix(0,0) = 1.0 / J(0,0);
}

template<> 
void PMLElementBase<1>::
compute_laplace_matrix_and_det(const DiagonalComplexMatrix& J, 
                               DenseComplexMatrix& laplace_matrix, 
                               std::complex<double>& detJ)
{
  const unsigned DIM = 1;

  detJ = J(0,0);  
  
  // resize and intialize result
  laplace_matrix.resize(DIM, DIM, std::complex<double>(0.0,0.0));

  laplace_matrix(0,0) = 1.0 / J(0,0);
}

template<> 
void PMLElementBase<1>::
compute_jacobian_inverse_and_det(const DiagonalComplexMatrix& J, 
                               DiagonalComplexMatrix& J_inv, 
                               std::complex<double>& J_det)
{
  const unsigned DIM = 1;
  
  J_det = J(0,0);  
  
  // resize and intialize result
  J_inv.resize(DIM, DIM, std::complex<double>(0.0,0.0));

  J_inv(0,0) = 1.0 / J(0,0);
}

template<> 
void PMLElementBase<2>::
compute_laplace_matrix_and_det(const DenseComplexMatrix& J, 
                               DenseComplexMatrix& laplace_matrix, 
                               std::complex<double>& detJ)
{
  const unsigned DIM = 2;
#ifdef PARANOID
  // check matrix dimensions are compatable and the DIM has been implemented
  if ( J.nrow() != DIM && J.ncol() != DIM )
  {
     throw OomphLibError("Input matrix must be 2x2", OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
  }
#endif

  detJ = J(0,0)*J(1,1) - J(1,0)*J(0,1);

  // resize and intialize result
  laplace_matrix.resize(DIM, DIM, 0.0);

  std::complex<double> invDetJ = 1.0/detJ;
  laplace_matrix(0,0) = (std::pow(J(0,1),2) + std::pow(J(1,1),2))*invDetJ; 
  laplace_matrix(1,1) = (std::pow(J(1,0),2) + std::pow(J(0,0),2))*invDetJ;
  laplace_matrix(0,1) = -(J(0,0)*J(0,1) + J(1,0)*J(1,1))*invDetJ;
  laplace_matrix(1,0) = laplace_matrix(0,1);
}

template<> 
void PMLElementBase<2>::
compute_laplace_matrix_and_det(const DiagonalComplexMatrix& J, 
                               DiagonalComplexMatrix& laplace_matrix, 
                               std::complex<double>& detJ)
{
  const unsigned DIM = 2;

  detJ = J(0,0) * J(1,1);

  // resize and intialize result
  laplace_matrix.resize(DIM, DIM, 0.0);

  laplace_matrix(0,0) = J(1,1) / J(0,0);
  laplace_matrix(1,1) = J(0,0) / J(1,1);
}

template<> 
void PMLElementBase<2>::
compute_laplace_matrix_and_det(const DiagonalComplexMatrix& J, 
                               DenseComplexMatrix& laplace_matrix, 
                               std::complex<double>& detJ)
{
  const unsigned DIM = 2;

  detJ = J(0,0) * J(1,1);

  // resize and intialize result
  laplace_matrix.resize(DIM, DIM, 0.0);

  laplace_matrix(0,0) = J(1,1) / J(0,0);
  laplace_matrix(1,1) = J(0,0) / J(1,1);
}

template<> 
void PMLElementBase<2>::
compute_jacobian_inverse_and_det(const DiagonalComplexMatrix& J, 
                               DiagonalComplexMatrix& J_inv, 
                               std::complex<double>& J_det)
{
  const unsigned DIM = 2;

  J_det = J(0,0) * J(1,1);

  // resize and intialize result
  J_inv.resize(DIM, DIM, 0.0);

  J_inv(0,0) = 1.0 / J(0,0);
  J_inv(1,1) = 1.0 / J(1,1);
}


template<> 
void PMLElementBase<3>::
compute_laplace_matrix_and_det(const DiagonalComplexMatrix& J, 
                               DiagonalComplexMatrix& laplace_matrix, 
                               std::complex<double>& detJ)
{
  const unsigned DIM = 3;

  detJ = J(0,0) * J(1,1) * J(2,2);

  // resize and intialize result
  laplace_matrix.resize(DIM, DIM, 0.0);

  laplace_matrix(0,0) = J(1,1) * J(2,2) / J(0,0);
  laplace_matrix(1,1) = J(0,0) * J(2,2) / J(1,1);
  laplace_matrix(2,2) = J(0,0) * J(1,1) / J(2,2);
}

template<> 
void PMLElementBase<3>::
compute_laplace_matrix_and_det(const DiagonalComplexMatrix& J, 
                               DenseComplexMatrix& laplace_matrix, 
                               std::complex<double>& detJ)
{
  const unsigned DIM = 3;

  detJ = J(0,0) * J(1,1) * J(2,2);

  // resize and intialize result
  laplace_matrix.resize(DIM, DIM, 0.0);

  laplace_matrix(0,0) = J(1,1) * J(2,2) / J(0,0);
  laplace_matrix(1,1) = J(0,0) * J(2,2) / J(1,1);
  laplace_matrix(2,2) = J(0,0) * J(1,1) / J(2,2);
}

template<> 
void PMLElementBase<3>::
compute_laplace_matrix_and_det(const DenseComplexMatrix& J, 
                               DenseComplexMatrix& laplace_matrix, 
                               std::complex<double>& detJ)
{
  const unsigned DIM = 3;

  detJ = J(0,0) * J(1,1) * J(2,2);

  // resize and intialize result
  laplace_matrix.resize(DIM, DIM, 0.0);

  throw OomphLibError(
    "Calculating Laplace matrix and det for 3x3 dense matrix not implemented",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
}

template<> 
void PMLElementBase<3>::
compute_jacobian_inverse_and_det(const DiagonalComplexMatrix& J, 
                               DiagonalComplexMatrix& J_inv, 
                               std::complex<double>& J_det)
{
  const unsigned DIM = 3;

  J_det = J(0,0) * J(1,1) * J(2,2);

  // resize and intialize result
  J_inv.resize(DIM, DIM, 0.0);

  J_inv(0,0) = 1.0 / J(0,0);
  J_inv(1,1) = 1.0 / J(1,1);
  J_inv(2,2) = 1.0 / J(2,2);
}

//=======================================================================
/// \short Compute n'J^{-1}, which when applied to the untransformed
/// grad operator gives us the transformed normal derivative
//=======================================================================
void TangentiallyDiscontinuousConformal2DPMLElement::compute_transformed_normal_derivative(
  const DenseComplexMatrix& J,
  const Vector<double>& n,
  Vector<std::complex<double> >& result)
{

  const unsigned DIM = 2;
#ifdef PARANOID
  unsigned long N = J.nrow();
  unsigned long n_col = J.ncol();
  unsigned long n_vec = n.size();
  // check matrix dimensions are compatable and the DIM has been implemented
  if ( N != n_col || N != n_vec || N > 2 )
  {
    throw OomphLibError(
      "Input vector must have size 2 and matrix must have size 2x2 ",
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // resize and intialize result
  result.resize(DIM, 0.0);

  const std::complex<double> detJ = J(0,0)*J(1,1) - J(1,0)*J(0,1);
  result[0] = (n[0]*J(1,1) - n[1]*J(0,1))/detJ;
  result[1] = (-n[0]*J(1,0) + n[1]*J(0,0))/detJ;

}

template <unsigned DIM>
void AxisAlignedPMLElement<DIM>::pml_transformation_jacobian(
  const unsigned& ipt,
  const Vector<double>& s,
  const Vector<double>& x,
  DiagonalComplexMatrix& pml_jacobian
)
{
  if (this->Pml_is_enabled)
  {  
    for (unsigned i=0; i<DIM; i++)
    {
      if(this->Pml_direction_active[i])
      {
        const double k = this->wavenumber();
        const double nu_i = nu(x,i);
        const double delta_i = delta(i);
        pml_jacobian(i,i) = Pml_mapping_pt->dtransformed_nu_dnu(nu_i, delta_i,
                                                                k);
      }
      else
      {
        pml_jacobian(i,i) = 1.0;
      }
    }
  }
  else
  {
    for (unsigned i=0; i<DIM; i++)
    {
      pml_jacobian(i,i) = 1.0;
    }
  }
}

template <unsigned DIM>
void AxisAlignedPMLElement<DIM>::pml_transformation_jacobian(
  const unsigned& ipt,
  const Vector<double>& s,
  const Vector<double>& x,
  DiagonalComplexMatrix& pml_jacobian,
  Vector<std::complex<double> >& transformed_x
)
{
  if (this->Pml_is_enabled)
  {  
    for (unsigned i=0; i<DIM; i++)
    {
      if(this->Pml_direction_active[i])
      {
        const double k = this->wavenumber();
        const double nu_i = nu(x,i);
        const double delta_i = delta(i);
        pml_jacobian(i,i) = Pml_mapping_pt->dtransformed_nu_dnu(nu_i, delta_i,
                                                                k);
        transformed_x[i] = Pml_inner_boundary[i]
                          + Pml_mapping_pt->transformed_nu(nu_i, delta_i, k);
      }
      else
      {
        pml_jacobian(i,i) = 1.0;
        transformed_x[i] = x[i];
      }
    }
  }
  else
  {
    for (unsigned i=0; i<DIM; i++)
    {
      pml_jacobian(i,i) = 1.0;
      transformed_x[i] = x[i];
    }
  }
}

void AnnularPMLElementBase::radial_to_cartesian_jacobian(
  const double& r,
  const double& theta,
  const std::complex<double>& rt,
  const std::complex<double>& drt_dr,
  DenseComplexMatrix& cartesian_jacobian) const
{
    DenseComplexMatrix radial_jacobian(2,2);
    radial_jacobian(0,0) = drt_dr;
    radial_jacobian(0,1) = 0.0;
    radial_jacobian(1,0) = 0.0;
    radial_jacobian(1,1) = 1.0;

    radial_to_cartesian_jacobian(r, theta, rt, radial_jacobian, cartesian_jacobian);
}

void AnnularPMLElementBase::radial_to_cartesian_jacobian(
  const double& r,
  const double& theta,
  const std::complex<double>& rt,
  const DenseComplexMatrix& radial_jacobian,
  DenseComplexMatrix& cartesian_jacobian) const
{
  
    // construct jacobian to move to transformed cartesian from transformed radial
    DenseComplexMatrix jacobian_xt_from_rt(2,2);
    jacobian_xt_from_rt(0,0) = cos(theta);
    jacobian_xt_from_rt(1,0) = sin(theta);
    jacobian_xt_from_rt(0,1) = -rt*sin(theta);
    jacobian_xt_from_rt(1,1) = rt*cos(theta);

    // construct jacobian to move to radial from cartesian
    DenseComplexMatrix jacobian_r_from_x(2,2);
    jacobian_r_from_x(0,0) = cos(theta);
    jacobian_r_from_x(0,1) = sin(theta);
    jacobian_r_from_x(1,0) = -sin(theta)/r;
    jacobian_r_from_x(1,1) = cos(theta)/r;

    // Multiply them all together
    // \f$(J_{\tilde x \tilde r}) (J_{\tilde r r}) (J_{\tilde r r})\f$

    // \f$(J_{\tilde x \tilde r}) (J_{\tilde r r})\f$
    DenseComplexMatrix temp(2,2);
    jacobian_xt_from_rt.multiply(radial_jacobian, temp);

    // Multiply previous by \f$(J_{\tilde r r})\f$, returning as mapping_jacobian
    temp.multiply(jacobian_r_from_x, cartesian_jacobian);
}

/// \short Implements an interface between the polar mapping and the Cartesian
/// output which the get residuals expects
void AnnularPMLElementBase::pml_transformation_jacobian(
  const unsigned& ipt,
  const Vector<double>& s,
  const Vector<double>& x,
  DenseComplexMatrix& cartesian_jacobian)
{
  // Get position from 
  const double theta = this->theta(s, x);
  const double r = this->radius(s, x);
  const double nu = this->nu(s, x);
  const double delta = this->delta();

  // Get the radial mapping jacobian and the transformed radial coordinate
  std::complex<double> dtr_dr;
  std::complex<double> tr;
  radial_transformation_jacobian(nu, delta, tr, dtr_dr);

  // Convert the radial Jacobian to a cartesian Jacobian
  radial_to_cartesian_jacobian(r, theta, tr, dtr_dr, cartesian_jacobian);
}

/// \short Implements an interface between the polar mapping and the Cartesian
/// output which the get residuals expects
void AnnularPMLElementBase::pml_transformation_jacobian(
  const unsigned& ipt,
  const Vector<double>& s,
  const Vector<double>& x,
  DenseComplexMatrix& cartesian_jacobian,
  Vector<std::complex<double> >& transformed_x)
{
  // Get position from 
  const double theta = this->theta(s, x);
  const double r = this->radius(s, x);
  const double nu = this->nu(s, x);
  const double delta = this->delta();

  // Get the radial mapping jacobian and the transformed radial coordinate
  std::complex<double> dtr_dr;
  std::complex<double> tr;
  radial_transformation_jacobian(nu, delta, tr, dtr_dr);

  // Convert the radial Jacobian to a cartesian Jacobian
  radial_to_cartesian_jacobian(r, theta, tr, dtr_dr, cartesian_jacobian);

  transformed_x[0] = tr * cos(theta);
  transformed_x[1] = tr * sin(theta);
}

void Conformal2DPMLElement::
get_pml_properties(const Vector<double>& s,
                   Vector<double>& x_inner,
                   Vector<double>& x_outer,
                   Vector<double>& p,
                   Vector<double>& dx_inner_dacross,
                   Vector<double>& dp_dacross,
                   double& delta)
{
  const unsigned N_NODE = this->nnode();
  const unsigned NU_S_INDEX = this->nu_s_index();
  const unsigned ACROSS_S_INDEX = this->across_s_index();
  const unsigned DIM = 2;

  // Vectors to hold local coordinates for the inner and outer parts of the PML
  Vector<double> s_inner(DIM);
  Vector<double> s_outer(DIM);

  // Get the values of nu normalised at the inner and outer boundaries of
  // this element from the element data
  const double nun_inner = Nu_at_s_min;
  const double nun_outer = Nu_at_s_max;

  // We use this to find the inner and outer position of the PML
  const double nun_det = nun_outer*(1.0-nun_inner) - nun_inner*(1.0-nun_outer);

  // Set LOCAL coordinate of point projected to the inner element boundary
  s_inner[ACROSS_S_INDEX] = s[ACROSS_S_INDEX];
  s_inner[NU_S_INDEX] = this->s_min();

  //Set up memory for the shape functions at the inside of the element
  Shape inner_psi(N_NODE);
  DShape inner_dpsi_ds(N_NODE, DIM);
  dshape_local(s_inner, inner_psi, inner_dpsi_ds);

  // Set LOCAL coordinate of point projected to the outer element boundary
  s_outer[ACROSS_S_INDEX] = s[ACROSS_S_INDEX];
  s_outer[NU_S_INDEX] = this->s_max();

  //Set up memory for the shape functions at the inside of the element
  Shape outer_psi(N_NODE);
  DShape outer_dpsi_ds(N_NODE, DIM);
  dshape_local(s_outer, outer_psi, outer_dpsi_ds);

  // Now project back and forward to the PML boundaries to get x_inner and
  // x_outer (we solve a 2x2 linear system twice) and them dacross
  // and delta which is thickness of PML at this point
  delta = 0.0;
  for(unsigned i=0; i<DIM; i++)
  {
    double x_el_inner_i = 0.0;
    double x_el_outer_i = 0.0;
    double dx_el_inner_dacross_i = 0.0;
    double dx_el_outer_dacross_i = 0.0;
    for(unsigned n=0; n<N_NODE; n++)
    {
      const double node_x = this->nodal_position(n,i);  
      x_el_inner_i += node_x*inner_psi(n);
      x_el_outer_i += node_x*outer_psi(n);
      dx_el_inner_dacross_i += node_x*inner_dpsi_ds(n, ACROSS_S_INDEX);
      dx_el_outer_dacross_i += node_x*outer_dpsi_ds(n, ACROSS_S_INDEX);
    }

    // Calculate the inner and outer coordinates of the PML
    x_inner[i] = (nun_outer*x_el_inner_i - nun_inner*x_el_outer_i)/nun_det;
    x_outer[i] = ( (1.0-nun_inner)*x_el_outer_i 
                  -(1.0-nun_outer)*x_el_inner_i )/nun_det;
    
    // Calculate unnormalised through-the-PML vector (we will normalise later)
    p[i] = x_outer[i] - x_inner[i];

    dx_inner_dacross[i] = ( nun_outer*dx_el_inner_dacross_i
                          - nun_inner*dx_el_outer_dacross_i)/nun_det;

    // Use dp_dacross to hold dx_outer_dacross, we will compute dp_dacross later
    dp_dacross[i] = ( (1.0-nun_inner)*dx_el_outer_dacross_i
                     -(1.0-nun_outer)*dx_el_inner_dacross_i )/nun_det;

    delta += p[i]*p[i];

    // faire: should we calculate p and dp_dacross in here??
    // Also, does delta change?
  }
  delta = sqrt(delta);

  double ddelta_dacross = 0.0;
  for(unsigned i=0; i<DIM; i++)
  {
    // Normalise the through-the-PML vector
    p[i] = p[i]/delta;
    // Recall that dp_dacross is currently holding dx_outer_dacross
    ddelta_dacross +=  p[i]*(dp_dacross[i] - dx_inner_dacross[i]);
  }

  // Finally calculate dp_dacross
  for(unsigned i=0; i<DIM; i++)
  {
    // Recall that dp_dacross is currently holding dx_outer_dacross
    dp_dacross[i] = (dp_dacross[i]-dx_inner_dacross[i]-ddelta_dacross*p[i])
                      /delta;
  }
}


void Conformal2DPMLElement::pml_transformation_jacobian(
  const unsigned& ipt,
  const Vector<double>& s,
  const Vector<double>& x,
  DenseComplexMatrix& jacobian
)
{
  const unsigned DIM = 2;

  double delta;
  Vector<double> x_inner(DIM);
  Vector<double> x_outer(DIM);
  Vector<double> p(DIM);
  Vector<double> dx_inner_dacross(DIM);
  Vector<double> dp_dacross(DIM);
  get_pml_properties(s, x_inner, x_outer, p, dx_inner_dacross, dp_dacross,
                     delta);
  
  const double NU = nu(s);
  
  const double k = wavenumber();
  std::complex<double> tnu = Pml_mapping_pt->transformed_nu(NU, delta, k);
  std::complex<double> dtnu_dnu = Pml_mapping_pt->dtransformed_nu_dnu(NU, delta,
                                                                      k);

  assemble_conformal_jacobian(p, dx_inner_dacross, dp_dacross, NU, tnu,
                              dtnu_dnu, jacobian);
}

void Conformal2DPMLElement::pml_transformation_jacobian(
  const unsigned& ipt,
  const Vector<double>& s,
  const Vector<double>& x,
  DenseComplexMatrix& jacobian,
  Vector<std::complex<double> >& transformed_x
)
{
  const unsigned DIM = 2;

  double delta;
  Vector<double> x_inner(DIM);
  Vector<double> x_outer(DIM);
  Vector<double> p(DIM);
  Vector<double> dx_inner_dacross(DIM);
  Vector<double> dp_dacross(DIM);
  get_pml_properties(s, x_inner, x_outer, p, dx_inner_dacross, dp_dacross,
                     delta);
  
  const double NU = nu(s);
  
  const double k = wavenumber();
  std::complex<double> tnu = Pml_mapping_pt->transformed_nu(NU, delta, k);
  std::complex<double> dtnu_dnu = Pml_mapping_pt->dtransformed_nu_dnu(NU, delta,
                                                                      k);

  assemble_conformal_jacobian(p, dx_inner_dacross, dp_dacross, NU, tnu,
                              dtnu_dnu, jacobian);
  
  transformed_x[0] = x_inner[0] + tnu*p[0]; 
  transformed_x[1] = x_inner[1] + tnu*p[1]; 
}


void Conformal2DPMLElement::assemble_conformal_jacobian(
  const Vector<double>& p,
  const Vector<double>& dx_inner_dacross,
  const Vector<double>& dp_dacross,
  const double& nu,
  const std::complex<double>& tnu,
  const std::complex<double>& dtnu_dnu,
  DenseComplexMatrix& mapping_jacobian)
{
  DenseComplexMatrix dtx_ds(2,2);
  dtx_ds(0,0) = dx_inner_dacross[0] + tnu*dp_dacross[0];
  dtx_ds(0,1) = dtnu_dnu*p[0];
  dtx_ds(1,0) = dx_inner_dacross[1] + tnu*dp_dacross[1];
  dtx_ds(1,1) = dtnu_dnu*p[1];

  DenseComplexMatrix ds_dx(2,2);
  double inv_det_dx_ds = 1.0/( (dx_inner_dacross[0] + dp_dacross[0]*nu)*p[1]
                              -(dx_inner_dacross[1] + dp_dacross[1]*nu)*p[0]);
  ds_dx(0,0) =  p[1]*inv_det_dx_ds;
  ds_dx(0,1) = -p[0]*inv_det_dx_ds;
  ds_dx(1,0) = -(dx_inner_dacross[1] + dp_dacross[1]*nu)*inv_det_dx_ds;
  ds_dx(1,1) =  (dx_inner_dacross[0] + dp_dacross[0]*nu)*inv_det_dx_ds;

  dtx_ds.multiply(ds_dx,mapping_jacobian);
}

void TangentiallyVaryingConformal2DPMLElement::assemble_conformal_jacobian(
  const Vector<double>& p,
  const Vector<double>& dx_inner_dacross,
  const Vector<double>& dp_dacross,
  const double& nu,
  const std::complex<double>& tnu,
  const std::complex<double>& dtnu_dnu,
  const std::complex<double>& dtnu_dacross,
  DenseComplexMatrix& mapping_jacobian)
{
  DenseComplexMatrix dtx_ds(2,2);
  dtx_ds(0,0) = dx_inner_dacross[0] + tnu*dp_dacross[0] + p[0]*dtnu_dacross;
  dtx_ds(0,1) = dtnu_dnu*p[0];
  dtx_ds(1,0) = dx_inner_dacross[1] + tnu*dp_dacross[1] + p[1]*dtnu_dacross;
  dtx_ds(1,1) = dtnu_dnu*p[1];

  DenseComplexMatrix ds_dx(2,2);
  double inv_det_dx_ds = 1.0/( (dx_inner_dacross[0] + dp_dacross[0]*nu)*p[1]
                              -(dx_inner_dacross[1] + dp_dacross[1]*nu)*p[0]);
  ds_dx(0,0) =  p[1]*inv_det_dx_ds;
  ds_dx(0,1) = -p[0]*inv_det_dx_ds;
  ds_dx(1,0) = -(dx_inner_dacross[1] + dp_dacross[1]*nu)*inv_det_dx_ds;
  ds_dx(1,1) =  (dx_inner_dacross[0] + dp_dacross[0]*nu)*inv_det_dx_ds;

  dtx_ds.multiply(ds_dx,mapping_jacobian);
}

template class AxisAlignedPMLElement<1>;
template class AxisAlignedPMLElement<2>;
template class AxisAlignedPMLElement<3>;

}
