//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1259 $
//LIC//
//LIC// $LastChangedDate: 2016-11-08 23:52:03 +0000 (Tue, 08 Nov 2016) $
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
// Driver for a specific 2D Helmholtz problem with
// perfectly matched layer treatment for the exterior boundaries

#include <fenv.h>

#include "math.h"
#include <complex>

#include <functional>

// Generic routines
#include "generic.h"

// The Helmholtz equations
#include "helmholtz.h"

// The pml Helmholtz equations
#include "pml_helmholtz.h"

// The meshes needed in the PML constructions
#include "meshes/triangle_mesh.h"
#include "meshes/rectangular_quadmesh.h"
#include "meshes/annular_mesh_with_pml.h"

// Get the Bessel functions
#include "oomph_crbond_bessel.h"

// Boundary ids for triangle mesh
#define INNER_SEMI_CIRCLE1 0
#define INNER_SEMI_CIRCLE2 1
#define MIDDLE_SEMI_CIRCLE1 2
#define MIDDLE_SEMI_CIRCLE2 3
#define OUTER_SEMI_CIRCLE1 4
#define OUTER_SEMI_CIRCLE2 5

// Boundary ids for structured mesh
#define INNER_CIRCLE 0
#define OUTER_CIRCLE 2
#define MIDDLE_CIRCLE 4

using namespace oomph;
using namespace std;

using namespace MathematicalConstants;
using namespace PowerSeriesPMLMappingFillIns;


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//===== start_of_namespace=============================================
/// Namespace for the Helmholtz problem parameters
//=====================================================================
namespace GlobalParameters
{

/// Wavenumber (also known as k), k=omega/c
double Wavenumber = 0.1;

/// Square of the wavenumber (also known as k^2)
double K_squared = Wavenumber * Wavenumber;

/// Inner radius of domain
double Inner_r = 1.0;

/// Outer radius of bulk, inner radius of PML
double Middle_r = 2.0;

/// Outer radius of PML
double Outer_r = 3.0;

/// Hankel function/Fourier mode to use in single hankel exact solution
int Single_hankel_order=1;

/// The number of segments on the outer boundary
unsigned N_node = 2;

/// The number of segments on the outer boundary
unsigned N_angular_elements = 16;

/// Target element area to feed to triangle
double Triangle_element_area = 0.01;

// Number of PML elements (for structured mesh)
unsigned N_pml = 10;

// Number of elements in the radial direction (for structured mesh)
unsigned N_radial_element = (unsigned) ceil(N_angular_elements* (Middle_r - Inner_r)/((Middle_r + Inner_r)*Pi));

/// The number of elements from command line, if it is zero then set auto
unsigned Set_num_n_radial_element = 0;

/// Perturbation to radial uniformness of mesh
double Mesh_perturbation = 0.0;

/// Set th exact solution on the outer boundary, useful if you want to test without PML
bool Set_exact_on_outer = false;

/// Should we do a test where we turn off the mapping in the PML?
bool Null_mapping_test = false;

/// Use a completely made up but simple mapping in the PML?
bool Dummy_mapping_test = false;

/// The mode we will test with
int Single_fourier_mode = 0;

/// Phase shift of the exact solution
double Phase_shift = 0.0;

/// Should suppress output of soln, exact solution and coarse solution?
bool Suppress_doc_solution = false;

// Decide which exact solution to use
bool Use_single_mode_solution = false;
bool Use_single_hankel_solution = false;
bool Use_exact_soln_from_file = false;

// Reads in coefficients from a file in the same directory as the driver code,
// where the real and imaginary parts of each coefficient are separated by a
// start_of_namespace
void read_coefs_from_file(const char filename[],
  Vector<std::complex<double> >& coef_vec)
{
  coef_vec.resize(0);
  ifstream file;
  file.open(filename);
  // Keep reading real and imaginary parts until end of file
  double coef_real, coef_imag;
  while(file >> coef_real && file >> coef_imag){
    coef_vec.push_back(std::complex<double>(coef_real,coef_imag));
  }
  file.close();
}

// -------------------------------------------------------------------------- //
//                           Solution from file
// -------------------------------------------------------------------------- //

Vector<std::complex<double> > Exact_soln_coefs_from_file(0);

// Reads in coefficients from a file called "exact_soln_coefficients.dat" in the
// same directory as the driver code, where the real and imaginary parts of each
// coefficient are separated by a start_of_namespace
void read_exact_soln_coefs_from_file(){

  read_coefs_from_file("exact_soln_coefficients.dat",Exact_soln_coefs_from_file);

  oomph_info << "Read in " << Exact_soln_coefs_from_file.size()
             << " exact solution coefficients from file" << std::endl;

}


/// \short Exact solution using coefficients form file using r and theta
std::complex<double> exact_soln_from_file_r_theta(const unsigned& theta_deriv,
  const double& r, const double& theta, bool normal_deriv=false)
{

  // If we are in the PML, and we aren't turning off the PML (null mapping test)
  // then return the "ideal" profile
  double pml_multiplier = 1.0;
  double r_local = r;
  if(r > Middle_r && !Null_mapping_test)
  {
    // Evaluate at inner PML boundary
    r_local = Middle_r;
    //  1 - \bar \nu
    pml_multiplier =  1.0 - (r - Middle_r)/(Outer_r - Middle_r);
  }

  // Argument for Bessel/Hankel functions
  double k = Wavenumber;
  double kr=k*r_local;

  if(r > Middle_r && Dummy_mapping_test)
  {
    // DenseComplexMatrix mapping_jacobian_radial(2,2);
    // std::complex<double> rt; // r transformed
    // Multiaxial_pml_mapping_pt->
    //   get_mapping_jacobian_radial(r,theta,K_squared,mapping_jacobian_radial,rt);
    double A = Dummy_nonaxisym_amplitude;
    double B = Dummy_mapping_complex_amplitude;
    double nun = (r - Middle_r)/(Outer_r - Middle_r);
    std::complex<double> rt = Middle_r + (1.0+A*sin(theta))*(nun+I*B*nun*nun);
    kr = k * real(rt);
    pml_multiplier = 1.0;
  }

  unsigned N_hankel_max = (Exact_soln_coefs_from_file.size()-1)/2;

  Vector<complex<double> > h(N_hankel_max+1), hp(N_hankel_max+1);

  // Evaluate Hankel at actual radius
  using namespace Hankel_functions_for_helmholtz_problem;
  Hankel_first(N_hankel_max,kr,h,hp);

  complex<double> u(0.0,0.0);

  // Sum over all coefficients and related Hankel functions
  // with a different Hankel if we want the normal derivative
  for (int n=-N_hankel_max; n<=int(N_hankel_max); n++)
  {
    double hankel_sign = 1.0;
    if ((n < 0) && (n % 2 != 0)) hankel_sign = -1.0;
    unsigned n_abs = unsigned(std::abs(n));
    u += Exact_soln_coefs_from_file[unsigned(N_hankel_max+n)]
         * pow(I*double(n),theta_deriv) * exp(I*theta*double(n))
         * hankel_sign * (normal_deriv ? k*hp[n_abs] : h[n_abs]);
  }

  return u*pml_multiplier;

} // end of exact_soln_from_file_r_theta

/// \short Exact solution using coefficients form file (vector returns real and
// imaginary parts).
void get_exact_soln_from_file(const Vector<double>& x, Vector<double>& u)
{
  // Switch to polar coordinates
  double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  double theta=atan2(x[1],x[0]);

  // Get solution using version with r and theta
  unsigned theta_deriv = 0;
  complex<double> u_ex = exact_soln_from_file_r_theta(theta_deriv,r,theta);

  // Get the real & imaginary part of the result
  u[0]=real(u_ex);
  u[1]=imag(u_ex);

} // end of get_exact_u_from_file


// Value at inner pml boundary
std::complex<double> exact_U_from_file(const int& theta_deriv, const double& theta){
  bool normal_deriv = false;
  return exact_soln_from_file_r_theta(theta_deriv, Middle_r, theta, normal_deriv);
}
// Normal derivative at inner pml boundary
std::complex<double> exact_N_from_file(const int& theta_deriv, const double& theta){
  bool normal_deriv = true;
  return exact_soln_from_file_r_theta(theta_deriv, Middle_r, theta, normal_deriv);
}


// -------------------------------------------------------------------------- //
//                           Zero solution
// -------------------------------------------------------------------------- //

/// \short Exact solution which is 0
void get_zero_solution(const Vector<double>& x, Vector<double>& u)
{
  u[0]=0.0;
  u[1]=0.0;
}

// -------------------------------------------------------------------------- //
//                           Single (Fourier) mode
// -------------------------------------------------------------------------- //

/// \short exact solution for a single fourier mode at a certain theta
complex<double> single_mode_theta(const int& theta_deriv, const double& theta)
{
  // Holder for solution
  complex<double> u;
  // Modification from shifting phase
  complex<double> u_shift = exp(I*Phase_shift);

  if(Single_fourier_mode == 0 && theta_deriv == 0)
  {
    u = u_shift * complex<double>(1.0,0.0);
  }
  else
  {
    u = u_shift * pow(double(Single_fourier_mode)*I,theta_deriv)*exp(Single_fourier_mode*theta*I);
  }

  return u;
}

/// \short exact solution for a single fourier mode at a vector x
complex<double> single_mode(const int& theta_deriv, const Vector<double>& x)
{
  double theta = atan2(x[1],x[0]);
  double r = sqrt(x[0]*x[0] + x[1]*x[1]);
  // r dependence give ssolution unit gradient at inner pml (Middle_r)
  return single_mode_theta(theta_deriv, theta)*(r - Middle_r + 1);
}

/// \short exact solution for a single fourier mode (0th derivative)
/// (vector returns real and imaginary parts).
void get_single_mode(const Vector<double>& x, Vector<double>& u)
{
  complex<double> u_ex = single_mode(0,x);
  u[0]=real(u_ex);
  u[1]=imag(u_ex);
}

// -------------------------------------------------------------------------- //
//                           Single Hankel mode
// -------------------------------------------------------------------------- //

std::complex<double> exact_single_hankel_r_theta(
  const int& theta_deriv, const double& r, const double& theta,
  bool normal_deriv=false)
{

  // If we are in the PML, and we aren't turning off the PML (null mapping test)
  // then return the "ideal" profile
  double pml_multiplier = 1.0;
  double r_local = r;
  if(r > Middle_r && !Null_mapping_test)
  {
    // Evaluate at inner PML boundary
    r_local = Middle_r;
    //  1 - \bar \nu
    pml_multiplier =  1.0 - (r - Middle_r)*(Outer_r - Middle_r);
  }

  // Make maths a bit nicer by caching wavenumber
  double k = Wavenumber;

  // Argument for Bessel/Hankel functions
  double kr=k*r_local;

  // Evaluate Bessel/Hankel functions
  complex <double > u_ex(0.0,0.0);

  // Number of Hankel functions to fetch
  unsigned N_hankel = (unsigned) abs(Single_hankel_order);

  Vector< std::complex<double> > h(N_hankel+1);
  Vector< std::complex<double> > hp(N_hankel+1);
  Vector< std::complex<double> > h_a(N_hankel+1);
  Vector< std::complex<double> > hp_a(N_hankel+1);

  // Evaluate Hankel at actual radius
  Hankel_functions_for_helmholtz_problem::Hankel_first(N_hankel,kr,h,hp);

  // Evaluate Hankel at inner (unit) radius for normalisation
  Hankel_functions_for_helmholtz_problem::Hankel_first(N_hankel,k*Middle_r,h_a,hp_a);


  // Solution without the Hankel bit that depends on r
  std::complex<double> u = exp(I*theta*double(Single_hankel_order))/h_a[N_hankel];
  // std::complex<double> u = exp(I*theta*double(Single_hankel_order));

  if(N_hankel != 0) u *= pow(I*double(N_hankel),double(theta_deriv));
  if(N_hankel == 0 && theta_deriv >= 1) u = 0;

  // For the normal derivative, we use hp otherwise just use regular hankel
  if(normal_deriv) u *= k*hp[N_hankel];
  else u *= h[N_hankel];

  return u*pml_multiplier;
} // end of get_exact_single_hankel_complex_r_theta

void get_exact_single_hankel(const Vector<double>& x, Vector<double>& u_vec)
{
  // Get the function itself (no deriv, 0) because this version is used for setting
  // the exact solution at nodal points
  double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  double theta=atan2(x[1],x[0]);
  std::complex<double> u = exact_single_hankel_r_theta(0,r,theta);

  // Assign real and imaginary parts to a vector as oomph-lib expects
  u_vec[0] = u.real();
  u_vec[1] = u.imag();
}

/// If no r is given, assume we are at the inner edge of PML
std::complex<double> exact_single_hankel_theta(const int& theta_deriv,const double& theta)
{
  double assumed_r = Middle_r;
  return exact_single_hankel_r_theta(theta_deriv,assumed_r,theta);
}

/// If no r is given, assume we are at the inner edge of PML
std::complex<double> exact_normal_single_hankel_theta(const int& theta_deriv,const double& theta)
{
  double assumed_r = Middle_r;
  bool normal_deriv = true;
  return exact_single_hankel_r_theta(theta_deriv,assumed_r,theta,normal_deriv);
}

void update_parameters(){
  oomph_info << "Updating parameters dependent on command line args" << std::endl;
  // If it's a valid ratio, use it
  if(Fourier_diff_modes_n_angular_elements_ratio > 0.0){
    Fourier_diff_modes = (unsigned) min(
      10.0,
      ceil(N_angular_elements/Fourier_diff_modes_n_angular_elements_ratio - 1.01)
    );
  }
  oomph_info << "Using "<< Fourier_diff_modes << " modes around boundary to differentiate the solution." << std::endl;

  if(Set_num_n_radial_element == 0)
  {
    N_radial_element = (unsigned) ceil(N_angular_elements*(Middle_r - Inner_r)/((Middle_r + Inner_r)*Pi));
    oomph_info << "Automatically determined "<< N_radial_element << " elements in the radial direction." << std::endl;
  }
  else
  {
    N_radial_element = Set_num_n_radial_element;
    oomph_info << "Manually set "<< N_radial_element << " elements in the radial direction." << std::endl;
  }

  K_squared = Wavenumber*Wavenumber;
  oomph_info << "Using k^2 = "<< K_squared << ". " << std::endl;

  if (Pml_full_integration)
  {
    Pml_integral_pt = new AnySizeAnisotropic2DGaussLegendre(N_node,Pml_integration_points);
  }
}



} // end of namespace

//========= start_of_generate_circle===================================
/// Helper function which returns a mesh containing all the elements in
/// Whole_mesh_pt which fall within a certain annular region. Note, this
/// does not construct the elements and nodes, just hold spointers to them
//=====================================================================
Mesh* create_annular_error_mesh(const Mesh* Whole_mesh_pt, const double& r_inner,const double& r_outer)
{

  // Create our mesh which will hold the error_mesh
  Mesh* Error_mesh_pt = new Mesh();

  // Holders for determining centre of element
  Vector<double> s(2,0.0);
  Vector<double> x(2,0.0);

  // Find all the elements in the annulus and add to new mesh
  unsigned n_element = Whole_mesh_pt->nelement();
  for(unsigned i=0;i<n_element;i++){
    // Get ith element of the Whole mesh
    GeneralisedElement* el_pt = Whole_mesh_pt->element_pt(i);
    // Get centre of element
    dynamic_cast<FiniteElement*>(el_pt)->interpolated_x(s,x);

    if(x[0]*x[0]+x[1]*x[1] > r_inner*r_inner && x[0]*x[0]+x[1]*x[1] < r_outer*r_outer){
      Error_mesh_pt->add_element_pt(el_pt);
    }
  }

  // Loop over all the nodes in the whole mesh
  unsigned n_node = Whole_mesh_pt->nnode();
  for (unsigned n=0; n<n_node; n++)
  {
    // Get pointer to node
    Node* nod_pt=Whole_mesh_pt->node_pt(n);

    // Extract nodal coordinates from node:
    nod_pt->position(x);

    if(x[0]*x[0]+x[1]*x[1] >= r_inner*r_inner && x[0]*x[0]+x[1]*x[1] <= r_outer*r_outer){
      Error_mesh_pt->add_node_pt(nod_pt);
    }
  }

  // Return our new mesh
  return Error_mesh_pt;

}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//========= start_of_problem_class=====================================
/// Problem class to demonstrate use of perfectly matched layers
/// for Helmholtz problems.
//=====================================================================
template<class ELEMENT>
class PMLProblem : public Problem
{

public:

/// Constructor
PMLProblem();

/// Destructor (empty)
~PMLProblem(){
};

/// Update the problem specs before solve (empty)
void actions_before_newton_solve(){
};

/// Update the problem specs after solve (empty)
void actions_after_newton_solve(){
};

/// \short Doc the solution. DocInfo object stores flags/labels for where the
/// output gets written to
void doc_solution(DocInfo& doc_info);


/// \short Create Helmholtz flux elements on boundary b of the Mesh pointed
/// to by mesh_pt and add them to the specified surface Mesh
void create_flux_elements(const unsigned &b, Mesh* const &mesh_pt,
                          Mesh* const & helmholtz_inner_boundary_mesh_pt);

void apply_boundary_conditions(const bool set_exact_on_outer = false);

void set_solution_everywhere(FiniteElement::SteadyExactSolutionFctPt fct = 0);

void enable_pmls();

// Load coefficients, do a newton solve, doc solution thenn dump coefficients
void solve_and_doc(DocInfo& doc_info);

// Holder for  solution (if it doesn't exist, set to 0)
FiniteElement::SteadyExactSolutionFctPt exact_fct_pt;

private:

/// Pointer to the "bulk" mesh
#ifdef TRIANGLE
TriangleMesh<ELEMENT>* Mesh_pt;
#else
PMLTwoDAnnularMesh<ELEMENT>* Mesh_pt;
#endif

/// \short Pointer to the mesh containing
/// the Helmholtz inner boundary condition elements
Mesh* Helmholtz_inner_boundary_mesh_pt;

/// \short Pointer to the mesh containing elements not in the PML (physical region)
Mesh* Bulk_mesh_pt;

/// \short Pointer to the mesh containing element in PML
Mesh* Pml_mesh_pt;

/// Trace file
ofstream Trace_file;

}; // end of problem class


//=======start_of_constructor=============================================
/// Constructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLProblem<ELEMENT>::PMLProblem()
{

  // Open trace file
  Trace_file.open("RESLT/trace.dat");

#ifdef TRIANGLE
  // Define a target segment size to be used when constructing inner and outer
  // boundaries
  double target_segment_size = 2.0 * Pi
                               / GlobalParameters::N_angular_elements;

  // Define the centre of all the circles (the origin)
  Vector<double> centre(2,0.0);


  // Outer boundary
  //---------------
  TriangleMeshClosedCurve* outer_boundary_pt=generate_circle(
    GlobalParameters::Outer_r, centre,
    OUTER_SEMI_CIRCLE1,
    OUTER_SEMI_CIRCLE2, target_segment_size);

  // Inner boundary
  //---------------
  TriangleMeshClosedCurve* inner_boundary_pt=generate_circle(
    GlobalParameters::Inner_r, centre,
    INNER_SEMI_CIRCLE1,
    INNER_SEMI_CIRCLE2, target_segment_size);

  // Middle boundary (inner PML boundary)
  //---------------
  TriangleMeshClosedCurve* middle_boundary_pt=generate_circle(
    GlobalParameters::Middle_r, centre,
    MIDDLE_SEMI_CIRCLE1,
    MIDDLE_SEMI_CIRCLE2, target_segment_size);

  // Create mesh
  // -----------------------
  // Use the TriangleMeshParameters object for helping on the manage
  // of the TriangleMesh parameters. The only parameter that needs to take
  // is the outer boundary.
  TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

  // Specify the internal closed boundaries
  Vector<TriangleMeshClosedCurve*> inner_boundary_pt_vec(2);
  inner_boundary_pt_vec[0] = inner_boundary_pt;
  inner_boundary_pt_vec[1] = middle_boundary_pt;
  triangle_mesh_parameters.internal_closed_curve_pt() = inner_boundary_pt_vec;

  // Target element size in bulk mesh (make it massive so that it just uses
  // segments)
  triangle_mesh_parameters.element_area() = GlobalParameters::Triangle_element_area;

  // Build "bulk" mesh
  Mesh_pt=new TriangleMesh<ELEMENT>(triangle_mesh_parameters);
#else

  Mesh_pt = new PMLTwoDAnnularMesh<ELEMENT>(
    GlobalParameters::N_angular_elements,
    GlobalParameters::N_radial_element, GlobalParameters::N_pml,
    GlobalParameters::Inner_r, GlobalParameters::Middle_r,
    GlobalParameters::Outer_r, GlobalParameters::Mesh_perturbation);

#endif
  // Create the main mesh
  add_sub_mesh(Mesh_pt);

  // Build the entire mesh from its submeshes
  build_global_mesh();

  /// Create seperate meshes to evaulate the error
  Bulk_mesh_pt = create_annular_error_mesh(Mesh_pt, GlobalParameters::Inner_r,
    GlobalParameters::Middle_r);
  Pml_mesh_pt = create_annular_error_mesh(Mesh_pt, GlobalParameters::Middle_r,
    GlobalParameters::Outer_r);

  //                           Exact solutions
  // --------------------------------------------------------------------------
  if(GlobalParameters::Use_single_hankel_solution)
  {
    oomph_info << "Using a " <<  GlobalParameters::Single_hankel_order
               << " order Hankel solution" << std::endl;
    // Set problem exact function pointers to the single hankel ones
    exact_fct_pt = GlobalParameters::get_exact_single_hankel;
  }
  else if(GlobalParameters::Use_single_mode_solution)
  {
    oomph_info << "Using a " <<  GlobalParameters::Single_fourier_mode
               << " Fourier mode solution (cone)" << std::endl;
    // Set problem exact function pointers to the single mode ones. In this case
    // U and N are intentionally the same
    exact_fct_pt = GlobalParameters::get_single_mode;
  }
  else if(GlobalParameters::Use_exact_soln_from_file)
  {
    oomph_info << "Using exact solution using Hankel coefficients from file" << std:: endl;
    // Set problem exact function pointer to the one which uses coefs from file
    GlobalParameters::read_exact_soln_coefs_from_file();
    exact_fct_pt = GlobalParameters::get_exact_soln_from_file;
  }
  else
  {
    throw OomphLibError("No exact solution set.",OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }


  //                         Boundary conditions
  // ------------------------------------------------------------------------


  if(GlobalParameters::Set_exact_everywhere) set_solution_everywhere(exact_fct_pt);

  apply_boundary_conditions(GlobalParameters::Set_exact_on_outer);

  }

  // Now that we have set up the transformation, enable the elements
  enable_pmls();

  // Setup equation numbering scheme
  cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor

bool on_inner_boundary(const double& r)
{
  return std::abs(r - GlobalParameters::Middle_r) < 1.0e-12;
}

bool on_inner_boundary(const Vector<double>& x)
{
  return std::abs(x[0]*x[0] + x[1]*x[1] - GlobalParameters::Middle_r*GlobalParameters::Middle_r) < 1.0e-12;
}

bool on_outer_boundary(const double& r)
{
  return std::abs(r - GlobalParameters::Outer_r) < 1.0e-12;
}

bool on_outer_boundary(const Vector<double>& x)
{
  return std::abs(x[0]*x[0] + x[1]*x[1] - GlobalParameters::Outer_r*GlobalParameters::Outer_r) < 1.0e-12;
}

double theta(const Vector<double>& x)
{
  return atan2(x[1],x[0]);
}

double nun(const double& r)
{
  return (r - GlobalParameters::Middle_r)/(GlobalParameters::Outer_r - GlobalParameters::Middle_r);
}

double nun(const Vector<double>& x)
{
  return nun(sqrt(x[0]*x[0]+x[1]*x[1]));
}


//=========================start_of_enable_pmls===========================
/// Locate PML and enable correct PML elements, pass in correct pointer
//========================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::enable_pmls()
{
  
  // Complete the build of all elements so they are fully functional
  unsigned n_element = this->mesh_pt()->nelement();
  for(unsigned e=0; e<n_element; e++)
  {
    // Upcast from GeneralisedElement to Pml Helmholtz bulk element
    PMLHelmholtzEquations<2> *el_pt =
      dynamic_cast<PMLHelmholtzEquations<2>*>(mesh_pt()->element_pt(e));

    if (el_pt!=0)
    {
      //Set the k_squared function pointer
      el_pt->k_squared_pt() = &GlobalParameters::K_squared;

      // Centre of element in local coordinates
#ifdef TRIANGLE
      Vector<double> s(2,0.3);
#else
      Vector<double> s(2,0.0);
#endif

      // Get the global coordinates of centre of element
      Vector<double> x(2,0.0);
      el_pt->interpolated_x(s,x);

      // If in PML region
      if (x[0]*x[0] + x[1]*x[1] > GlobalParameters::Middle_r*GlobalParameters::Middle_r)
      {
        // Enable PML with dummy boundary values (we do not use the inbuilt
        // mechanism for calculating distance and direction in PML)
        // el_pt->enable_pml(1,0.0,1.0);

        // el_pt->pml_mapping_pt() = 0;
      }
    }
  }
}

//==================start_of_on_boundary functions========================
/// Little helper functions which tell us if we are on a certain boundary
//========================================================================
bool on_inner_boundary(const unsigned& b)
{
#ifdef TRIANGLE
  if(b == INNER_SEMI_CIRCLE1 || b == INNER_SEMI_CIRCLE2) return true;
#else
  if(b == INNER_CIRCLE) return true;
#endif
  else return false;
}

bool on_middle_boundary(const unsigned& b)
{
#ifdef TRIANGLE
  if(b == MIDDLE_SEMI_CIRCLE1 || b == MIDDLE_SEMI_CIRCLE2) return true;
#else
  if(b == MIDDLE_CIRCLE) return true;
#endif
  else return false;
}

bool on_outer_boundary(const unsigned& b)
{
#ifdef TRIANGLE
  if(b == OUTER_SEMI_CIRCLE1 || b == OUTER_SEMI_CIRCLE2) return true;
#else
  if(b == OUTER_CIRCLE) return true;
#endif
  else return false;
}

//==================start_of_apply_boundary_conditions====================
/// Apply boundary conditions
//========================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::apply_boundary_conditions(const bool set_exact_on_outer)
{

  // Boundary conditions are set on the surface of the circle
  // as a constant nonzero Dirichlet boundary condition
  unsigned n_bound = Mesh_pt->nboundary();

  for(unsigned b=0; b<n_bound; b++)
  {
    // Set dirichlet 0 BC on outer boundary if we aren't setting exact there
    if (on_outer_boundary(b) && !set_exact_on_outer)
    {
      unsigned n_node = Mesh_pt->nboundary_node(b);
      for (unsigned n=0; n<n_node; n++)
      {

        Node* nod_pt=Mesh_pt->boundary_node_pt(b,n);
        nod_pt->pin(0);
        nod_pt->pin(1);

        nod_pt->set_value(0,0.0);
        nod_pt->set_value(1,0.0);
      }
    }

    // Set dirichlet exact boundary conditions on inner boundary (and outer if asked to)
    if (on_inner_boundary(b) || (set_exact_on_outer && on_outer_boundary(b)))
    {
      unsigned n_node = Mesh_pt->nboundary_node(b);
      for (unsigned n=0; n<n_node; n++)
      {
        // Get pointer to node
        Node* nod_pt=Mesh_pt->boundary_node_pt(b,n);

        // Extract nodal coordinates from node:
        Vector<double> pos(2);
        nod_pt->position(pos);

        // Compute the value of the exact solution at the nodal point
        Vector<double> u(2);
        (*exact_fct_pt)(pos,u);

        nod_pt->pin(0);
        nod_pt->pin(1);

        // Assign the value
        nod_pt->set_value(0,u[0]);
        nod_pt->set_value(1,u[1]);

      }
    }
  }

} // end of apply_boundary_conditions

//==================start_of_apply_boundary_conditions====================
/// Set a solution everywhere without pinning the nodes, defaults to the
/// exact solution assigned to the problem if none is provided
//========================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::set_solution_everywhere(FiniteElement::SteadyExactSolutionFctPt fct_pt)
{
  // If no function provided, use the one assigned to the problem
  if(fct_pt == 0) fct_pt = this->exact_fct_pt;

  // Loop over all the nodes in the whole mesh
  unsigned n_node = Mesh_pt->nnode();
  for (unsigned n=0; n<n_node; n++)
  {
    // Get pointer to node
    Node* nod_pt=Mesh_pt->node_pt(n);

    // Extract nodal coordinates from node:
    Vector<double> pos(2);
    nod_pt->position(pos);

    // Compute the value of the exact solution at the nodal point
    Vector<double> u(2);
    (*fct_pt)(pos,u);

    // Assign the value
    nod_pt->set_value(0,u[0]);
    nod_pt->set_value(1,u[1]);

  }
} // end of apply_boundary_conditions

//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{

  ofstream some_file,some_file2, csv_file;
  char filename[100];

  // Number of plot points for fine and coarse solns
  unsigned npts=5;
  unsigned coarse_npts=2;

  if(!GlobalParameters::Suppress_doc_solution)
  {
    // Output solution
    //-----------------
    sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    Mesh_pt->output(some_file,npts);
    some_file.close();

    // Output coarse solution
    //-----------------
    sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    Mesh_pt->output(some_file,coarse_npts);
    some_file.close();

    // Output exact solution
    //----------------------
    sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    Mesh_pt->output_fct(some_file,npts,exact_fct_pt);
    some_file.close();

    // Output exact solution in bulk
    //----------------------
    sprintf(filename,"%s/exact_bulk_soln%i.dat",doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    Bulk_mesh_pt->output_fct(some_file,npts,exact_fct_pt);
    some_file.close();
  }

  double error,norm;
  sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->compute_error(some_file,exact_fct_pt,error,norm);
  some_file.close();

  // Doc L2 error and norm of solution
  oomph_info << std::endl;
  oomph_info << "Bulk error: " << sqrt(error) << std::endl;
  oomph_info << "Bulk norm: " << sqrt(norm) << std::endl;
  oomph_info << "Normalised bulk error: " << sqrt(error/norm) << std::endl;

  // Write norm of exact bulk, normalised error in bulk and norm of pml solution to trace file
  // Note that pml_error is actually the norm of the pml because we use the zero exact solution
  Trace_file << sqrt(norm) << ", "<< sqrt(error/norm) << std::endl;

} // end of doc

template<class ELEMENT> void PMLProblem<ELEMENT>::solve_and_doc(DocInfo& doc_info)
{
  // Problem is actually linear, no need to recompute residuals before and after
  // problem level Newton solve
  this->problem_is_nonlinear(false);

  // Solve the problem
  if(!GlobalParameters::Disable_solve) newton_solve();

  // Save results to file
  doc_solution(doc_info);
}

//==========start_of_main=================================================
/// Solve 2D Helmholtz problem
//========================================================================
int main(int argc, char **argv)
{

  // For 128 elements around boundary we get error, increase tolerance to ignore
  ToleranceForVertexMismatchInPolygons::Tolerable_error = 10e-12;

  //Read the parameters for the problem
  //-----------------------------------
  // Set up all the command line arguments linking them to global parameters, to
  // see what each does, see the definitions in GlobalParameters
  {
    CommandLineArgs::setup(argc,argv);

    using namespace CommandLineArgs;
    using namespace GlobalParameters;

    // Mesh options
    specify_command_line_flag("--inner_r",&Inner_r);
    specify_command_line_flag("--middle_r",&Middle_r);
    specify_command_line_flag("--outer_r",&Outer_r);
    specify_command_line_flag("--n_angular_elements",&N_angular_elements);
    specify_command_line_flag("--n_radial_elements",&Set_num_n_radial_element);
    specify_command_line_flag("--mesh_perturbation",&Mesh_perturbation);
    specify_command_line_flag("--triangle_element_area",&Triangle_element_area);
    specify_command_line_flag("--n_pml",&N_pml);
    specify_command_line_flag("--n_node",&N_node);

    // Exact solution options
    specify_command_line_flag("--set_exact_everywhere");
    specify_command_line_flag("--pin_bulk");
    specify_command_line_flag("--exact_soln_from_file");
    specify_command_line_flag("--single_hankel",&Single_hankel_order);

    // Diff mesh options
    specify_command_line_flag("--use_exact_in_mapping");
    specify_command_line_flag("--calculate_fourier_coefs_from_exact");
    specify_command_line_flag("--fourier_mode_ratio",&Fourier_diff_modes_n_angular_elements_ratio);
    specify_command_line_flag("--fourier_diff_modes",&Fourier_diff_modes);
    specify_command_line_flag("--dump_fourier_coefs");
    specify_command_line_flag("--read_fourier_coefs");

    specify_command_line_flag("--calculate_hankel_series_from_exact");

    // Mapping options
    specify_command_line_flag("--newton_mapping");
    specify_command_line_flag("--conformal_newton_mapping");
    specify_command_line_flag("--n_trunc",&N_trunc);
    specify_command_line_flag("--fill_in",&Fill_in);
    specify_command_line_flag("--debug_mapping");
    specify_command_line_flag("--dummy_nonaxisym_amplitude",&Dummy_nonaxisym_amplitude);
    specify_command_line_flag("--impose_straight_line_in_pml");

    // Other options
    specify_command_line_flag("--fourier_full_integration");
    specify_command_line_flag("--pml_integration_points",&Pml_integration_points);
    specify_command_line_flag("--wavenumber",&Wavenumber);
    specify_command_line_flag("--disable_solve");
    specify_command_line_flag("--suppress_doc_solution");

    // Run specific tests (these will overide other settings if required)
    specify_command_line_flag("--test_single_fourier_mode",&Single_fourier_mode);
    specify_command_line_flag("--test_diff");
    specify_command_line_flag("--null_mapping_test");
    specify_command_line_flag("--dummy_mapping_test");
    specify_command_line_flag("--objective_test",&Output_objective_theta_samples);
    specify_command_line_flag("--objective_test_theta_min",&Output_objective_theta_min);
    specify_command_line_flag("--objective_test_theta_max",&Output_objective_theta_max);

    // Parse command line
    parse_and_assign();

    // Read if any of the boolean flags have been specified (shouldn't bools just be like the others?)
    // Exact solution options
    Set_exact_everywhere=command_line_flag_has_been_set("--set_exact_everywhere");
    Pin_bulk=command_line_flag_has_been_set("--pin_bulk");
    Use_exact_soln_from_file=command_line_flag_has_been_set("--exact_soln_from_file");

    // Diff mesh options
    Use_exact_in_mapping=command_line_flag_has_been_set("--use_exact_in_mapping");
    Calculate_fourier_coefs_from_exact=command_line_flag_has_been_set("--calculate_fourier_coefs_from_exact");
    Dump_coefs_flag=command_line_flag_has_been_set("--dump_fourier_coefs");
    Read_coefs_flag=command_line_flag_has_been_set("--read_fourier_coefs");
    if(command_line_flag_has_been_set("--fourier_diff_modes"))
    {
      Fourier_diff_modes_n_angular_elements_ratio = -1.0; // Invalid ratio => Manual
    }

    Calculate_hankel_series_from_exact=command_line_flag_has_been_set("--calculate_hankel_series_from_exact");

    // Mapping options
    Use_newton_transformation=command_line_flag_has_been_set("--newton_mapping");
    Use_conformal_newton_transformation=command_line_flag_has_been_set("--conformal_newton_mapping");
    Debug_mapping=command_line_flag_has_been_set("--debug_mapping");
    Impose_straight_line_in_pml=command_line_flag_has_been_set("--impose_straight_line_in_pml");

    // Other options
    Fourier_diffmesh_full_integration=
      command_line_flag_has_been_set("--fourier_full_integration");
    Pml_full_integration=
      command_line_flag_has_been_set("--pml_integration_points");
    Disable_solve=command_line_flag_has_been_set("--disable_solve");
    Suppress_doc_solution=command_line_flag_has_been_set("--suppress_doc_solution");

    // Run specific tests (these will overide other settings if required)
    if (command_line_flag_has_been_set("--single_hankel"))
    {
      Use_single_hankel_solution=true;
    }

    if(command_line_flag_has_been_set("--test_single_mode"))
    {
      Use_single_mode_solution=true;
    }

    if(command_line_flag_has_been_set("--test_diff"))
    {
      Disable_solve = true;
      doc_fourier = true;
      Calculate_fourier_coefs_from_exact = true;
    }
    if(command_line_flag_has_been_set("--null_mapping_test"))
    {
      Null_mapping_test = true;
      Fill_in = NULL_MAPPING;
      Set_exact_on_outer = true;
    }
    if(command_line_flag_has_been_set("--dummy_mapping_test"))
    {
      Dummy_mapping_test = true;
      Set_exact_on_outer = true;
    }
    if(command_line_flag_has_been_set("--objective_test"))
    {
      Dummy_mapping_test = true;
      Disable_solve = true;
      Use_conformal_newton_transformation = true;
      Output_objective_samples = true;
    }

    // Doc what has actually been specified on the command line
    doc_specified_flags();

    // Update any global parameters which depend on ones we have just set
    update_parameters();
  }

  // Create label for output
  //------------------------
  DocInfo doc_info;

  // Set output directory
  doc_info.set_directory("RESLT");

  // Construct problem with specific number of nodes (order) and triangle/quad
#ifdef TRIANGLE
  if(GlobalParameters::N_node == 2)
  {
    PMLProblem<TPMLHelmholtzElement<2, 2,
      AnnularFromCartesianPMLElement> >  problem;
    problem.solve_and_doc(doc_info);
  }
  if(GlobalParameters::N_node == 3)
  {
    PMLProblem<TPMLHelmholtzElement<2, 3,
      AnnularFromCartesianPMLElement> >  problem;
    problem.solve_and_doc(doc_info);
  }
  if(GlobalParameters::N_node == 4)
  {
    PMLProblem<TPMLHelmholtzElement<2, 4,
      AnnularFromCartesianPMLElement> >  problem;
    problem.solve_and_doc(doc_info);
  }
#else
  if(GlobalParameters::N_node == 2)
  {
    PMLProblem<RefineableQPMLHelmholtzElement<2, 2,
      AnnularFromCartesianPMLElement> >  problem;
    problem.solve_and_doc(doc_info);
  }
  if(GlobalParameters::N_node == 3)
  {
    PMLProblem<RefineableQPMLHelmholtzElement<2, 3,
      AnnularFromCartesianPMLElement> >  problem;
    problem.solve_and_doc(doc_info);
  }
  if(GlobalParameters::N_node == 4)
  {
    PMLProblem<RefineableQPMLHelmholtzElement<2, 4,
      AnnularFromCartesianPMLElement> >  problem;
    problem.solve_and_doc(doc_info);
  }
#endif


  return 0;
}
