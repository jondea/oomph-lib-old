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
#include <cmath>

#include <functional>

// Generic routines
#include "generic.h"

// The Helmholtz equations
#include "helmholtz.h"

// Boundary ids for triangle mesh
#define FISH_BOUNDARY_ID1 0
#define FISH_BOUNDARY_ID2 1
#define FISH_BOUNDARY_ID3 2
#define FISH_BOUNDARY_ID4 3
#define FISH_BOUNDARY_ID5 4
#define FISH_BOUNDARY_ID6 5
#define FISH_BOUNDARY_ID7 6
#define FISH_BOUNDARY_ID8 7
#define OUTER_BULK_BOUNDARY_ID1 8
#define OUTER_BULK_BOUNDARY_ID2 9

// Outer PML boundary ID
#define OUTER_PML_BOUNDARY 2

// The pml Helmholtz equations
#include "pml_helmholtz.h"

// The meshes needed in the PML constructions
#include "meshes/triangle_mesh.h"
#include "meshes/rectangular_quadmesh.h"

// Get the Bessel functions
#include "oomph_crbond_bessel.h"

#include "meshes/wrap_pml.h"

// #include "mesh_utils.h"

using namespace oomph;
using namespace MathematicalConstants;

//=========================================================================
///    2D Geometric Object in the shape of a simple fish
//        ______
//       /_     \/|
//         \      |
//        _/      |
//       \______/\|
//=========================================================================
class FishGeomObject : public GeomObject
{
public:

  /// \short Constructor, default fish is roughly the size of the unit circle
  /// use scale parameter to make uniformly bigger and smaller
  FishGeomObject(const double scale = 1.0, const bool puffer_mode=false):
    GeomObject(1,2), Puffer_mode(puffer_mode)
  {
    Mouth_depth_ratio = scale*0.5;
    Body_height = scale*1.0;
    Body_width = scale*1.5;
    Tail_height = scale*1.0;
    Tail_length = scale*0.15;
    Theta_lip = 2.7;
    Theta_tail = 0.6;
  }

  /// Broken copy constructor
  FishGeomObject(const FishGeomObject& dummy)
  {
    BrokenCopy::broken_copy("FishGeomObject");
  }

  /// Broken assignment operator
  void operator=(const FishGeomObject&)
  {
    BrokenCopy::broken_assign("FishGeomObject");
  }

  /// Destructor:  Empty
  ~FishGeomObject(){}

  /// \short Position Vector at Lagrangian coordinate zeta
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
    if(std::abs(zeta[0]) > 4.0)
    {
      oomph_info << "zeta is" << zeta[0] << std::endl;
      throw OomphLibError("Zeta must be between -4.0 and 4.0.",
        OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    if(Puffer_mode)
    {
      r[0] = Body_height*cos(Pi*zeta[0]/4.0);
      r[1] = Body_height*sin(Pi*zeta[0]/4.0);
      return;
    }

    // Positive zeta is top half, negative zeta is bottom half
    // Starts at mouth and goes round to end of tail
         if(std::abs(zeta[0]) <= 1.0) mouth(zeta, r);
    else if(std::abs(zeta[0]) <= 2.0) body(zeta, r);
    else if(std::abs(zeta[0]) <= 3.0) tail_slope(zeta, r);
    else if(std::abs(zeta[0]) <= 4.0) tail_back(zeta, r);
  }

  /// Mouth of fish (straight line)
  void mouth(const Vector<double>& zeta, Vector<double>& r) const
  {
    double mouth_hinge_x = -(1.0-Mouth_depth_ratio)*Body_width/2.0;
    double mouth_hinge_y = 0.0;
    double lip_x = Body_width/2.0*cos(Theta_lip);
    double lip_y = Body_height/2.0*sin(Theta_lip);

    double t = (std::abs(zeta[0]) - 0.0);

    r[0] = mouth_hinge_x + (lip_x-mouth_hinge_x)*t;
    r[1] = mouth_hinge_y + (lip_y-mouth_hinge_y)*t;

    // Negative z means lower, so flip y coord
    r[1] = copysign(r[1], zeta[0]);
  }

  double mouth_arclength() const
  {
    double mouth_hinge_x = -(1.0-Mouth_depth_ratio)*Body_width/2.0;
    double mouth_hinge_y = 0.0;
    double lip_x = Body_width/2.0*cos(Theta_lip);
    double lip_y = Body_height/2.0*sin(Theta_lip);

    return sqrt(   (mouth_hinge_x-lip_x)*(mouth_hinge_x-lip_x)
                 + (mouth_hinge_y-lip_y)*(mouth_hinge_y-lip_y) );
  }

  // Body of fish (circular arc)
  void body(const Vector<double>& zeta, Vector<double>& r) const
  {
    double theta = (std::abs(zeta[0]) - 1.0)*(Theta_tail - Theta_lip) + Theta_lip;
    r[0] =  Body_width/2.0*cos(theta);
    r[1] =  Body_height/2.0*sin(theta);

    // Negative z means lower, so flip y coord
    r[1] = copysign(r[1], zeta[0]);
  }

  double body_arclength() const
  {
    // Approximation of perimeter of half an ellipse by Ramanujan
    double a = Body_width /2.0;
    double b = Body_height/2.0;
    return (Theta_lip-Theta_tail)/2.0*(3.0*(a+b)-sqrt((3.0*a+b)*(a+3.0*b)));
  }

  /// Slope of the tail of the fish (straight line)
  void tail_slope(const Vector<double>& zeta, Vector<double>& r) const
  {
    double start_x = Body_width /2.0*cos(Theta_tail);
    double start_y = Body_height/2.0*sin(Theta_tail);
    double end_x = Body_width/2.0 + Tail_length;
    double end_y = Tail_height/2.0;

    double t = (std::abs(zeta[0]) - 2.0);

    r[0] = start_x + (end_x - start_x)*t;
    r[1] = start_y + (end_y - start_y)*t;

    // Negative z means lower, so flip y coord
    r[1] = copysign(r[1], zeta[0]);
  }

  double tail_slope_arclength() const
  {
    double start_x = Body_width /2.0*cos(Theta_tail);
    double start_y = Body_height/2.0*sin(Theta_tail);
    double end_x = Body_width/2.0 + Tail_length;
    double end_y = Tail_height/2.0;

    return sqrt(   (start_x-end_x)*(start_x-end_x)
                 + (start_y-end_y)*(start_y-end_y) );
  }

  /// Back of the tail of the fish (straight line)
  void tail_back(const Vector<double>& zeta, Vector<double>& r) const
  {
    double start_x = Body_width/2.0 + Tail_length;
    double start_y = Tail_height/2.0;
    double end_x = start_x;
    double end_y = 0.0;

    double t = (std::abs(zeta[0]) - 3.0);

    r[0] = start_x + (end_x - start_x)*t;
    r[1] = start_y + (end_y - start_y)*t;

    // Negative z means lower, so flip y coord
    r[1] = copysign(r[1], zeta[0]);
  }

  double tail_back_arclength() const
  {
    double start_x = Body_width/2.0 + Tail_length;
    double start_y = Tail_height/2.0;
    double end_x = start_x;
    double end_y = 0.0;

    return sqrt(   (start_x-end_x)*(start_x-end_x)
                 + (start_y-end_y)*(start_y-end_y) );
  }

  double arclength(const double& z1, const double& z2) const
  {
    if(Puffer_mode) return Pi*std::abs(z2 - z1)/4.0;

    // Loop over each section, increasing z and adding the contributions into a
    double a = 0.0;
    double z = z1;
    while(z < z2)
    {
      if(z < -3.0)
      {
        a += (min(-3.0, z2) - z)*this->tail_back_arclength();
        z = -3.0;
      }
      else if(z < -2.0)
      {
        a += (min(-2.0, z2) - z)*this->tail_slope_arclength();
        z = -2.0;
      }
      else if(z < -1.0)
      {
        a += (min(-1.0, z2) - z)*this->body_arclength();
        z = -1.0;
      }
      else if(z < 0.0)
      {
        a += (min(0.0, z2) - z)*this->mouth_arclength();
        z = 0.0;
      }
      else if(z < 1.0)
      {
        a += (min(1.0, z2) - z)*this->mouth_arclength();
        z = 1.0;
      }
      else if(z < 2.0)
      {
        a += (min(2.0, z2) - z)*this->body_arclength();
        z = 2.0;
      }
      else if(z < 3.0)
      {
        a += (min(3.0, z2) - z)*this->tail_slope_arclength();
        z = 3.0;
      }
      else if(z < 4.0)
      {
        a += (min(4.0, z2) - z)*this->tail_back_arclength();
        z = 4.0;
      }
    }
    return a;
  }

private:

  /// How far mouth goes into fish, 0 is no mouth, 1 takes mouth to origin
  double Mouth_depth_ratio;

  /// Height (top to bottom in y) of the body of the fish
  double Body_height;

  /// Width (left to right in x)  of the body of the fish
  double Body_width;

  /// Height (top to bottom in y) of the tail of the fish
  double Tail_height;

  /// Distance the tail sticks out past the circle (part of which forms body)
  double Tail_length;

  /// Angle the lip (end of mouth) makes the with origin
  double Theta_lip;

  /// Angle the intersection of the tail and body makes the with origin
  double Theta_tail;

  /// Turn into circle
  bool Puffer_mode;

};


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

/// Angle of the incoming planar wave which is scattered off the fish
double Incoming_wave_angle = 0.9;

/// Amplitude of the incoming planar wave which is scattered off the fish
double Incoming_wave_amplitude = 1.0;

/// Incoming wave solution
std::complex<double> Incoming_wave(Vector<double>& x)
{
  double k_x = cos(Incoming_wave_angle)*Wavenumber;
  double k_y = sin(Incoming_wave_angle)*Wavenumber;
  return Incoming_wave_amplitude*exp(I*(k_x*x[0] + k_y*x[1]));
}

/// Size of fish, roughly the inner radius of domain
double Fish_scale = 1.0;

/// Turn the fish into a circle (like a puffer fish)
bool Puffer_mode = true;

/// Outer radius of bulk, inner radius of PML
double Bulk_outer_r = 2.0;

/// Outer radius of PML
double Pml_outer_r = 3.0;

/// Number of nodes per 1D element, order of elements + 1
unsigned N_node = 2;

/// Rough target element length used to generate triangles
double Element_length = 0.1;

/// Target element area to feed to triangle, assumes an equilateral triangle
double Element_area = (sqrt(3.0)/4.0)*Element_length*Element_length;

/// Maximum number of times to adapt solution before stopping
unsigned Max_adapt = 0;

/// Use a DtN mapping to impos outgoing BCs
bool Use_dtn = false;

// Number of PML elements through the thickness of the PML
unsigned N_pml = 1;

/// Use a the scale free Bermudez PML Mapping, if false we use the default
bool Use_sfb = false;

/// Should we disable Newton solve (just set up problem and doc solution)
bool Disable_solve = false;

/// Should we suppress output of solution?
bool Short_doc = false;

bool Dump_restart_file = false;

bool Calculate_error_from_restart_file = false;

MeshAsGeomObject* Exact_soln_mesh_pt = 0;

template<class ELEMENT>
void get_u_from_mesh(const Vector<double>& x, Vector<double>& u)
{
  GeomObject* geom_pt;
  Vector<double> s(2,0.0);
  GlobalParameters::Exact_soln_mesh_pt->locate_zeta(x, geom_pt, s);

  ELEMENT* el_pt = dynamic_cast<ELEMENT*>(geom_pt);

  //Find number of nodes
  const unsigned n_node = el_pt->nnode();

  //Local shape function
  Shape psi(n_node);

  //Find values of shape function
  el_pt->shape(s,psi);

  //Get the index at which the helmholtz unknown is stored
  const unsigned u_nodal_index_real = el_pt->u_index_helmholtz().real();
  const unsigned u_nodal_index_imag = el_pt->u_index_helmholtz().imag();

  //Loop over the local nodes and sum
  u[0] = 0.0;
  u[1] = 0.0;
  for(unsigned l=0;l<n_node;l++)
  {
    //Add to the interpolated value
    u[0] += el_pt->nodal_value(l,u_nodal_index_real)*psi[l];
    u[1] += el_pt->nodal_value(l,u_nodal_index_imag)*psi[l];
  }
}

UniaxialPMLMapping* Sfb_mapping_pt = new ScaleFreeBermudezPMLMapping;

void update_parameters()
{
  // Calculate k^2 from k and print to console
  K_squared = Wavenumber*Wavenumber;
  oomph_info << "k^2 = " << K_squared << std::endl;

  // Calculate target element area from element length and print to console
  Element_area = (sqrt(3.0)/4.0)*Element_length*Element_length;
  oomph_info << "Target element area = " << Element_area << std::endl;

}

} // end of namespace

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//========= start_of_problem_class=====================================
/// Problem class to demonstrate use of perfectly matched layers
/// for Helmholtz problems.
//=====================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
class PMLProblem : public Problem
{

public:

/// Constructor
PMLProblem(const bool use_dtn, const bool enable_doc = true);

/// Destructor (empty)
~PMLProblem(){
};

/// Update the problem specs before solve (empty)
void actions_before_newton_solve(){
  if (this->Use_dtn)
  {
    Dtn_mesh_pt->setup_gamma();
  }
};

/// Update the problem specs after solve (empty)
void actions_after_newton_solve(){
};

/// Recompute gamma integral before checking Newton residuals
void actions_before_newton_convergence_check()
{
  if (this->Use_dtn)
  {
    Dtn_mesh_pt->setup_gamma();
  }
}

/// Actions before adapt: Wipe the PML meshes
void actions_before_adapt();

/// Actions after adapt: Rebuild the PML meshes
void actions_after_adapt();

/// Wipe the PML/DtN meshes
void actions_before_read_unstructured_meshes()
{
  this->actions_before_adapt();
};

/// Do nothing because we don't use the DtN/PML meshes
void actions_after_read_unstructured_meshes(){};

/// \short Doc the solution
void doc_solution();

/// \short Create Helmholtz flux elements on boundary b of the Mesh pointed
/// to by mesh_pt and add them to the specified surface Mesh
void create_flux_elements(const unsigned &b, Mesh* const &mesh_pt,
                          Mesh* const & helmholtz_inner_boundary_mesh_pt);

void apply_inner_boundary_conditions();

void set_k_squared_in_bulk();

void create_and_setup_pmls();

/// \short Create BC elements on boundary b of the Mesh pointed
/// to by bulk_mesh_pt and add them to the specified survace Mesh
void create_dtn_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
  Mesh* const & dtn_mesh_pt);

void create_and_setup_dtn();

// Load coefficients, do a newton solve then doc solution
void solve_and_doc();

/// Pointer to the "bulk" mesh
RefineableTriangleMesh<BULK_ELEMENT>* Bulk_mesh_pt;

private:

/// \short Pointer to the mesh containing element in PML
Mesh* Pml_mesh_pt;

/// \short Pointer to mesh containing the DtN boundary condition elements
HelmholtzDtNMesh<BULK_ELEMENT>* Dtn_mesh_pt;

/// Trace file
ofstream Trace_file;

/// DocInfo object stores flags/labels for where the / output gets written to
DocInfo Doc_info;

/// Flag whether or not this problem will use a DtN mapping, this is member data
/// so that we can construct a separate problem for exact soln
bool Use_dtn;

/// Flag which decides whether we run the doc_solution function
bool Enable_doc;

}; // end of problem class

//=======start_of_constructor=============================================
/// Constructor for Helmholtz problem
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
PMLProblem<BULK_ELEMENT, PML_ELEMENT>::PMLProblem(const bool use_dtn,
  const bool enable_doc):
    Pml_mesh_pt(0), Dtn_mesh_pt(0), Doc_info(), Use_dtn(use_dtn),
    Enable_doc(enable_doc)
{
  // Set output directory
  Doc_info.set_directory("RESLT");

  // Open trace file
  Trace_file.open("RESLT/trace.dat");

  // Define the centre of all the circles (the origin)
  Vector<double> centre(2,0.0);

  // Fishy
  //---------------
  // Create circle representing inner boundary
  FishGeomObject* fish_pt=new FishGeomObject(GlobalParameters::Fish_scale,
                                             GlobalParameters::Puffer_mode);

  Vector<TriangleMeshCurveSection*> fish_section_pt(8);

  // The intrinsic coordinates for the beginning and end of the section
  double s_start = -4.0;
  double s_end = s_start + 1.0;

  unsigned n_segments;

  // Loop over the 8 sections of the fish
  unsigned b_start = FISH_BOUNDARY_ID1;
  unsigned i = 0;
  for(unsigned b=b_start; b<b_start+8; b++)
  {

    // Number of segments on eighth of fish
    n_segments = ceil(fish_pt->arclength(s_start,s_end)
                      /GlobalParameters::Element_length);
    //
    oomph_info << "Construcing fish segment number " << b << " with "
               << n_segments << " segments"  << std::endl;
    fish_section_pt[i] =
      new TriangleMeshCurviLine(fish_pt, s_start, s_end, n_segments, b);
    s_start = s_end;
    s_end += 1.0;
    i++;
  }

  // Combine fish sections
  TriangleMeshClosedCurve* inner_boundary_pt =
    new TriangleMeshClosedCurve(fish_section_pt, centre);


  // Outer boundary
  //---------------
  // Create circle representing outer boundary
  Circle* circle_pt=new Circle(centre[0], centre[1],
                               GlobalParameters::Bulk_outer_r);

  Vector<TriangleMeshCurveSection*> boundary_line_pt(2);

  // Number of segments on half of circle
  n_segments = max(2.0, ceil(Pi*GlobalParameters::Bulk_outer_r
                              /GlobalParameters::Element_length));

  // The intrinsic coordinates for the beginning and end of the curve
  s_start = 0.0;
  s_end   = Pi;
  boundary_line_pt[0]=
    new TriangleMeshCurviLine(circle_pt, s_start, s_end, n_segments,
                              OUTER_BULK_BOUNDARY_ID1);

  // The intrinsic coordinates for the beginning and end of the curve
  s_start = Pi;
  s_end   = 2.0*Pi;
  boundary_line_pt[1]=
    new TriangleMeshCurviLine(circle_pt, s_start, s_end, n_segments,
                              OUTER_BULK_BOUNDARY_ID2);

  // Combine two semi-circles and return
  TriangleMeshClosedCurve* outer_boundary_pt =
    new TriangleMeshClosedCurve(boundary_line_pt,centre);

  // Create mesh
  // -----------------------
  // Use the TriangleMeshParameters object for helping on the manage
  // of the TriangleMesh parameters. The only parameter that needs to take
  // is the outer boundary.
  TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

  // Specify the internal closed boundaries
  Vector<TriangleMeshClosedCurve*> bulk_boundary_pt_vec(1);
  bulk_boundary_pt_vec[0] = inner_boundary_pt;
  triangle_mesh_parameters.internal_closed_curve_pt() = bulk_boundary_pt_vec;

  // Target element size in bulk mesh (make it massive so that it just uses
  // segments)
  triangle_mesh_parameters.element_area() = GlobalParameters::Element_area;

  // Build "bulk" mesh
  Bulk_mesh_pt=new RefineableTriangleMesh<BULK_ELEMENT>(triangle_mesh_parameters);

  // Disable projection as problem is linear (so only wastes time)
  Bulk_mesh_pt->disable_projection();

  // Tiny minimum element size when refining to pick up sharp geometry
  Bulk_mesh_pt->min_element_size() = 1.0e-6;

  // Create the main mesh
  add_sub_mesh(Bulk_mesh_pt);

  // Create/set error estimator
  Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

  // Set k squared pointer in the bulk
  set_k_squared_in_bulk();

  if(this->Use_dtn)
  {
    // Create DtN mesh and elements
    create_and_setup_dtn();
  }
  else
  {
    // Create PML mesh, create and set mapping and set k squared
    create_and_setup_pmls();
  }

  // Build the entire mesh from its submeshes
  build_global_mesh();

  // Boundary conditions
  apply_inner_boundary_conditions();

  // Setup equation numbering scheme
  cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;

  // Now that equation numbers have been set up, compute gamma integral
  if(this->Use_dtn) Dtn_mesh_pt->setup_gamma();

} // end of constructor


//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of face elements
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::actions_before_adapt()
{
  // Save results to file
  if(Enable_doc) doc_solution();

  // Remove PML meshes because we do not adapt them
  if(Pml_mesh_pt != 0)
  {
    delete Pml_mesh_pt;
    Pml_mesh_pt=0;
  }

  // Remove DtN mesh because we do not adapt them
  if(Dtn_mesh_pt != 0)
  {
    delete Dtn_mesh_pt;
    Dtn_mesh_pt=0;
  }

  // Rebuild the Problem's global mesh from its various sub-meshes
  // but first flush all its submeshes
  flush_sub_meshes();

  // Then add the triangular mesh back
  add_sub_mesh(Bulk_mesh_pt);

  //  Rebuild the global mesh such that it now stores
  // the triangular mesh only
  rebuild_global_mesh();

}// end of actions_before_adapt


//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the face element meshes
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::actions_after_adapt()
{
  // Set k squared in all bulk elements
  set_k_squared_in_bulk();

  if(this->Use_dtn)
  {
    // Create DtN mesh and elements
    create_and_setup_dtn();
  }
  else
  {
    // Create PML mesh, create and set mapping and set k squared
    create_and_setup_pmls();
  }

  // Apply boundary conditions on inner bulk
  apply_inner_boundary_conditions();

  // Rebuild the Problem's global mesh from its various sub-meshes
  rebuild_global_mesh();

}// end of actions_after_adapt

//=======================start_of_set_k_squared_in_bulk=========================
/// Loop over bulk elements and set k_squared pointer
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::set_k_squared_in_bulk()

{
  // Complete the build of all elements so they are fully functional
  unsigned n_element = Bulk_mesh_pt->nelement();
  for(unsigned e=0; e<n_element; e++)
  {
    // Upcast from GeneralisedElement to Pml Helmholtz bulk element
    BULK_ELEMENT *el_pt =
      dynamic_cast<BULK_ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    if (el_pt!=0)
    {
      //Set the k_squared function pointer
      el_pt->k_squared_pt() = &GlobalParameters::K_squared;
    }
  }
}

//====================start_of_create_and_setup_pmls======================
/// Create wrapped PML around outer boundary, and set up with wavenumber
/// and mapping pointer and pin outer PML boundary at zero
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::create_and_setup_pmls()
{

  // Use helper function to build a wrapped PML mesh around outer bulk boundary
  Vector<unsigned> Outer_boundary_ids(2);
  Outer_boundary_ids[0] = OUTER_BULK_BOUNDARY_ID1;
  Outer_boundary_ids[1] = OUTER_BULK_BOUNDARY_ID2;
  Pml_mesh_pt =
    WrapPMLHelper::create_wrapped_pml_mesh<BULK_ELEMENT, PML_ELEMENT>(
      Bulk_mesh_pt, Outer_boundary_ids, GlobalParameters::N_pml);

  // Add newly created PML mesh to global mesh
  add_sub_mesh(Pml_mesh_pt);

  // Loop over all elements, setting k_squared and mapping
  unsigned n_element = Pml_mesh_pt->nelement();
  for(unsigned e=0; e<n_element; e++)
  {
    // Upcast from GeneralisedElement to Pml Helmholtz bulk element
    PML_ELEMENT *el_pt =
      dynamic_cast<PML_ELEMENT*>(Pml_mesh_pt->element_pt(e));

    if (el_pt!=0)
    {
      //Set the k_squared function pointer
      el_pt->k_squared_pt() = &GlobalParameters::K_squared;

      el_pt->pml_mapping_pt() = GlobalParameters::Sfb_mapping_pt;
    }
  }

  // Pin outer boundary nodes of PML at zero
  unsigned n_node = Pml_mesh_pt->nboundary_node(OUTER_PML_BOUNDARY);
  for (unsigned n=0; n<n_node; n++)
  {
    // Get pointer to node
    Node* nod_pt=Pml_mesh_pt->boundary_node_pt(OUTER_PML_BOUNDARY, n);

    // Pin all dofs
    nod_pt->pin(0);
    nod_pt->pin(1);

    // Assign the value
    nod_pt->set_value(0,0.0);
    nod_pt->set_value(1,0.0);
  }

}

//====================start_of_create_and_setup_dtn=======================
/// Create DtN mesh and add elements around the boundary, then setup
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::create_and_setup_dtn()
{
  // Choose number of fourier modes to use based on number of boundary elements
  unsigned n_el_outer_boundaries = 0;
  n_el_outer_boundaries += Bulk_mesh_pt->nboundary_element(OUTER_BULK_BOUNDARY_ID1);
  n_el_outer_boundaries += Bulk_mesh_pt->nboundary_element(OUTER_BULK_BOUNDARY_ID2);
  unsigned n_fourier = std::min(20u, n_el_outer_boundaries/4);

  // Pointer to mesh containing the Helmholtz outer boundary condition
  // elements. Specify outer radius and number of Fourier terms to be
  // used in gamma integral
  Dtn_mesh_pt = new HelmholtzDtNMesh<BULK_ELEMENT>(
    GlobalParameters::Bulk_outer_r, n_fourier);

  // Create outer boundary elements from all elements that are
  // adjacent to the outer boundary , but add them to a separate mesh.
  create_dtn_elements(OUTER_BULK_BOUNDARY_ID1, Bulk_mesh_pt, Dtn_mesh_pt);
  create_dtn_elements(OUTER_BULK_BOUNDARY_ID2, Bulk_mesh_pt, Dtn_mesh_pt);

  // Add Dtn mesh to global mesh
  add_sub_mesh(Dtn_mesh_pt);

  // Loop over the flux elements to pass pointer to DtN
  // BC for the outer boundary
  unsigned n_element=Dtn_mesh_pt->nelement();
  for(unsigned e=0;e<n_element;e++)
  {
    // Upcast from GeneralisedElement to Helmholtz flux element
    HelmholtzDtNBoundaryElement<BULK_ELEMENT> *el_pt =
      dynamic_cast< HelmholtzDtNBoundaryElement<BULK_ELEMENT>*>(
        Dtn_mesh_pt->element_pt(e));

    // Set pointer to the mesh that contains all the boundary condition
    // elements on this boundary
    el_pt->set_outer_boundary_mesh_pt(Dtn_mesh_pt);
  }

}

//============start_of_create_dtn_elements==============================
/// Create outer DtN elements on the b-th boundary of
/// the Mesh object pointed to by bulk_mesh_pt and add the elements
/// to the Mesh object pointed to by dtn_mesh_pt.
//===========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::
create_dtn_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                    Mesh* const & dtn_mesh_pt)
{
  // Loop over the bulk elements adjacent to boundary b?
  unsigned n_element = bulk_mesh_pt->nboundary_element(b);
  for(unsigned e=0;e<n_element;e++)
  {
    // Get pointer to the bulk element that is adjacent to boundary b
    BULK_ELEMENT* bulk_elem_pt = dynamic_cast<BULK_ELEMENT*>(
      bulk_mesh_pt->boundary_element_pt(b,e));

    //Find the index of the face of element e along boundary b
    int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

    // Build the corresponding outer flux element
    HelmholtzDtNBoundaryElement<BULK_ELEMENT>* flux_element_pt = new
      HelmholtzDtNBoundaryElement<BULK_ELEMENT>(bulk_elem_pt,face_index);

    //Add the flux boundary element to the  helmholtz_outer_boundary_mesh
    dtn_mesh_pt->add_element_pt(flux_element_pt);
  }
} // end of create_dtn_elements


//===============start_of_apply_inner_boundary_conditions=================
/// Apply boundary conditions
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::apply_inner_boundary_conditions()
{

  // Boundary conditions are set on the surface of the circle
  // as a constant nonzero Dirichlet boundary condition
  unsigned n_bound = Bulk_mesh_pt->nboundary();

  for(unsigned b=0; b<n_bound; b++)
  {

    // Set dirichlet exact boundary conditions on inner boundary
    if (b == FISH_BOUNDARY_ID1 || b == FISH_BOUNDARY_ID2 ||
        b == FISH_BOUNDARY_ID3 || b == FISH_BOUNDARY_ID4 ||
        b == FISH_BOUNDARY_ID5 || b == FISH_BOUNDARY_ID6 ||
        b == FISH_BOUNDARY_ID7 || b == FISH_BOUNDARY_ID8)
    {
      unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
      for (unsigned n=0; n<n_node; n++)
      {
        // Get pointer to node
        Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

        // Extract nodal coordinates from node:
        Vector<double> pos(2);
        nod_pt->position(pos);

        // Compute the value of the exact solution at the nodal point
        std::complex<double> u = - GlobalParameters::Incoming_wave(pos);

        nod_pt->pin(0);
        nod_pt->pin(1);

        // Assign the value
        nod_pt->set_value(0,u.real());
        nod_pt->set_value(1,u.imag());

      }
    }
  }
} // end of apply_inner_boundary_conditions

//=====================start_of_doc=======================================
/// Doc the solution: Doc_info contains labels/output directory etc.
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::doc_solution()
{

  ofstream file;
  char filename[100];

  // Number of plot points for fine and coarse solns
  unsigned npts=GlobalParameters::N_node;
  unsigned coarse_npts=2;

  if(!GlobalParameters::Short_doc)
  {
    // Output solution
    //-----------------
    sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
            Doc_info.number());
    file.open(filename);
    Bulk_mesh_pt->output(file,npts);
    file.close();

    // Output coarse solution
    //-----------------
    sprintf(filename,"%s/coarse_soln%i.dat",Doc_info.directory().c_str(),
            Doc_info.number());
    file.open(filename);
    Bulk_mesh_pt->output(file,coarse_npts);
    file.close();

    if (!this->Use_dtn)
    {
      // Output PML solution
      //-----------------
      sprintf(filename,"%s/pml_soln%i.dat",Doc_info.directory().c_str(),
      Doc_info.number());
      file.open(filename);
      Pml_mesh_pt->output(file,npts);
      file.close();
    }

  }

  // Output number of bulk nodes
  //----------------------------
  sprintf(filename,"%s/bulk_nodes%i.dat",Doc_info.directory().c_str(),
          Doc_info.number());
  file.open(filename);
  file << Bulk_mesh_pt->nnode() << std::endl;
  file.close();

  if(GlobalParameters::Calculate_error_from_restart_file)
  {
    // Construct problem, then read mesh/data from restart file
    bool exact_use_dtn = true;
    bool exact_enable_doc = false; // Disable doc_solution to avoid recursion
    PMLProblem<BULK_ELEMENT, PML_ELEMENT>* Exact_soln_problem_pt =
      new PMLProblem<BULK_ELEMENT, PML_ELEMENT>(exact_use_dtn, exact_enable_doc);
    std::ifstream exact_restart_file("exact_soln_restart_file.dat");
    bool unsteady = false;
    Exact_soln_problem_pt->actions_before_read_unstructured_meshes();
    Exact_soln_problem_pt->read(exact_restart_file, unsteady);
    exact_restart_file.close();

    GlobalParameters::Exact_soln_mesh_pt = new MeshAsGeomObject(Exact_soln_problem_pt->Bulk_mesh_pt);

    // Calculate error
    double error,norm;
    sprintf(filename,"%s/restart_error%i.dat",Doc_info.directory().c_str(),
    Doc_info.number());
    file.open(filename);
    Bulk_mesh_pt->compute_error(file,GlobalParameters::get_u_from_mesh<BULK_ELEMENT>,error,norm);
    file.close();

    // Now we are done we can delete the MeshAsGeomObject and the Problem
    delete GlobalParameters::Exact_soln_mesh_pt;
    GlobalParameters::Exact_soln_mesh_pt = 0;
    delete Exact_soln_problem_pt;

    // Doc L2 error and norm of solution
    oomph_info << std::endl << "Error calculated using restart file" << std::endl;
    oomph_info << "Norm of error   : " << sqrt(error) << std::endl;
    oomph_info << "Norm of solution: " << sqrt(norm) << std::endl;
    oomph_info << "Relative error: " << sqrt(error/norm) << std::endl
               << std::endl;

    // Write relative error to a file for easy access
    sprintf(filename,"%s/relative_error%i.dat",Doc_info.directory().c_str(),
      Doc_info.number());
    file.open(filename);
    file << sqrt(error/norm) << std::endl;
    file.close();
    
    Trace_file << sqrt(norm) << ", "<< sqrt(error/norm) << std::endl;
  }

  if(GlobalParameters::Dump_restart_file)
  {
    sprintf(filename,"%s/soln_restart_file%i.dat",Doc_info.directory().c_str(),
      Doc_info.number());
    file.open(filename);
    file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    this->Enable_doc = false;
    this->actions_before_adapt();
    this->Enable_doc = true;
    this->dump(file);
    file.close();
  }

  Doc_info.number()++;
} // end of doc


template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::solve_and_doc()
{
  // Problem is actually linear, no need to recompute residuals before and after
  // problem level Newton solve.
  this->problem_is_nonlinear(true);

  // Solve the problem
  if(!GlobalParameters::Disable_solve) newton_solve(GlobalParameters::Max_adapt);

  // Save results to file
  if(Enable_doc) doc_solution();
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
    specify_command_line_flag("--fish_scale",&Fish_scale);
    specify_command_line_flag("--bulk_outer_r",&Bulk_outer_r);
    specify_command_line_flag("--outer_r",&Pml_outer_r);
    specify_command_line_flag("--el_length",&Element_length);
    specify_command_line_flag("--n_pml",&N_pml);
    specify_command_line_flag("--use_sfb");
    specify_command_line_flag("--n_node",&N_node);
    specify_command_line_flag("--puffer_mode");
    specify_command_line_flag("--max_adapt",&Max_adapt);

    // Other options
    specify_command_line_flag("--wavenumber",&Wavenumber);
    specify_command_line_flag("--disable_solve");
    specify_command_line_flag("--short_doc");
    specify_command_line_flag("--use_dtn");
    specify_command_line_flag("--calculate_error_from_file");
    specify_command_line_flag("--dump_problem");

    // Parse command line
    parse_and_assign();

    Puffer_mode=command_line_flag_has_been_set("--puffer_mode");

    // Other options
    Disable_solve=command_line_flag_has_been_set("--disable_solve");
    Short_doc=command_line_flag_has_been_set("--short_doc");
    Use_dtn=command_line_flag_has_been_set("--use_dtn");
    Use_sfb=command_line_flag_has_been_set("--use_sfb");
    Calculate_error_from_restart_file=command_line_flag_has_been_set("--calculate_error_from_file");
    Dump_restart_file=command_line_flag_has_been_set("--dump_problem");

    // Doc what has actually been specified on the command line
    doc_specified_flags();

    // Update any global parameters which depend on ones we have just set
    update_parameters();
  }

  // Construct problem with specific number of nodes (order) and triangle/quad
  if(GlobalParameters::N_node == 2)
  {
    PMLProblem< ProjectableHelmholtzElement< THelmholtzElement<2,2> >,
      QPMLHelmholtzElement<2,2,Conformal2DPMLElement> >
        problem(GlobalParameters::Use_dtn);
    problem.solve_and_doc();
  }
  if(GlobalParameters::N_node == 3)
  {
    PMLProblem< ProjectableHelmholtzElement< THelmholtzElement<2,3> >,
      QPMLHelmholtzElement<2,3,Conformal2DPMLElement> >
        problem(GlobalParameters::Use_dtn);
    problem.solve_and_doc();
  }
  if(GlobalParameters::N_node == 4)
  {
    PMLProblem< ProjectableHelmholtzElement< THelmholtzElement<2,4> >,
      QPMLHelmholtzElement<2,4,Conformal2DPMLElement> >
        problem(GlobalParameters::Use_dtn);
    problem.solve_and_doc();
  }


  return 0;
}
