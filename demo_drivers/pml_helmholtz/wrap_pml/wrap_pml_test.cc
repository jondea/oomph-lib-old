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
// Driver for a specific 2D Helmholtz problem with
// perfectly matched layer treatment for the exterior boundaries

#include <fenv.h>

#include "math.h"
#include <complex>

// Generic routines
#include "generic.h"

// The pml Helmholtz equations
#include "pml_helmholtz.h"

// The meshes needed in the PML constructions
#include "meshes/triangle_mesh.h"
#include "meshes/rectangular_quadmesh.h"

#include "meshes/wrap_pml.h"

using namespace oomph;
using namespace std;

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//===== start_of_namespace=============================================
/// Namespace for the Helmholtz problem parameters
//=====================================================================
namespace GlobalParameters
{

/// Wavenumber (also known as k), k=omega/c
double Wavenumber = sqrt(50.0);

/// Square of the wavenumber (also known as k^2)
double K_squared = Wavenumber * Wavenumber;

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
PMLProblem();

/// Destructor (empty)
~PMLProblem(){
}

/// Update the problem specs before solve (empty)
void actions_before_newton_solve(){
}

/// Update the problem specs after solve (empty)
void actions_after_newton_solve(){
}

/// \short Doc the solution. DocInfo object stores flags/labels for where the
/// output gets written to
void doc_solution(DocInfo& doc_info);

/// Create PML meshes
void create_pml_meshes();

// Apply boundary conditions
void apply_boundary_conditions();

#ifdef ADAPTIVE

/// Actions before adapt: Wipe the PML meshes
void actions_before_adapt();

/// Actions after adapt: Rebuild the PML meshes
void actions_after_adapt();

#endif

private:

#ifdef ADAPTIVE

/// Pointer to the refineable "bulk" mesh
RefineableTriangleMesh<BULK_ELEMENT>* Bulk_mesh_pt;

#else

/// Pointer to the "bulk" mesh
TriangleMesh<BULK_ELEMENT>* Bulk_mesh_pt;

#endif

/// Ids of the outer boundaries
Vector<unsigned> Outer_boundary_id;

/// Pointer to the PML mesh
Mesh* Pml_mesh_pt;

/// Trace file
ofstream Trace_file;

}; // end of problem class

//=========================================================================
/// Straight 1D line in 2D space
//=========================================================================
class MyStraightLine : public GeomObject
{
public:
  /// Constructor:  Pass start and end points
  MyStraightLine(const Vector<double>& r_start, const Vector<double>& r_end,
    const double& zeta_start, const double& zeta_end)
   :  GeomObject(1,2), R_start(r_start), R_end(r_end), Zeta_start(zeta_start),
      Zeta_end(zeta_end)
   { }


  /// Broken copy constructor
  MyStraightLine(const MyStraightLine& dummy)
   {
    BrokenCopy::broken_copy("MyStraightLine");
   }

  /// Broken assignment operator
  void operator=(const MyStraightLine&)
   {
    BrokenCopy::broken_assign("MyStraightLine");
   }


  /// Destructor:  Empty
  ~MyStraightLine(){}

  /// \short Position Vector at Lagrangian coordinate zeta
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
    // Fraction along the line
    double fraction = (zeta[0]-Zeta_start)/(Zeta_end - Zeta_start);

    // Position Vector
    r[0] = R_start[0]+(R_end[0]-R_start[0])*fraction;
    r[1] = R_start[1]+(R_end[1]-R_start[1])*fraction;
  }

private:

  /// Start point of line
  Vector<double> R_start;

  /// End point of line
  Vector<double> R_end;

  /// Start zeta value
  double Zeta_start;

  /// End zeta value
  double Zeta_end;

};



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////





//=======start_of_constructor=============================================
/// Constructor for Helmholtz problem
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
PMLProblem<BULK_ELEMENT, PML_ELEMENT>::PMLProblem(): Outer_boundary_id(4)
{

  // Open trace file
  Trace_file.open("RESLT/trace.dat");

  // Create circle representing inner boundary
  double a=0.2;
  double x_c=0.0;
  double y_c=0.0;
  Circle* inner_circle_pt=new Circle(x_c, y_c, a);

  // Outer boundary
  //---------------
  TriangleMeshClosedCurve* outer_boundary_pt=0;

  Vector<TriangleMeshCurveSection*> outer_boundary_line_pt(4);
  Vector<MyStraightLine*> outer_boundary_geom_line_pt(4);

  // Each polyline only has three vertices, provide storage for their
  // coordinates
  Vector<double> start_vertex(2);
  Vector<double> end_vertex(2);

  // First polyline
  start_vertex[0]=-2.0;
  start_vertex[1]=-2.0;
  end_vertex[0]=2.0;
  end_vertex[1]=-2.0;

  // Build the 1st boundary line
  Outer_boundary_id[0]=2;
  double s_start = 0.0;
  double s_end   = 1.0;
  outer_boundary_geom_line_pt[0] = new MyStraightLine(start_vertex, end_vertex,
                                                      s_start, s_end);
  outer_boundary_line_pt[0]=
    new TriangleMeshCurviLine(outer_boundary_geom_line_pt[0], s_start, s_end, 1,
                              Outer_boundary_id[0]);

  // Second boundary polyline
  start_vertex[0]=2.0;
  start_vertex[1]=-2.0;
  end_vertex[0]=2.0;
  end_vertex[1]=2.0;

  // Build the 2nd boundary polyline
  Outer_boundary_id[1]=3;
  s_start = s_end;
  s_end   = s_start + 1.0;
  outer_boundary_geom_line_pt[1] = new MyStraightLine(start_vertex, end_vertex,
                                                      s_start, s_end);
  outer_boundary_line_pt[1]=
    new TriangleMeshCurviLine(outer_boundary_geom_line_pt[1], s_start, s_end, 1,
                              Outer_boundary_id[1]);

  // Third boundary polyline
  start_vertex[0]=2.0;
  start_vertex[1]=2.0;
  end_vertex[0]=-2.0;
  end_vertex[1]=2.0;

  // Build the 3rd boundary polyline
  Outer_boundary_id[2]=4;
  s_start = s_end;
  s_end   = s_start + 1.0;
  outer_boundary_geom_line_pt[2] = new MyStraightLine(start_vertex, end_vertex,
                                                      s_start, s_end);
  outer_boundary_line_pt[2]=
    new TriangleMeshCurviLine(outer_boundary_geom_line_pt[2], s_start, s_end, 1,
                              Outer_boundary_id[2]);

  // Fourth boundary polyline
  start_vertex[0]=-2.0;
  start_vertex[1]=2.0;
  end_vertex[0]=-2.0;
  end_vertex[1]=-2.0;

  // Build the 4th boundary polyline
  Outer_boundary_id[3]=5;
  s_start = s_end;
  s_end   = s_start + 1.0;
  outer_boundary_geom_line_pt[3] = new MyStraightLine(start_vertex, end_vertex,
                                                      s_start, s_end);
  outer_boundary_line_pt[3]=
    new TriangleMeshCurviLine(outer_boundary_geom_line_pt[3], s_start, s_end, 1,
                              Outer_boundary_id[3]);


  // Create the triangle mesh closed curve for outer boundary
  Vector<double> inside_coords(2);
  inside_coords[0]=1.0;
  inside_coords[1]=1.0;
  outer_boundary_pt=new TriangleMeshClosedCurve(outer_boundary_line_pt,
                                                inside_coords);

  // Inner circular boundary
  //------------------------
  unsigned n_segments = 8;
  Vector<TriangleMeshCurveSection*> inner_boundary_line_pt(2);

  // The intrinsic coordinates for the beginning and end of the curve
  s_start = 0.0;
  s_end   = MathematicalConstants::Pi;
  unsigned boundary_id = 0;
  inner_boundary_line_pt[0]=
    new TriangleMeshCurviLine(inner_circle_pt,
                              s_start,
                              s_end,
                              n_segments,
                              boundary_id);

  // The intrinsic coordinates for the beginning and end of the curve
  s_start = MathematicalConstants::Pi;
  s_end   = 2.0*MathematicalConstants::Pi;
  boundary_id = 1;
  inner_boundary_line_pt[1]=
    new TriangleMeshCurviLine(inner_circle_pt,
                              s_start,
                              s_end,
                              n_segments,
                              boundary_id);


  // Combine to hole
  //----------------
  Vector<TriangleMeshClosedCurve*> hole_pt(1);
  Vector<double> hole_coords(2);
  hole_coords[0]=0.0;
  hole_coords[1]=0.0;
  hole_pt[0]=new TriangleMeshClosedCurve(inner_boundary_line_pt,hole_coords);


  // Use the TriangleMeshParameters object for helping on the manage
  // of the TriangleMesh parameters. The only parameter that needs to take
  // is the outer boundary.
  TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

  // Specify the closed curve using the TriangleMeshParameters object
  triangle_mesh_parameters.internal_closed_curve_pt() = hole_pt;

  // Target element size in bulk mesh
  triangle_mesh_parameters.element_area() = 0.5;

#ifdef ADAPTIVE

  // Build adaptive "bulk" mesh
  Bulk_mesh_pt=new RefineableTriangleMesh<BULK_ELEMENT>(triangle_mesh_parameters);

  // Create/set error estimator
  Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

  // Choose error tolerances to force some uniform refinement
  Bulk_mesh_pt->min_permitted_error()=0.00004;
  Bulk_mesh_pt->max_permitted_error()=0.0001;

#else

  // Build "bulk" mesh
  Bulk_mesh_pt=new TriangleMesh<BULK_ELEMENT>(triangle_mesh_parameters);

#endif

  // Create the main triangular mesh
  add_sub_mesh(Bulk_mesh_pt);

  // Create PML meshes and add them to the global mesh
  create_pml_meshes();

  // Build the entire mesh from its submeshes
  build_global_mesh();

  // Let's have a look where the boundaries are
  Bulk_mesh_pt->output("global_mesh.dat");
  Bulk_mesh_pt->output_boundaries("global_mesh_boundary.dat");

  ofstream outfile("global_mesh_boundary_coordinates.dat");
  for (unsigned i = 0; i < 4; i++)
  {
    Bulk_mesh_pt->output_boundary_coordinates(i,outfile);
  }
  outfile.close();

  // Complete the build of all elements so they are fully functional
  unsigned n_element = Bulk_mesh_pt->nelement();
  for(unsigned e=0; e<n_element; e++)
  {
    // Upcast from GeneralisedElement to Helmholtz bulk element
    PMLHelmholtzEquationsBase<2> *el_pt =
      dynamic_cast<PMLHelmholtzEquationsBase<2>*>(Bulk_mesh_pt->element_pt(e));

    //Set the k_squared double pointer
    el_pt->k_squared_pt() = &GlobalParameters::K_squared;
  }

  // Apply boundary conditions
  apply_boundary_conditions();

  // Setup equation numbering scheme
  cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor

#ifdef ADAPTIVE

//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of face elements
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::actions_before_adapt()
{
  // Before adapting the added PML meshes must be removed
  // as they are not refineable and are to be rebuilt from the
  // newly refined triangular mesh
  delete PML_right_mesh_pt;
  PML_right_mesh_pt=0;

  // Rebuild the Problem's global mesh from its various sub-meshes
  // but first flush all its submeshes
  flush_sub_meshes();

  // Then add the triangular mesh back
  add_sub_mesh(Bulk_mesh_pt);

  //  Rebuild the global mesh such that it now stores
  // the triangular mesh only
  rebuild_global_mesh();

} // end of actions_before_adapt


//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the face element meshes
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::actions_after_adapt()
{

  // Build PML meshes  and add them to the global mesh
  create_pml_meshes();

  // Build the entire mesh from its submeshes
  rebuild_global_mesh();

  // Complete the build of all elements so they are fully functional

  // Loop over the entire mesh elements to set up element-specific
  // things that cannot be handled by constructor
  unsigned n_element = Bulk_mesh_pt->nelement();

  for(unsigned e=0; e<n_element; e++)
  {
    // Upcast from GeneralisedElement to PMLHelmholtz bulk element
    PMLHelmholtzEquationsBase<2> *el_pt =
      dynamic_cast<PMLHelmholtzEquationsBase<2>*>(Bulk_mesh_pt->element_pt(e));

    //Set the frequency function pointer
    el_pt->k_squared_pt() = &GlobalParameters::K_squared;
  }

  // Re-apply boundary conditions
  apply_boundary_conditions();

} // end of actions_after_adapt

#endif

//============start_of_create_pml_meshes======================================
/// Create PML meshes and add them to the problem's sub-meshes
//============================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::create_pml_meshes()
{
  // Use helper function to build a wrapped pml mesh
  Pml_mesh_pt =
    WrapPMLHelper::create_wrapped_pml_mesh<BULK_ELEMENT,PML_ELEMENT>(
      Bulk_mesh_pt, Outer_boundary_id, 10);

  // Add to problem as sub mesh
  add_sub_mesh(Pml_mesh_pt);

} // end of create_pml_meshes


//==================start_of_apply_boundary_conditions====================
/// Apply boundary conditions
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::apply_boundary_conditions()
{

  // Boundary conditions are set on the surface of the circle
  // as a constant nonzero Dirichlet boundary condition
  unsigned n_bound = Bulk_mesh_pt->nboundary();

  for(unsigned b=0; b<n_bound; b++)
  {
    unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
    for (unsigned n=0; n<n_node; n++)
    {
      if ((b==0) || (b==1))
      {
        Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,n);
        nod_pt->pin(0);
        nod_pt->pin(1);

        nod_pt->set_value(0,0.1);
        nod_pt->set_value(1,0.0);
      }
    }
  }

} // end of apply_boundary_conditions


//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class BULK_ELEMENT, class PML_ELEMENT>
void PMLProblem<BULK_ELEMENT, PML_ELEMENT>::doc_solution(DocInfo& doc_info)
{

  ofstream some_file,some_file2;
  char filename[100];

  // Number of plot points
  unsigned npts = 3;
  unsigned npts_coarse=2;

  // Output solution
  //-----------------
  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file,npts);
  some_file.close();

  // Output coarse solution
  //-----------------------
  sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file,npts_coarse);
  some_file.close();


  // Output solution within pml domains
  //-----------------------------------
  sprintf(filename,"%s/pml_soln%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Pml_mesh_pt->output(some_file,npts);
  some_file.close();

  // Output coarse solution within pml domains
  //------------------------------------------
  sprintf(filename,"%s/pml_coarse_soln%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Pml_mesh_pt->output(some_file,npts_coarse);
  some_file.close();


  // Write norm of solution to trace file
  double norm=0.0;
  Bulk_mesh_pt->compute_norm(norm);
  Trace_file << norm << std::endl;

  // //Do animation of Helmholtz solution
  // //-----------------------------------
  // unsigned nstep=40;
  // for (unsigned i=0;i<nstep;i++)
  //  {
  //   sprintf(filename,"%s/helmholtz_animation%i_frame%i.dat",
  //           doc_info.directory().c_str(),
  //           doc_info.number(),i);
  //   some_file.open(filename);
  //   double phi=2.0*MathematicalConstants::Pi*double(i)/double(nstep-1);
  //   unsigned nelem=Bulk_mesh_pt->nelement();
  //   for (unsigned e=0;e<nelem;e++)
  //    {
  //     BULK_ELEMENT* el_pt=dynamic_cast<BULK_ELEMENT*>(
  //      Bulk_mesh_pt->element_pt(e));
  //     el_pt->output_real(some_file,phi,npts);
  //    }
  //   some_file.close();
  //  }

} // end of doc

//==========start_of_main=================================================
/// Solve 2D Helmholtz problem
//========================================================================
int main(int argc, char **argv)
{
  CommandLineArgs::setup(argc,argv);

  // Don't complain if an element is inverted
  FiniteElement::Accept_negative_jacobian = true;

  //Set up the problem
  //------------------

#ifdef ADAPTIVE

  // Set up the problem with projectable 2D six-node elements from the
  // TPMLHelmholtzElement family.
  PMLProblem<ProjectablePMLHelmholtzElement
             <TPMLHelmholtzElement<2,3,Conformal2DPMLElement> >,
             QPMLHelmholtzElement<2,3,Conformal2DPMLElement> > problem;

  // Set up the problem with 2D ten-node elements from the
  // TPMLHelmholtzElement family.
  // PMLProblem<ProjectablePMLHelmholtzElement
  //  <TPMLHelmholtzElement<2,4> > > problem;

#else

  // Set up the problem with 2D six-node elements from the
  // TPMLHelmholtzElement family.
  PMLProblem<TPMLHelmholtzElement<2,3,Conformal2DPMLElement>,
             QPMLHelmholtzElement<2,3,Conformal2DPMLElement> > problem;

  // Set up the problem with 2D ten-node elements from the
  // TPMLHelmholtzElement family.
  //   PMLProblem<TPMLHelmholtzElement<2,4> >  problem;

#endif

  // Create label for output
  //------------------------
  DocInfo doc_info;

  // Set output directory
  doc_info.set_directory("RESLT");


#ifdef ADAPTIVE

  // Max. number of adaptations
  unsigned max_adapt=1;

  // Solve the problem with the adaptive Newton's method, allowing
  // up to max_adapt mesh adaptations after every solve.
  // problem.newton_solve(max_adapt);

#else

  // Solve the problem with Newton's method
  // problem.newton_solve();

#endif

  //Output solution
  problem.doc_solution(doc_info);

} //end of main
