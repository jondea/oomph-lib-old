//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
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
#ifndef OOMPH_LAGRANGE_ENFORCED_FLOW_PRECONDITIONERS_HEADER
#define OOMPH_LAGRANGE_ENFORCED_FLOW_PRECONDITIONERS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// oomphlib headers
#include "../generic/matrices.h"
#include "../generic/assembly_handler.h"
#include "../generic/problem.h"
#include "../generic/block_preconditioner.h"
#include "../generic/preconditioner.h"
#include "../generic/SuperLU_preconditioner.h"
#include "../generic/matrix_vector_product.h"
#include "navier_stokes_elements.h"
#include "refineable_navier_stokes_elements.h"
#include "navier_stokes_preconditioners.h"

namespace oomph
{

namespace Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper
{
  /// \short CG with diagonal preconditioner for the Lagrange multiplier
  /// subsidiary linear systems.
  extern Preconditioner* get_lagrange_multiplier_preconditioner();
}


//==========================================================================
/// \short The preconditioner for the Lagrange multiplier constrained 
/// Navier-Stokes equations. The velocity components are constrained by 
/// Lagrange multiplier, which are applied via OOMPH-LIB's FACE elements.
/// 
/// A Vector of meshes is taken, each mesh contains a different type of
/// block preconditionable element. Each element must not only classify it's 
/// own degrees of freedom but also the associated DOF from the 'bulk' 
/// element.
/// 
/// The first mesh in the Vector Mesh_pt is assumed to be the 'bulk' mesh.
/// The rest are assumed to contain FACEELMENTS applying the required 
/// constraint.
/// 
/// Thus the most general block structure (in 3D) is: 
/// 
///  0 1 2 3   4 5 6 7  8  ..x   x+0 x+1 x+2 x+3 x+4 
/// [u v w p] [u v w l1 l2 ...] [u   v   w   l1  l2 ...] ... 
///   Bulk       Surface 1             Surface 2         ... 
/// 
/// where the dof types in [] are the dof types in each mesh.
/// It is assumed that in all surface mesh (after the bulk mesh), the first
/// spatial dimension number of dof types are the constrained velocity.
/// 
/// Consider the case of imposing parallel outflow (3 constrained velocity
/// dof types and 2 lagrange multiplier dof types) and tangential flow (3
/// constrained velocity dof types and 1 lagrange multiplier dof type)
/// along two different boundaries in 3D. The resulting natural block dof
/// type structure is: 
/// [0 1 2 3] [4  5  6   7   8 ] [9  10 11 12 ]
/// [u v w p] [up vp wp Lp1 Lp2] [ut vt wt Lt1]
/// 
/// Given that we know the spatial dimension of the problem, this information
/// can be conveniently stored in a Vector N_doftype_in_mesh = [4, 5, 4]. This
/// Vector will be used to re-order the dof types to group together the
/// velocity, pressure, then lagrange dof types like so: 
/// 
///  0 4  9  1 5  10  2 6  11    3    7   8  12   
/// [u up ut v vp vt  w wp wt ] [p] [Lp1 Lp2 Lt1] 
///
///    0 4  9  1 5  10  2 6  11  3  7   8  12   
///    u up ut v vp vt  w wp wt  p Lp1 Lp2 Lt1  
///  ..... this is too hard to do without Latex....
///
/// We use the preconditioner in the form... check my first year report...
/// 
/// Giving rise to the blocked Jacobian:
/// F G^t
///      L
/// D   
///  L
/// 
/// Here F is the momentum block, G the discrete gradient operator,
/// and D the discrete divergence operator. (For unstabilised elements, 
/// we have D = G^T and in much of the literature the divergence matrix is 
/// denoted by B.) The L blocks
//==========================================================================
class LagrangeEnforcedflowPreconditioner 
  : public BlockPreconditioner<CRDoubleMatrix>
{
  public:

  /// \short This preconditioner includes the option to use subsidiary 
  /// operators other than SuperLUPreconditioner for this problem. 
  /// This is the typedef of a function that should return an instance
  /// of a subsidiary preconditioning operator.  This preconditioner is 
  /// responsible for the destruction of the subsidiary preconditioners.
  typedef Preconditioner* (*SubsidiaryPreconditionerFctPt)();

  /// Constructor - sets the defaults for control flags
  LagrangeEnforcedflowPreconditioner():BlockPreconditioner<CRDoubleMatrix>()
  {
    // Null the pointers:
    
    // The Navier Stokes preconditioner pointer.
    Navier_stokes_preconditioner_pt = 0;
    
    // Null the function pointer which is used to create the subsidiary
    // preconditioners for the W block(s)
    Lagrange_multiplier_subsidiary_preconditioner_function_pt = 0;

    // Set the vector of preconditioner pointers for the W block(s) to zero.
    Lagrange_multiplier_preconditioner_pt.resize(0,0);

    // By default, the Lagrange multiplier block is preconditioned by the 
    // SuperLU preconditioner. To be honest, it does not matter what we set
    // here, it is reset later in the function setup(...).
    Lagrange_multiplier_preconditioner_is_block_preconditioner = false;

    // By default, the Navier Stokes block is preconditioned by the 
    // Least Square Commutator (LSC) preconditioner. To be honest, it foes not
    // matter what we set here, since it is reset in the function setup(...).
    Navier_stokes_preconditioner_is_block_preconditioner = true;

    // flag to indicate to use SuperLU or not.
    Using_superlu_ns_preconditioner = true;

    // Initially, there are no meshes set!
    My_mesh_pt.resize(0,0);

    // The number of meshes
    My_nmesh = 0;

    // The number of DOF types within the meshes.
    My_ndof_types_in_mesh.resize(0,0);

    // flag to indicate LSC preconditioner
    Use_default_norm_of_f_scaling = true;
    Scaling_sigma = 0.0;
    Scaling_sigma_multiplier = 1.0;
    N_lagrange_doftypes = 0;
    N_fluid_doftypes = 0;

    N_velocity_doftypes = 0;

    Doc_time = false;

    Doc_prec = false;

    Doc_linear_solver_info_pt = 0;

    Use_diagonal_w_block = true;

    Mapping_info_calculated = false;

    Label_pt = 0;

    Doc_prec_directory_pt = 0;

    Do_matcat_test = false;
    Do_vec_test = false;
  }

  /// destructor
  virtual ~LagrangeEnforcedflowPreconditioner()
  {
    this->clean_up_memory();
  }

  /// Broken copy constructor
  LagrangeEnforcedflowPreconditioner (const LagrangeEnforcedflowPreconditioner&)
  {
    BrokenCopy::broken_copy("LagrangeEnforcedflowPreconditioner");
  }

  /// Broken assignment operator
  void operator=(const LagrangeEnforcedflowPreconditioner&)
  {
    BrokenCopy::broken_assign(" LagrangeEnforcedflowPreconditioner");
  }

  /// Setup method for the LagrangeEnforcedflowPreconditioner.
  void setup();

  /// Add a scalar to each of the diagonal entry of a matrix.
  void add_scaling_to_diag(double &Scaling, CRDoubleMatrix *&block_pt);

  /// Extract the diagonal entries of a matrix.
  void get_diag(CRDoubleMatrix *&block_pt, Vector<double>& diag);

  /// Element-wise addition of two matrices.
  void add_matrices(CRDoubleMatrix *&block_pt1, CRDoubleMatrix *&block_pt2);

  /// Use the diagonal approximation for the W block.
  void use_diagonal_w_block() {Use_diagonal_w_block = true;}

  /// Use block diagonal W block.
  void use_block_diagonal_w_block() {Use_diagonal_w_block  = false;}

  ///Enable documentation of time
  void enable_doc_prec() {Doc_prec = true;}

  ///Disable documentation of time
  void disable_doc_prec() {Doc_prec = false;}

  ///Enable documentation of time
  void enable_doc_time() {Doc_time = true;}

  ///Disable documentation of time
  void disable_doc_time() {Doc_time = false;}

  /// Set  label
  void set_label_pt(std::string* label_pt) 
  {
    Label_pt = label_pt;
  }

  /// 
  void set_doc_prec_directory_pt(std::string* doc_prec_directory_pt) 
  {
    Doc_prec_directory_pt = doc_prec_directory_pt;
  }

  /// \short Apply the preconditioner.
  /// r is the residual (rhs), z will contain the solution.
  void preconditioner_solve(const DoubleVector& r, DoubleVector& z)
  {
    // Working vectors.
    DoubleVector temp_vec;
    DoubleVector another_temp_vec;
    DoubleVector yet_another_temp_vec;

    // First we solve all the w blocks:
    if(Lagrange_multiplier_preconditioner_is_block_preconditioner)
    {
      for (unsigned l_i = 0; l_i < N_lagrange_doftypes; l_i++) 
      {
        Lagrange_multiplier_preconditioner_pt[l_i]->preconditioner_solve(r,z);
      }
    }
    else
    {

      // Loop through all of the Lagrange multipliers
      for(unsigned l_i = 0; l_i < N_lagrange_doftypes; l_i++)
      {
        // Get the block type of block l_i
        const unsigned l_ii = N_fluid_doftypes + l_i;

        // Extract the block
        this->get_block_vector(l_ii,r,temp_vec);

        Lagrange_multiplier_preconditioner_pt[l_i]
          ->preconditioner_solve(temp_vec,another_temp_vec);

        const unsigned vec_nrow_local = another_temp_vec.nrow_local();
        double* vec_values_pt = another_temp_vec.values_pt();
        
        for (unsigned i = 0; i < vec_nrow_local; i++) 
        {
          vec_values_pt[i] = vec_values_pt[i]*Scaling_sigma;
        }
        this->return_block_vector(l_ii,another_temp_vec,z);
        
        temp_vec.clear();
        another_temp_vec.clear();
      }
    }

    // Now solve the Navier-Stokes block.

    // At this point, all vectors are cleared.
    if(Using_superlu_ns_preconditioner)
    {
      // Get the concatenated fluid vector.
      Vector<unsigned> fluid_block_indices(N_fluid_doftypes,0);
      for (unsigned b = 0; b < N_fluid_doftypes; b++) 
      {
        fluid_block_indices[b] = b;
      }

      this->get_concatenated_block_vector(fluid_block_indices,r,temp_vec);

      // temp_vec contains the (concatenated) fluid rhs.
      Navier_stokes_preconditioner_pt
        ->preconditioner_solve(temp_vec,another_temp_vec);

      temp_vec.clear();

      // Now return it.
      this->return_concatenated_block_vector(fluid_block_indices,
                                             another_temp_vec,z);

      another_temp_vec.clear();
    }
    else
    {
      // This is a BlockPreconditioner
      Navier_stokes_preconditioner_pt->preconditioner_solve(r,z);
    }
  } // end of preconditioner_solve

  /// Set the meshes, the first mesh must be the fluid mesh
  void set_meshes(const Vector<Mesh*> &mesh_pt)
  {
    // There should be at least two meshes for this preconditioner.
    const unsigned nmesh = mesh_pt.size();

#ifdef PARANOID
    if(nmesh < 2)
    {
      std::ostringstream err_msg;
      err_msg << "There should be at least two meshes.\n";
      throw OomphLibError(err_msg.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Check that all pointers are not null
    for(unsigned mesh_i = 0; mesh_i < nmesh; mesh_i++)
    {
      if (mesh_pt[mesh_i]==0)
      {
        std::ostringstream err_msg;
        err_msg << "The pointer mesh_pt[" << mesh_i << "] is null.\n";
        throw OomphLibError(err_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }

    // We assume that the first mesh is the Navier Stokes "bulk" mesh. 
    // To check this, the elemental dimension must be the same as the 
    // nodal (spatial) dimension.
    //
    // We store the elemental dimension i.e. the number of local coordinates
    // required to parametrise its geometry.
    const unsigned elemental_dim = mesh_pt[0]->elemental_dimension();

    // The dimension of the nodes in the first element in the (supposedly)
    // bulk mesh.
    const unsigned nodal_dim = mesh_pt[0]->nodal_dimension();

    // Check if the first mesh is the "bulk" mesh.
    // Here we assume only one mesh contains "bulk" elements.
    // All subsequent meshes contain block preconditionable elements which
    // re-classify the bulk velocity dofs to constrained velocity dofs.
    if (elemental_dim != nodal_dim) 
    {
      std::ostringstream err_msg;
      err_msg << "In the first mesh, the elements have elemental dimension of "
              << elemental_dim << ", with a nodal dimension of "
              << nodal_dim << ".\n"
              << "The first mesh does not contain 'bulk' elements.\n"
              << "Please re-order your mesh_pt vector.\n";

      throw OomphLibError(err_msg.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    
    // Set the number of meshes 
    this->set_nmesh(nmesh);
    
    // Set the meshes
    for(unsigned mesh_i = 0; mesh_i < nmesh; mesh_i++)
     {
      this->set_mesh(mesh_i,mesh_pt[mesh_i]);
     }

    // We also store the meshes and number of meshes locally in this class.
    // This is slightly redundant, since we always have it stored in the upper
    // most master block preconditioner. But at the moment there is no 
    // mapping/look up scheme between master and subsidiary block 
    // preconditioners for meshes.
    // So if this is a subsidiary block preconditioner, we don't know which of
    // the master's meshes belong to us. We need this information to set up
    // look up lists in the function setup(...).
    // So we store them local to this class.
    My_mesh_pt = mesh_pt;
    My_nmesh = nmesh;
  } // EoFunc set_meshes


  /// \short Access function to the Scaling sigma of the preconditioner
  double& scaling_sigma()
  {
    Use_default_norm_of_f_scaling = false;
    return Scaling_sigma;
  }
  
  /// \short Function to get the scaling Sigma of the preconditioner
  double scaling_sigma() const
  {
    return Scaling_sigma;
  }

  /// \short Access function to the Scaling sigma of the preconditioner
  double& scaling_sigma_multiplier()
  {
    return Scaling_sigma_multiplier;
  }

  /// \short Function to get the scaling Sigma of the preconditioner
  double scaling_sigma_multiplier() const
  {
    return Scaling_sigma_multiplier;
  } 

  /// Use default scaling?
  void use_default_norm_of_f_scaling()
  {
    Use_default_norm_of_f_scaling = true;
  }

  /// Function to set a new momentum matrix preconditioner (inexact solver)
  void set_navier_stokes_lsc_preconditioner(
      Preconditioner* new_ns_preconditioner_pt = 0)
  {
    // Check if pointer is non-zero.
    if(new_ns_preconditioner_pt == 0)
    {
      std::ostringstream warning_stream;
      warning_stream << "WARNING: \n"
                     << "The LSC preconditioner point is null.\n" 
                     << "Using default (SuperLU) preconditioner.\n" 
                     << std::endl;
      OomphLibWarning(warning_stream.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
      
      Navier_stokes_preconditioner_pt = 0;
      Using_superlu_ns_preconditioner = true;
    }
    else
    {
      // If the default SuperLU preconditioner has been used
      // clean it up now...
      if (Using_superlu_ns_preconditioner 
          && Navier_stokes_preconditioner_pt != 0)
      {
        delete Navier_stokes_preconditioner_pt;
      }
      
      Navier_stokes_preconditioner_pt = new_ns_preconditioner_pt;
      Using_superlu_ns_preconditioner = false;
    }
  }

  ///\short Function to (re-)set momentum matrix preconditioner (inexact
  /// solver) to SuperLU
  void set_superlu_preconditioner_for_navier_stokes_block()
  {
    if (!Using_superlu_ns_preconditioner)
    {
      delete Navier_stokes_preconditioner_pt;
      Navier_stokes_preconditioner_pt = new SuperLUPreconditioner;
      Using_superlu_ns_preconditioner = true;
    }
  }

  /// Maybe I need to remove this.
  void set_doc_linear_solver_info_pt
    (DocLinearSolverInfo* doc_linear_solver_info_pt)
  {
    Doc_linear_solver_info_pt = doc_linear_solver_info_pt;
  }

  /// \short By default the Lagrange multiplier subsidiary systems are 
  /// preconditioner with SuperLUPreconditioner. For a different 
  /// preconditioner, pass a function to this 
  /// method returning a different subsidiary operator.
  void set_lagrange_multiplier_subsidiary_preconditioner(
      SubsidiaryPreconditionerFctPt prec_fn)
  {
    Lagrange_multiplier_subsidiary_preconditioner_function_pt = prec_fn;
  }

  /// \short Clears the memory.
  void clean_up_memory();

  bool Do_matcat_test;
  bool Do_vec_test;

  private:

  /// 
  std::string* Label_pt;

  /// 
  std::string *Doc_prec_directory_pt;

  /// 
  bool Doc_time;

  /// 
  bool Doc_prec;

  /// 
  bool Mapping_info_calculated;

  /// \short the Scaling_sigma variable of this preconditioner
  double Scaling_sigma;
  double Scaling_sigma_multiplier;

  /// 
  bool Use_default_norm_of_f_scaling;

  /// The Lagrange multiplier subsidiary preconditioner function pointer
  SubsidiaryPreconditionerFctPt 
    Lagrange_multiplier_subsidiary_preconditioner_function_pt;

  /// Same W solvers are used in both exact and LSC.
  /// Pointer to the 'preconditioner' for the W matrix
  Vector<Preconditioner*> Lagrange_multiplier_preconditioner_pt;

  /// Boolean to indicate if the Lagrange multiplier preconditioner is a 
  /// block preconditioner or not.
  bool Lagrange_multiplier_preconditioner_is_block_preconditioner;

  /// Pointer to the 'preconditioner' for the Navier-Stokes block
  Preconditioner* Navier_stokes_preconditioner_pt;

  /// Boolean to indicate if the preconditioner for the Navier Stokes block
  /// is a block preconditioner or not.
  bool Navier_stokes_preconditioner_is_block_preconditioner;

  /// flag to indicate whether the default NS preconditioner is used
  bool Using_superlu_ns_preconditioner;

  /// Bool to use diagonal or block diagonal W block.
  bool Use_diagonal_w_block;

  /// the re-arraned doftypes: velocity, pressure, lagrange.
  Vector<unsigned> Doftype_list_vpl;
  
  /// the re-arraned doftypes: bulk, constrained, pressure, lagrange.
  Vector<unsigned> Subsidiary_list_bcpl;

  /// \short Storage for the meshes. 
  /// The first mesh must always be the Navier-Stokes (bulk) mesh.
  Vector<Mesh*> My_mesh_pt;

  /// 
  Vector<unsigned> My_ndof_types_in_mesh;

  /// \short Store the number of meshes. 
  /// This is required to create the look up lists.
  unsigned My_nmesh;

  /// \short The number of lagrange multiplier dof types
  unsigned N_lagrange_doftypes;

  /// \short The number of fluid dof types (including pressure)
  unsigned N_fluid_doftypes;

  /// \short The number of velocity dof types.
  unsigned N_velocity_doftypes;

  bool Replace_all_f_blocks;



  /// \short Pointer to Doc_linear_solver_info.
  /// used for book keeping purposes.
  DocLinearSolverInfo* Doc_linear_solver_info_pt;

}; // end of LagrangeEnforcedflowPreconditioner class




/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////









  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////



// void LagrangeEnforcedflowPreconditioner::assemble_inv_press_and_veloc_mass_matrix_diagonal
}
#endif
