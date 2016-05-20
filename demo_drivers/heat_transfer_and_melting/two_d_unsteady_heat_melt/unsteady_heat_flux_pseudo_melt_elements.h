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
// Header file for elements that are used to apply prescribed flux
// boundary conditions to the UnsteadyHeat equations, ensuring
// that boundary temperature never exceeds the melting temperature.
// The required flux (applied in addition to whatever is externally
// imposed) is just "sucked into nothing" here -- in the actual
// melt elements, this flux will be what drives the melt process.
#ifndef OOMPH_UNSTEADY_HEAT_FLUX_PSEUDO_MELT_ELEMENTS_HEADER
#define OOMPH_UNSTEADY_HEAT_FLUX_PSEUDO_MELT_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif


//Standard libray headers
#include <cmath>

// oomph-lib includes
#include "generic/Qelements.h"


namespace oomph
{

//======================================================================
/// \short A class for elements that allow the imposition of an 
/// applied flux on the boundaries of UnsteadyHeat elements.
/// "Temperature" is limited by melt-like construction --
/// just the right amount of flux is "withdrawn" (if required) 
/// to limit the temperature to the melt temperature.
/// The element geometry is obtained from the  FaceGeometry<ELEMENT> 
/// policy class. 
/// 
/// This is demo code -- the real "melt elements" which actually
/// move the surface according to the Stefan condition are implemented
/// elsewhere using the same functionality.
//======================================================================
template <class ELEMENT>
class UnsteadyHeatFluxPseudoMeltElement : public virtual FaceGeometry<ELEMENT>, 
public virtual FaceElement
{
 
public:

 /// \short Function pointer to the prescribed-flux function fct(x,f(x)) -- 
 /// x is a Vector! 
 typedef void (*UnsteadyHeatPrescribedFluxFctPt)
  (const double& time, const Vector<double>& x, double& flux);

 /// \short Constructor, takes the pointer to the "bulk" element and the 
 /// index of the face to be created, and (optionally) the face element ID.
 UnsteadyHeatFluxPseudoMeltElement(FiniteElement* const &bulk_el_pt, 
                             const int &face_index, const unsigned &id=0);
 
 /// Broken copy constructor
 UnsteadyHeatFluxPseudoMeltElement(const 
                                   UnsteadyHeatFluxPseudoMeltElement& dummy) 
  { 
   BrokenCopy::broken_copy("UnsteadyHeatFluxPseudoMeltElement");
  } 
 
 /// Access function for the prescribed-flux function pointer
 UnsteadyHeatPrescribedFluxFctPt& flux_fct_pt() {return Flux_fct_pt;}
 
 /// Compute the element residual vector
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
 {
  //Call the generic residuals function with flag set to 0
  //using a dummy matrix argument
  fill_in_generic_residual_contribution_ust_heat_flux(
   residuals,GeneralisedElement::Dummy_matrix,0);
 }
 
 /// Specify the value of nodal zeta from the face geometry
 /// \short The "global" intrinsic coordinate of the element when
 /// viewed as part of a geometric object should be given by
 /// the FaceElement representation, by default (needed to break
 /// indeterminacy if bulk element is SolidElement)
 double zeta_nodal(const unsigned &n, const unsigned &k,           
                   const unsigned &i) const 
 {return FaceElement::zeta_nodal(n,k,i);}     

 /// Output function
 void output(std::ostream &outfile)
 {
  unsigned nplot=5;
  output(outfile,nplot);
 }
 
 /// \short Output function
 void output(std::ostream &outfile, const unsigned &n_plot)
 {
  // Locally cache the index at which the variable is stored
  const unsigned u_index_ust_heat = U_index_ust_heat;
  
  // Get continuous time from timestepper of first node
  double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
    
  //Find out how many nodes there are
  const unsigned n_node = nnode();
  
  //Set up memory for the shape functions
  Shape psi(n_node);
  
  // Local and global coordinates
  Vector<double> s(Dim-1);
  Vector<double> interpolated_x(Dim);

  // Tecplot header info
  outfile << this->tecplot_zone_string(n_plot);
  
  // Loop over plot points
  unsigned num_plot_points=this->nplot_points(n_plot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {
    // Get local coordinates of plot point
    this->get_s_plot(iplot,n_plot,s);
    
    // Outer unit normal
    Vector<double> unit_normal(Dim);
    outer_unit_normal(s,unit_normal);
    
    //Find the shape functions
    shape(s,psi);
    
   // Melt flux
   double melt_flux=0.0;
   
   // Temperature
   double u=0.0;
   
   //Initialise to zero
   for(unsigned i=0;i<Dim;i++)
    {
     interpolated_x[i] = 0.0;
    }
   
   //Calculate stuff
   for(unsigned l=0;l<n_node;l++) 
    {
     // get the node pt
     Node* nod_pt = node_pt(l);
     
     // Cast to a boundary node
     BoundaryNodeBase *bnod_pt = 
      dynamic_cast<BoundaryNodeBase*>(nod_pt);
     
     // Get the index of the first nodal value associated with
     // this FaceElement
     unsigned first_index=
      bnod_pt->index_of_first_value_assigned_by_face_element(Melt_id);
     
     // Melt rate
     melt_flux+=nod_pt->value(first_index)*psi[l];
          
     // Temperature
     u+=nod_pt->value(u_index_ust_heat)*psi[l]; 
     
     //Loop over directions
     for(unsigned i=0;i<Dim;i++)
      {
       interpolated_x[i] += nodal_position(l,i)*psi[l];
      }
    }
   
   //Get the imposed flux
   double flux;
   get_flux(time,interpolated_x,flux);
   
   //Output the x,y,..
   for(unsigned i=0;i<Dim;i++) 
    {
     outfile << interpolated_x[i] << " ";
    }
   
   // Output imposed flux and melt flux
   outfile << flux << " ";
   outfile << melt_flux << " ";

   // Output normal
   for(unsigned i=0;i<Dim;i++) 
    {
     outfile << unit_normal[i] << " ";
    } 
    outfile << std::endl;
  }
 }

 /// C-style output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(FILE* file_pt) {FaceGeometry<ELEMENT>::output(file_pt);}

 /// \short C-style output function -- forward to broken version in 
 /// FiniteElement until somebody decides what exactly they want to plot 
 /// here...
 void output(FILE* file_pt, const unsigned &n_plot)
  {FaceGeometry<ELEMENT>::output(file_pt,n_plot);}



 /// \short Plot landscape of residuals for j-th node
 void plot_residual_landscape(std::ostream &landscape_outfile,
                              std::ostream &soln_outfile,
                              Problem* problem_pt,
                              const unsigned& j)
 {
  // Which node?
  Node* nod_pt=this->node_pt(j);
   
  // Cast to a boundary node
  BoundaryNodeBase *bnod_pt =
   dynamic_cast<BoundaryNodeBase*>(nod_pt);
  
  // Get the index of the first nodal value associated with
  // this FaceElement
  unsigned first_index=
   bnod_pt->index_of_first_value_assigned_by_face_element(Melt_id);
  
  // Backup values
  double u_back=nod_pt->value(U_index_ust_heat);
  double m_back=nod_pt->value(first_index);

  int global_eqn_u = nod_pt->eqn_number(U_index_ust_heat);
  int global_eqn_m = nod_pt->eqn_number(first_index);

  if (global_eqn_u<0) abort();
  if (global_eqn_m<0) abort();

  // Output actual solution
  DoubleVector residuals;
  problem_pt->get_residuals(residuals);
  soln_outfile << u_back << " "
               << m_back << " "
               << residuals[global_eqn_u] << " "
               << residuals[global_eqn_m]<< std::endl;

  // Range of plot values
  double u_min=melt_temperature()-5.0;
  double u_max=melt_temperature()+5.0;
  double m_min=-5.0;
  double m_max=5.0;

  unsigned nplot=100;
  landscape_outfile << "ZONE I=" << nplot << ", J=" << nplot << std::endl;
  for (unsigned i=0;i<nplot;i++)
   {
    double x0=u_min+(u_max-u_min)*double(i)/double(nplot-1);
    for (unsigned j=0;j<nplot;j++)
     {
      double x1=m_min+(m_max-m_min)*double(j)/double(nplot-1);

      nod_pt->set_value(U_index_ust_heat,x0);
      nod_pt->set_value(first_index,x1);

      problem_pt->get_residuals(residuals);
      landscape_outfile << x0 << " "
                        << x1 << " "
                        << residuals[global_eqn_u] << " "
                        << residuals[global_eqn_m] << std::endl;
     }
   }

  // Reset values
  nod_pt->set_value(U_index_ust_heat,u_back);
  nod_pt->set_value(first_index,m_back);
 }


 /// Melt temperature
 double melt_temperature()
 {
  if (Melt_temperature_pt==0)
   {
    return 0.0;
   }
  else
   {
    return *Melt_temperature_pt;
   }
 }

 /// Pointer to (non-default) melt temperature
 double*& melt_temperature_pt()
 {
  return Melt_temperature_pt;
 }
  
 
 /// Switch off melting by pinning melting dofs
 void disable_melting()
 {
  unsigned n_node=nnode();
  for(unsigned l=0;l<n_node;l++) 
   {
    // get the node pt
    Node* nod_pt = node_pt(l);
    
    // Cast to a boundary node
    BoundaryNodeBase *bnod_pt = 
     dynamic_cast<BoundaryNodeBase*>(nod_pt);
    
    // Get the index of the first nodal value associated with
    // this FaceElement
    unsigned first_index=
     bnod_pt->index_of_first_value_assigned_by_face_element(Melt_id);
    
    // Pin
    nod_pt->pin(first_index);
   }
 }
  
  protected:
 
 /// \short Function to compute the shape and test functions and to return 
 /// the Jacobian of mapping between local and global (Eulerian)
 /// coordinates
 inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
  const
 {
  //Find number of nodes
  unsigned n_node = nnode();
  
  //Get the shape functions
  shape(s,psi);
  
  //Set the test functions to be the same as the shape functions
  for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}
  
  //Return the value of the jacobian
  return J_eulerian(s);
 }
 
 /// \short Function to calculate the prescribed flux at a given spatial
 /// position and at a given time
 void get_flux(const double& time, const Vector<double>& x, double& flux)
 {
  //If the function pointer is zero return zero
  if(Flux_fct_pt == 0)
   {
    flux=0.0;
   }
  //Otherwise call the function
  else
   {
    (*Flux_fct_pt)(time,x,flux);
   }
 }
 
 
  private:
 
 
 /// \short Compute the element residual vector.
 /// flag=1(or 0): do (or don't) compute the Jacobian as well. 
 void fill_in_generic_residual_contribution_ust_heat_flux(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  unsigned flag);
 
 /// Function pointer to the (global) prescribed-flux function
 UnsteadyHeatPrescribedFluxFctPt Flux_fct_pt;
 
 ///The spatial dimension of the problem
 unsigned Dim;
 
 ///The index at which the unknown is stored at the nodes
 unsigned U_index_ust_heat;
 
 /// Id of first value assigned by this face element
 unsigned Melt_id;
 
 /// Pointer to non-default melt temperature
 double* Melt_temperature_pt;
 
};

//////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////// 

//===========================================================================
/// Constructor, takes the pointer to the "bulk" element and the 
/// index of the face to be created, and (optionally) the face element ID.
//===========================================================================
template<class ELEMENT>
UnsteadyHeatFluxPseudoMeltElement<ELEMENT>::
UnsteadyHeatFluxPseudoMeltElement(FiniteElement* const &bulk_el_pt, 
                        const int &face_index, const unsigned &id) : 
 FaceGeometry<ELEMENT>(), FaceElement()
{ 
 // Use default melt temperature
 Melt_temperature_pt=0;

#ifdef PARANOID
 {
  //Check that the element is not a refineable 3d element
  ELEMENT* elem_pt = new ELEMENT;
  //If it's three-d
  if(elem_pt->dim()==3)
   {
    //Is it refineable
    if(dynamic_cast<RefineableElement*>(elem_pt))
     {
      //Issue a warning
      OomphLibWarning(
       "This flux element will not work correctly if nodes are hanging\n",
       "UnsteadyHeatFluxPseudoMeltElement::Constructor",
       OOMPH_EXCEPTION_LOCATION);
     }
   }
 }
#endif
  
 // Let the bulk element build the FaceElement, i.e. setup the pointers 
 // to its nodes (by referring to the appropriate nodes in the bulk
 // element), etc.
 bulk_el_pt->build_face_element(face_index,this);
  
 // Melt flux stored at nodes
 Vector<unsigned> n_additional_values(nnode(),1); 
 
 // Now add storage and set the map containing
 // the position of the first entry of this face element's
 // additional values.
 Melt_id=id; 
 add_additional_values(n_additional_values,Melt_id);
 
 // Initialise the prescribed-flux function pointer to zero
 Flux_fct_pt = 0;
 
 // Extract the dimension of the problem from the dimension of 
 // the first node
 Dim = this->node_pt(0)->ndim();
 
 //Set up U_index_ust_heat. Initialise to zero, which probably won't change
 //in most cases, oh well, the price we pay for generality
 U_index_ust_heat = 0;
 
 //Cast to the appropriate UnsteadyHeatEquation so that we can
 //find the index at which the variable is stored
 //We assume that the dimension of the full problem is the same
 //as the dimension of the node, if this is not the case you will have
 //to write custom elements, sorry
 switch(Dim)
  {
   //One dimensional problem
  case 1:
  {
   UnsteadyHeatEquations<1>* eqn_pt = 
    dynamic_cast<UnsteadyHeatEquations<1>*>(bulk_el_pt);
   //If the cast has failed die
   if(eqn_pt==0)
    {
     std::string error_string =
      "Bulk element must inherit from UnsteadyHeatEquations.";
     error_string += 
      "Nodes are one dimensional, but cannot cast the bulk element to\n";
     error_string += "UnsteadyHeatEquations<1>\n.";
     error_string += 
      "If you desire this functionality, you must implement it yourself\n";
     
     throw OomphLibError(error_string,
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   //Otherwise read out the value
   else
    {
     //Read the index from the (cast) bulk element
     U_index_ust_heat = eqn_pt->u_index_ust_heat();
    }
  }
  break;
  
  //Two dimensional problem
  case 2:
  {
   UnsteadyHeatEquations<2>* eqn_pt = 
    dynamic_cast<UnsteadyHeatEquations<2>*>(bulk_el_pt);
   //If the cast has failed die
   if(eqn_pt==0)
    {
     std::string error_string =
      "Bulk element must inherit from UnsteadyHeatEquations.";
     error_string += 
      "Nodes are two dimensional, but cannot cast the bulk element to\n";
     error_string += "UnsteadyHeatEquations<2>\n.";
     error_string += 
      "If you desire this functionality, you must implement it yourself\n";
     
     throw OomphLibError(error_string,
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
     //Read the index from the (cast) bulk element.
     U_index_ust_heat = eqn_pt->u_index_ust_heat();
    }
  }
  break;
  
  //Three dimensional problem
  case 3:
  {
   UnsteadyHeatEquations<3>* eqn_pt = 
    dynamic_cast<UnsteadyHeatEquations<3>*>(bulk_el_pt);
   //If the cast has failed die
   if(eqn_pt==0)
    {
     std::string error_string =
      "Bulk element must inherit from UnsteadyHeatEquations.";
     error_string += 
      "Nodes are three dimensional, but cannot cast the bulk element to\n";
     error_string += "UnsteadyHeatEquations<3>\n.";
     error_string += 
      "If you desire this functionality, you must implement it yourself\n";
     
     throw OomphLibError(error_string,
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
     
    }
   else
    {
     //Read the index from the (cast) bulk element.
     U_index_ust_heat = eqn_pt->u_index_ust_heat();
    }
  }
  break;
  
  //Any other case is an error
  default:
   std::ostringstream error_stream; 
   error_stream <<  "Dimension of node is " << Dim 
                << ". It should be 1,2, or 3!" << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
   break;
  }
 
}
 
 
 
//===========================================================================
/// Compute the element's residual vector and the Jacobian matrix.
//===========================================================================
template<class ELEMENT>
void UnsteadyHeatFluxPseudoMeltElement<ELEMENT>::
fill_in_generic_residual_contribution_ust_heat_flux(
 Vector<double> &residuals, DenseMatrix<double> &jacobian, 
 unsigned flag)
{

 // Code isn't quite complete yet...
 if (flag==1) 
  {
   throw OomphLibError("Jacobian not yet implemented.\n",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 //Find out how many nodes there are
 const unsigned n_node = nnode();
   
 // Get continuous time from timestepper of first node
 double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
  
 //Set up memory for the shape and test functions
 Shape psif(n_node), testf(n_node);
   
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();
   
 //Set the Vector to hold local coordinates
 Vector<double> s(Dim-1);
   
 //Integer to store the local equation and unknowns
 int local_eqn=0;
   
 // Locally cache the index at which the variable is stored
 const unsigned u_index_ust_heat = U_index_ust_heat;
   
 //Loop over the integration points
 //--------------------------------
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
     
   //Assign values of s
   for(unsigned i=0;i<(Dim-1);i++)
    {
     s[i] = integral_pt()->knot(ipt,i);
    }

   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Find the shape and test functions and return the Jacobian
   //of the mapping
   double J = shape_and_test(s,psif,testf);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   //Need to find position to feed into flux function
   Vector<double> interpolated_x(Dim);
   
   // Melt flux
   double melt_flux=0.0;
   
   // Temperature
   double u=0.0;

   //Initialise to zero
   for(unsigned i=0;i<Dim;i++)
    {
     interpolated_x[i] = 0.0;
    }
   
   //Calculate quantities
   for(unsigned l=0;l<n_node;l++) 
    {
     // get the node pt
     Node* nod_pt = node_pt(l);

     // Cast to a boundary node
     BoundaryNodeBase *bnod_pt = 
      dynamic_cast<BoundaryNodeBase*>(nod_pt);
     
     // Get the index of the first nodal value associated with
     // this FaceElement
     unsigned first_index=
      bnod_pt->index_of_first_value_assigned_by_face_element(Melt_id);
     
     // Melt rate
     melt_flux+=nod_pt->value(first_index)*psif[l];
     
     // Temperature
     u+=nod_pt->value(u_index_ust_heat)*psif[l]; 
     
     //Loop over components
     for(unsigned i=0;i<Dim;i++)
      {
       interpolated_x[i] += nodal_position(l,i)*psif[l];
      }
    }
   
   //Get the imposed flux
   double flux=0.0;
   get_flux(time,interpolated_x,flux);
   
   //Now add to the appropriate equations
   
   //Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     local_eqn = nodal_local_eqn(l,u_index_ust_heat);
     /*IF it's not a boundary condition*/
     if(local_eqn >= 0)
      {
       //Add the prescribed flux terms
       residuals[local_eqn] -= (flux-melt_flux)*testf[l]*W;
      }
    }
  }

 // Collocation for melt rate:
 //---------------------------

 //Loop over the nodes
 for(unsigned l=0;l<n_node;l++)
  {
   // get the node pt
   Node* nod_pt = node_pt(l);
   
   // Cast to a boundary node
   BoundaryNodeBase *bnod_pt =
    dynamic_cast<BoundaryNodeBase*>(nod_pt);
   
   // Get the index of the first nodal value associated with
   // this FaceElement
   unsigned first_index=
    bnod_pt->index_of_first_value_assigned_by_face_element(Melt_id);
   local_eqn = nodal_local_eqn(l,first_index);
   
   /*IF it's not a boundary condition*/
   if(local_eqn >= 0)
    {
     double u=nod_pt->value(u_index_ust_heat);
     double melt_flux=nod_pt->value(first_index);
     
     // Piecewise linear variation
     if ((melt_temperature()-u)>melt_flux)
      {
       residuals[local_eqn]+=melt_flux;
      }
     else
      {
       residuals[local_eqn]+=(melt_temperature()-u);
      }
    }
  }

}
}

#endif
