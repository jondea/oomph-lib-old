#ifndef OOMPH_ODE_ELEMENTS_H
#define OOMPH_ODE_ELEMENTS_H

#include "../generic/oomph_definitions.h"
#include "../generic/oomph_utilities.h"

#include "../generic/matrices.h"
#include "../generic/Vector.h"
#include "../generic/elements.h"
#include "../generic/timesteppers.h"

namespace oomph
{

 /// Element for integrating an initial value ODE
 class ODEElement : public GeneralisedElement
 {

 public:


  /// Constructor: Pass timestepper and a solution function pointer
  ODEElement(TimeStepper* timestepper_pt,
             SolutionFunctorBase* exact_solution_pt)
  {
   Exact_solution_pt = exact_solution_pt;

   Vector<double> exact = this->exact_solution(0);
   unsigned nvalue = exact.size();

   add_internal_data(new Data(timestepper_pt, nvalue));

   Use_fd_jacobian = false;
  }


  virtual ~ODEElement() {}

  unsigned nvalue() const
  {
   return internal_data_pt(0)->nvalue();
  }

  /// Get residuals
  void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
   // Get pointer to one-and-only internal data object
   Data* dat_pt = internal_data_pt(0);

   // Get it's values
   Vector<double> u(nvalue(), 0.0);
   dat_pt->value(u);

   // Get timestepper
   TimeStepper* timestepper_pt = dat_pt->time_stepper_pt();

   // Get continuous time
   double t = timestepper_pt->time();

   Vector<double> deriv = derivative_function(t, u);
   for(unsigned j=0, nj=deriv.size(); j<nj; j++)
    {
     // Get dudt approximation from timestepper
     double dudt = timestepper_pt->time_derivative(1, dat_pt, j);

     // Residual is difference between the exact derivative and our
     // timestepper's derivative estimate.
     residuals[j] = deriv[j] - dudt;
    }
  }

  void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                        DenseMatrix<double>& jacobian)
  {
   // Get residuals
   fill_in_contribution_to_residuals(residuals);

   if(Exact_solution_pt->have_jacobian() && !Use_fd_jacobian)
    {
     // get df/du jacobian
     double t = internal_data_pt(0)->time_stepper_pt()->time();
     Vector<double> dummy, u(nvalue(), 0.0);
     internal_data_pt(0)->value(u);
     Exact_solution_pt->jacobian(t, dummy, u, jacobian);

     // We need jacobian of residual = f - dudt so subtract diagonal
     // (dudt)/du term.
     const double a = internal_data_pt(0)->time_stepper_pt()->weight(1,0);
     const unsigned n = nvalue();
     for(unsigned i=0; i<n; i++)
      {
       jacobian(i, i) -= a;
      }
    }
   else
    {
     // Use FD for jacobian
     GeneralisedElement::fill_in_jacobian_from_internal_by_fd
      (residuals, jacobian, true);
    }
  }

  void fill_in_contribution_to_mass_matrix(Vector<double>& residuals,
                                           DenseMatrix<double>& mm)
  {
   fill_in_contribution_to_residuals(residuals);
   for(unsigned j=0, nj=nvalue(); j<nj; j++)
    {
     mm(j, j) = 1;
    }
  }

  /// Exact solution
  Vector<double> exact_solution(const double& t) const
  {
#ifdef PARANOID
   if(Exact_solution_pt == 0)
    {
     throw OomphLibError("No exact solution function",
                         OOMPH_EXCEPTION_LOCATION,
                         OOMPH_CURRENT_FUNCTION);
    }
#endif
   Vector<double> dummy_x;
   return (*Exact_solution_pt)(t, dummy_x);
  }

  /// Exact solution
  Vector<double> derivative_function(const double& t,
                                     const Vector<double>& u)
  {
#ifdef PARANOID
   if(Exact_solution_pt == 0)
    {
     throw OomphLibError("No derivative function",
                         OOMPH_EXCEPTION_LOCATION,
                         OOMPH_CURRENT_FUNCTION);
    }
#endif
   Vector<double> dummy_x;
   return Exact_solution_pt->derivative(t, dummy_x, u);
  }

  SolutionFunctorBase* Exact_solution_pt;

  bool Use_fd_jacobian;
 };

}

#endif