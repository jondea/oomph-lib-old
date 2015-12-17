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
//Header file for a class that is used to extract the faces of bulk
//elements and represent them as geometric objects, primarily for use
//in FSI problems

//Include guards to prevent multiple inclusion of the header
#ifndef OOMPH_FACE_ELEMENT_AS_GEOMETRIC_OBJECT_HEADER
#define OOMPH_FACE_ELEMENT_AS_GEOMETRIC_OBJECT_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include<algorithm>

//Include the geometric object header file
#include "geom_objects.h"
#include "shape.h"
#include "multi_domain.h"

namespace oomph
{

//=======================================================================
/// Class that is used to create FaceElement from bulk elements and to
/// provide these FaceElement with a geometric object representation.
/// The local coordinates of the FaceElements are used as the intrinisic
/// coordinates for its GeomObject representation.
/// 
/// These elements are used primarily to set up the interaction terms
/// in FSI problems and are expected to be created from meshes so 
/// that they lie on a particular boundary of the mesh.
//=======================================================================
 template<class ELEMENT>
 class FaceElementAsGeomObject : public virtual FaceGeometry<ELEMENT>,
                                 public virtual FaceElement,
                                 public virtual ElementWithExternalElement
 {
 public:

  ///\short Constructor which takes a pointer to a "bulk" element,
  /// to which this element is attached. The face index, indicates
  /// the face (of the bulk element) that is to be constructed.
  /// Note that this element tends to be constructed  
  /// by the doubly templated Mesh::build_face_mesh() and therefore
  /// has to have the same interface as the the generic FaceElement
  /// constructor. Hence the boundary number (within the mesh) on which this 
  /// element is located must be setup afterwards!
  FaceElementAsGeomObject(FiniteElement* const &element_pt, 
                          const int &face_index) :
   FaceGeometry<ELEMENT>(), FaceElement(), 
   //The geometric object has an intrinsic dimension one less than 
   //the "bulk" element, but the actual dimension of the problem remains 
   //the same
   //GeomObject(element_pt->dim()-1,element_pt->nodal_dimension()),
   ElementWithExternalElement()
   { 
    //Attach the geometrical information to the element. N.B. This function
    //also assigns nbulk_value from the required_nvalue of the bulk element
    element_pt->build_face_element(face_index,this);
    GeomObject::set_nlagrangian_and_ndim(element_pt->dim()-1,
                                         element_pt->nodal_dimension());
   }


  /// Broken copy constructor
  FaceElementAsGeomObject(const FaceElementAsGeomObject&) 
   { 
    BrokenCopy::broken_copy("FaceElementAsGeomObject");
   } 
 
  /// Broken assignment operator
  void operator=(const FaceElementAsGeomObject&) 
   {
    BrokenCopy::broken_assign("FaceElementAsGeomObject");
   }
 
  /// \short The "global" intrinsic coordinate of the element when
  /// viewed as part of a geometric object should be given by
  /// the FaceElement representation, by default
  double zeta_nodal(const unsigned &n, const unsigned &k,           
                    const unsigned &i) const 
   {return FaceElement::zeta_nodal(n,k,i);}     


  /// \short How many items of Data does the shape of the object depend on?
  /// None! We're dealing with a pure geometric FiniteElement!
  unsigned ngeom_data() const {return 0;}
 
  /// \short Return pointer to the j-th Data item that the object's 
  /// shape depends on. Object doesn't depend on any geom Data
  /// so we die if this gets called.
  Data* geom_data_pt(const unsigned& j)
   {
    std::ostringstream error_message;
    error_message 
     << "FaceElementAsGeomObject::geom_data_pt() is deliberately broken\n"
     << "as it does not depend on any geometric Data" << std::endl;
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
    // Dummy return
    return 0;
   }

  /// Override fill in contribution to jacobian, nothing should be done
  void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                        DenseMatrix<double> &jacobian) 
   {
    std::ostringstream warn_message;
    warn_message
     << "Warning: You have just called the empty function \n"
     << 
     "fill_in_contribution_to_jacobian() for a FaceElementAsGeometricObject.\n"
     << 
     "These Elements should only be used to setup interactions, so should\n"
     << 
     "not be included in any jacobian calculations\n";
   
    OomphLibWarning(
     warn_message.str(),
     "FaceElementAsGeometricObject::fill_in_contribution_to_jacobian()",
     OOMPH_EXCEPTION_LOCATION);
   }

  /// \short Function to describe the local dofs of the element. The ostream 
  /// specifies the output stream to which the description 
  /// is written; the string stores the currently 
  /// assembled output that is ultimately written to the
  /// output stream by Data::describe_dofs(...); it is typically
  /// built up incrementally as we descend through the
  /// call hierarchy of this function when called from 
  /// Problem::describe_dofs(...)
  void describe_local_dofs(std::ostream& out,
                           const std::string& current_string) const
   {
    // Call the ElementWithExternalElement's describe function
    ElementWithExternalElement::describe_local_dofs(out,current_string);
   }

  /// Unique final overrider needed for assign_all_generic_local_eqn_numbers
  void assign_all_generic_local_eqn_numbers(const bool &store_local_dof_pt)
   {
    // Call the ElementWithExternalElement's assign function
    ElementWithExternalElement::assign_all_generic_local_eqn_numbers(
     store_local_dof_pt);
   }
 

 };


 //============================================================================
 /// \short A class to do comparison of the elements by lexicographic
 /// ordering, based on the boundary coordinates at the element's first node. 
 //============================================================================
 template<class ELEMENT>
 class CompareBoundaryCoordinate
 {
 public:
   
  ///The actual comparison operator
  int operator() (GeneralisedElement* const &element1_pt,
                  GeneralisedElement* const &element2_pt)
   {
    //OK Dynamic cast the elements
    FaceElementAsGeomObject<ELEMENT> *cast_element1_pt = 
     dynamic_cast<FaceElementAsGeomObject<ELEMENT>*>(element1_pt);
    FaceElementAsGeomObject<ELEMENT> *cast_element2_pt = 
     dynamic_cast<FaceElementAsGeomObject<ELEMENT>*>(element2_pt);

#ifdef PARANOID
    if (cast_element1_pt==0)
     {
      std::ostringstream error_message;
      error_message 
       << "Failed to cast element1_pt to a FaceElementAsGeomObject"
       << std::endl;
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }

    if (cast_element2_pt==0)
     {
      std::ostringstream error_message;
      error_message 
       << "Failed to cast element2_pt to a FaceElementAsGeomObject"
       << std::endl;
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif


    // Warning that this still needs to be generalised to higher
    // dimensions (don't want to implement it until I can test it
    // -- at the moment, the ordering isn't particularly important
    // anyway...
//      if (cast_element1_pt->dim()!=1)
//       {
//        std::ostringstream warn_message;
//        warn_message
//         << "Warning: Ordering of elements is currently based on their \n"
//         << "zero-th surface coordinate. This may not be appropriate for\n"
//         << cast_element1_pt->dim() << "-dimensional elements. \n";
//         OomphLibWarning(warn_message.str(),
//                         "CompareBoundaryCoordinate::()",
//                         OOMPH_EXCEPTION_LOCATION);
//       }


    return 
     cast_element1_pt->zeta_nodal(0,0,0) < 
     cast_element2_pt->zeta_nodal(0,0,0);
   }
 };
 
 


}

#endif
