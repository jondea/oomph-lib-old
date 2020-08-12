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
#ifndef OOMPH_WRAP_PML_HEADER
#define OOMPH_WRAP_PML_HEADER

#include "simple_rectangular_quadmesh.h"

namespace oomph
{
namespace WrapPMLHelper
{
  /// Sorter which compares boundary coordinate of two nodes along boundary b
  /// Called as BoundaryNodeSorter(boundary_id)(node1_pt, node2_pt)
  struct BoundaryNodeSorter {
    BoundaryNodeSorter(int b) { this->b = b; }

    template <class NODE>
    bool operator () (NODE* node1_pt, NODE* node2_pt)
    {
      Vector<double> zeta1(1);
      dynamic_cast<BoundaryNode<NODE>*>(node1_pt)->
        get_coordinates_on_boundary(b,zeta1);

      Vector<double> zeta2(1);
      dynamic_cast<BoundaryNode<NODE>*>(node2_pt)->
        get_coordinates_on_boundary(b,zeta2);

      return zeta1 < zeta2;
    }

    int b;
  };

  /// Calculate an interior value of the boundary coord of a face element by
  /// summing the value of all its nodes. Could also be called a centre
  template <class ELEMENT>
  double face_element_interior_boundary_coord(
    const FaceElementAsGeomObject<ELEMENT>* face_el,
    const unsigned& b)
  {
    Vector<double> zeta_vec(1);
    double zeta = 0.0;
    unsigned nnode = face_el->nnode();
    for(unsigned i=0; i<nnode; i++)
    {
      face_el->node_pt(i)->get_coordinates_on_boundary(b,zeta_vec);
      zeta += zeta_vec[0];
    }
    return zeta/double(nnode);
  }

  /// Sorter which compares an interior point of two face elements on a boundary
  /// specified by b.
  /// Called as FaceElementSorter(boundary_id)(face_el1_pt, face_el2_pt)
  struct FaceElementSorter {
    FaceElementSorter(int b) { this->b = b; }

    template <class ELEMENT>
    bool operator () (
      FaceElementAsGeomObject<ELEMENT>* face_el1_pt,
      FaceElementAsGeomObject<ELEMENT>* face_el2_pt)
    {
      // Compare an interior boundary coordinate point of the two face elementss
      return (   face_element_interior_boundary_coord(face_el1_pt, b)
               < face_element_interior_boundary_coord(face_el2_pt, b) );
    }

    int b;
  };

  template <class ELEMENT>
  bool is_face_el_oriented_same_as_boundary_coord(
    FaceElementAsGeomObject<ELEMENT>* face_el_pt, const unsigned& b)
  {
    // Get zeta at start of face element
    Vector<double> zeta_start(1);
    face_el_pt->node_pt(0)->get_coordinates_on_boundary(b,zeta_start);

    // Get zeta at end of face element
    unsigned nnode = face_el_pt->nnode();
    Vector<double> zeta_end(1);
    face_el_pt->node_pt(nnode-1)->get_coordinates_on_boundary(b,zeta_end);

    // Orientation is same as boundary coord if start is less than end
    return zeta_start[0] < zeta_end[0];
  }

  /// Make mesh truly periodic by attaching two x ends of mesh together
  /// Note, assumes Quad elements with
  void stitch_x_boundaries(Mesh* mesh_pt, const unsigned n_el_x,
    const unsigned n_pml, const unsigned nnode_1d)
  {
    // Loop over the rows in y direction
    for(unsigned i=0; i<n_pml; i++)
    {
      // Get pointers to the leftmost and rightmost elements on this row
      FiniteElement* leftmost_el_pt =
        dynamic_cast<FiniteElement*>(mesh_pt->element_pt(i*n_el_x));
      FiniteElement* rightmost_el_pt =
        dynamic_cast<FiniteElement*>(mesh_pt->element_pt((i+1)*n_el_x-1));

      // Loop over nodes along the x-boundaries
      for(unsigned j=0; j<nnode_1d; j++)
      {
        // Mark the rightmost node as obsolete for deletion
        rightmost_el_pt->node_pt(nnode_1d-1+j*nnode_1d)->set_obsolete();

        // Change pointer in rightmost elements to the leftmost node
        rightmost_el_pt->node_pt(nnode_1d-1+j*nnode_1d) =
          leftmost_el_pt->node_pt(j*nnode_1d);
      }
    }

    // Delete the actual nodes and references to them in the mesh and elements
    mesh_pt->prune_dead_nodes();
  }  

  template<class BULK_ELEMENT, class PML_ELEMENT>
  Mesh* create_wrapped_pml_mesh(Mesh* bulk_mesh_pt,
                                const Vector<unsigned>& boundary_ids,
                                const unsigned& n_pml,
                                const bool hang_interior_nodes = false,
                                TimeStepper* time_stepper_pt=
                                &Mesh::Default_TimeStepper)
  {
    oomph_info << "Warning: WrapPMLHelper::create_wrapped_pml_mesh assumes"
               << "that boundaries are " << std::endl << "connected and have a "
               << "strictly increasing boundary coordinate" << std::endl;

    // Number of boundaries passed in
    unsigned n_boundary = boundary_ids.size();

    // Storage for constructed face elements, their orientations and number
    Vector<FaceElementAsGeomObject<BULK_ELEMENT>*> face_element_pt(0);
    std::vector<bool> face_element_orientation(0);
    unsigned nface_element = 0;

    // Storage for pointers to nodes on boundaries and number
    Vector<Node*> boundary_node_pt(0);
    unsigned nboundary_node = 0;

    // Loop over all boundaries
    for(unsigned ib=0; ib<n_boundary; ib++)
    {
      // Cache boundary id
      unsigned b = boundary_ids[ib];


      //      Face elements
      //--------------------------

      // Storage for and number of face elements on boundary b
      unsigned nface_element_on_b = bulk_mesh_pt->nboundary_element(b);
      Vector<FaceElementAsGeomObject<BULK_ELEMENT>*>
        face_element_on_b_pt(nface_element_on_b);

      // Construct all face elements on boundary b and store
      for(unsigned e=0; e<nface_element_on_b; e++)
      {
        // Get pointer to the bulk element that is adjacent to boundary b
        FiniteElement* bulk_elem_pt =
          bulk_mesh_pt->boundary_element_pt(b,e);

        // Find the index of the face of element e along boundary b
        int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

        // Build the corresponding face element and store pointer in vector
        face_element_on_b_pt[e] = new FaceElementAsGeomObject<BULK_ELEMENT>(
                                        bulk_elem_pt, face_index);
      }

      // Sort the face elements based on an interior boundary coord
      std::sort(face_element_on_b_pt.begin(), face_element_on_b_pt.end(),
                FaceElementSorter(b));

      // Append the sorted face elements to the big vector and work out their
      // orientation
      unsigned nface_element_before_b = nface_element;
      nface_element += nface_element_on_b;
      face_element_pt.resize(nface_element);
      face_element_orientation.resize(nface_element);
      for(unsigned e=0; e<nface_element_on_b; e++)
      {
        // Store face element pointers from this boundary in big vector
        face_element_pt[nface_element_before_b+e] = face_element_on_b_pt[e];

        // Determine orientation
        face_element_orientation[nface_element_before_b+e] =
          is_face_el_oriented_same_as_boundary_coord(
            face_element_pt[nface_element_before_b+e], b);
      }

      //      Face elements
      //--------------------------

      // Store all nodes on boundary b
      unsigned nboundary_node_on_b = bulk_mesh_pt->nboundary_node(b);
      Vector<Node*> boundary_node_on_b_pt(nboundary_node_on_b);
      for(unsigned n=0; n<nboundary_node_on_b; n++)
      {
        boundary_node_on_b_pt[n] = dynamic_cast<Node*>(bulk_mesh_pt->boundary_node_pt(b,n));
      }

      // Sort nodes by their boundary coordinates
      std::sort(boundary_node_on_b_pt.begin(), boundary_node_on_b_pt.end(),
                BoundaryNodeSorter(b));

      // Append the sorted nodes to the big vector
      // (except the last node on the boundary, which is present on the next b)
      unsigned nboundary_node_before_b = nboundary_node;
      nboundary_node += nboundary_node_on_b-1;
      boundary_node_pt.resize(nboundary_node);
      for(unsigned e=0; e<nboundary_node_on_b-1; e++)
      {
        boundary_node_pt[nboundary_node_before_b+e] = boundary_node_on_b_pt[e];
      }
    }

    // Cache various properties from representative face element
    double face_s_min = face_element_pt[0]->s_min();
    double face_s_max = face_element_pt[0]->s_max();

    // Now we know how many elements we need, create the PML mesh
    Mesh* pml_mesh_pt = new RectangularQuadMesh<PML_ELEMENT >(
                              nface_element, n_pml, 1.0, 1.0);

    // Get a representative element from PML mesh and cache various properties
    PML_ELEMENT* repr_pml_element_pt =
      dynamic_cast<PML_ELEMENT*>(pml_mesh_pt->element_pt(0));
    double pml_s_min = repr_pml_element_pt->s_min();
    double pml_s_max = repr_pml_element_pt->s_max();
    unsigned nnode_1d = repr_pml_element_pt->nnode_1d();

    // Make truly periodic
    stitch_x_boundaries(pml_mesh_pt, nface_element, n_pml, nnode_1d);

    // Enable PML elements
    // -------------------
    unsigned through_the_pml_index = 1;
    unsigned n_element = pml_mesh_pt->nelement();
    for(unsigned e=0; e<n_element; e++)
    {
      // Upcast from GeneralisedElement to Pml Helmholtz bulk element
      Conformal2DPMLElement *el_pt =
      dynamic_cast<Conformal2DPMLElement*>(pml_mesh_pt->element_pt(e));
#ifdef PARANOID
      if (el_pt==0)
      {
        throw OomphLibError("PML element must inherit from"
                            "Conformal2DPMLElement to use wrapped PML helper",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Get local coordinate of this knot point
      Vector<double> s(2,0.0);

      // Get undeformed position through the PML at innermost point of element
      s[through_the_pml_index] = pml_s_min;
      double nun_inner = el_pt->interpolated_x(s, through_the_pml_index);

      // Get undeformed position through the PML at outermost point of element
      s[through_the_pml_index] = pml_s_max;
      double nun_outer = el_pt->interpolated_x(s, through_the_pml_index);

      // Enable PML in correct direction and note the inner and outer values
      // of through the PML coordinate at the inner and outer part of element
      el_pt->enable_pml(through_the_pml_index, nun_inner, nun_outer);
    }

    //   Set PML node positions
    // ---------------------------

    // Start off as if we were just in the last face element
    unsigned e = nface_element-1;
    unsigned previous_e = e;
    // ...we go the opposite direction w.r.t. the pml elements
    unsigned inner_pml_e = 0;
    unsigned previous_inner_pml_e = inner_pml_e;

    // Loop over boundary nodes
    for(unsigned i=0; i<nboundary_node; i++)
    {

      // Store previous element indices
      previous_e = e;
      previous_inner_pml_e = inner_pml_e;

      // If node is at end of element, increment face element counter by 1
      if(i%(nnode_1d-1) == 0)
      {
        e = (e+1)%nface_element;
        //...we go the opposite direction w.r.t. the pml elements
        inner_pml_e = (inner_pml_e-1+nface_element)%nface_element;
      }

      Vector<double> unit_normal(2);

      // Find outer unit normal, if on a boundary between two elements, take
      // the average of the two limits
      if(i%(nnode_1d-1) == 0) // On an element boundary
      {
        Vector<double> s(1, face_s_max);

        // At the max or min, depending on the element orientation
        if(face_element_orientation[e]) s[0] = face_s_min;
        else                            s[0] = face_s_max;

        // Evaluate unit normal
        Vector<double> unit_normal1(2);
        face_element_pt[e]->outer_unit_normal(s, unit_normal1);

        // At the max or min, depending on the element orientation
        if(face_element_orientation[previous_e]) s[0] = face_s_max;
        else                                     s[0] = face_s_min;

        // Evaluate unit normal
        Vector<double> unit_normal2(2);
        face_element_pt[previous_e]->outer_unit_normal(s, unit_normal2);

        // Find average unit normal
        unit_normal[0] = (unit_normal1[0] + unit_normal2[0])/2.0;
        unit_normal[1] = (unit_normal1[1] + unit_normal2[1])/2.0;
      }
      else // Interior of element, just find unit norml
      {
        // Find node number, changes if orientation is reversed
        unsigned face_el_node_index = i%(nnode_1d-1);
        if(!face_element_orientation[e])
        {
          face_el_node_index = nnode_1d - 1 - i%(nnode_1d-1);
        }

        // Find local coord of node and find unit normal at this point
        Vector<double> s(1);
        face_element_pt[e]->local_coordinate_of_node(face_el_node_index,s);
        face_element_pt[e]->outer_unit_normal(s, unit_normal);
      }

      unsigned pml_e = inner_pml_e;
      unsigned previous_pml_e = previous_inner_pml_e;
      PML_ELEMENT* pml_el_pt =
        dynamic_cast<PML_ELEMENT*>(pml_mesh_pt->element_pt(pml_e));

      // Find position of the boundary node in bulk
      Vector<double> x_inner(2);
      boundary_node_pt[i]->position(x_inner);

      // Loop over elements through the PML (y-direction)
      for(unsigned j=0; j<n_pml; j++)
      {
        // Shift up one row
        pml_e = inner_pml_e + j*nface_element;
        previous_pml_e = previous_inner_pml_e + j*nface_element;
        pml_el_pt = dynamic_cast<PML_ELEMENT*>(
          pml_mesh_pt->element_pt(pml_e));

        // Loop over nodes through the PML in this element
        for(unsigned k=0; k<nnode_1d; k++)
        {
          // Index of node in the pml_element, remember PML element goes in
          // opposite direction to boundary nodes (indexed by i)
          unsigned n = nnode_1d-1 - i%(nnode_1d-1) + k*nnode_1d;
          Node* pml_node_pt = pml_el_pt->node_pt(n);

          // At inner PML boundary
          if(j==0 && k==0)
          {
            // Mark point at inner PML boundary as obsolete
            pml_node_pt->set_obsolete();

            // Set pointer in PML element to the node in the bulk
            pml_el_pt->node_pt(n) = boundary_node_pt[i];

            // If node is in two elements, change node pointer there too
            if(n == (nnode_1d-1))
            {
              PML_ELEMENT* pml_el_pt2 =
                dynamic_cast<PML_ELEMENT*>(
                  pml_mesh_pt->element_pt(previous_pml_e));

              pml_el_pt2->node_pt(0) = boundary_node_pt[i];
            }
          }
          else if(k!=0) // In rest of PML, move nodes (k=0 => already done)
          {
            // Find the through the PML coordinate nu
            double nu = pml_node_pt->x(through_the_pml_index);

            // Set PML node position to x_inner + nu*unit_normal
            pml_node_pt->x(0) = x_inner[0] + nu*unit_normal[0];
            pml_node_pt->x(1) = x_inner[1] + nu*unit_normal[1];

            // If we want to hang_interior_nodes and we aren't at the last node
            if(hang_interior_nodes && !(j==n_pml-1 && k==nnode_1d-1))
            {
              // Make a new HangInfo object
              HangInfo* hanginfo_pt = new HangInfo(2);

              // Inner master is the innermost corresponding node in the inner
              // pml element
              PML_ELEMENT* inner_el_pt =
                dynamic_cast<PML_ELEMENT*>(
                  pml_mesh_pt->element_pt(inner_pml_e));

              Node* inner_master_node_pt = inner_el_pt->
                          node_pt(nnode_1d-1 - i%(nnode_1d-1));

              // Add outer master node to new HangInfo
              hanginfo_pt->set_master_node_pt(0, inner_master_node_pt, 1-nu);

              // Outer master is the outermost corresponding node in the outer
              // pml element
              PML_ELEMENT* outer_el_pt =
                dynamic_cast<PML_ELEMENT*>(
                  pml_mesh_pt->element_pt((n_pml-1)*nface_element+inner_pml_e));

              Node* outer_master_node_pt = outer_el_pt->
                node_pt(nnode_1d-1 - i%(nnode_1d-1) + (nnode_1d-1)*nnode_1d);

              // Add outer master node to new HangInfo
              hanginfo_pt->set_master_node_pt(1, outer_master_node_pt, nu);

              // Geometrically pin the hanging node
              pml_node_pt->set_hanging_pt(hanginfo_pt,-1);

            }
          }

        } // End loop over nodes through the PML in this element
      } // End loop over elements through the PML (y-direction)

    } // End loop over boundary nodes

    // Delete storage for nodes which have been replaced by boundary nodes
    pml_mesh_pt->prune_dead_nodes();

    // Delete face elements we constructed
    for(unsigned i=0; i<nface_element; i++)
    {
      delete face_element_pt[i];
      face_element_pt[i]=0;
    }

    return pml_mesh_pt;
  }
}


// // Map out a representative element, work out where nodes are
// ASSOCIATED_PML_QUAD_ELEMENT* rep_element = pml_mesh_pt->element_pt(0);
// unsigned nnode_1d = rep_element->nnode_1d();
// Vector<unsigned> outer_node_index(nnode_1d);
// for(unsigned i=0; i<nnode_1d; i++)
// {
//   Vector<double>s_outer(2);
//   rep_element->local_coordinate_of_node(i, s_outer);
//   s_outer[1] = rep_element->s_max();
//   for(unsigned j=nnode_1d*nnode_1d-1; j>=nnode_1d; j--)
//   {
//     Vector<double>s(2);
//     rep_element->local_coordinate_of_node(j, s);
//     s_outer[1] = rep_element->s_max();
//     if(abs(s_outer[0]-s[0])<1.0e-12 && abs(s_outer[1]-s[1])<1.0e-12)
//     {
//       outer_node_index[i] = j;
//       break;
//     }
//   }
// }
}
#endif
