#pragma once
#include <Global.hpp>

#include <myUtils.hpp>
#include <LinAlg.hpp>

#include "Geo2D_ConvexPolygon.hpp"

// geometrical object of type convex polygon
namespace MyFem::Mesh2D {
  
using namespace LinAlg;

class Mesh {
 public:
  Mesh();
 protected:
  Array<Vectord> nodeCoords_; // Array of (x, y) coords for global node Ids
  Array<Array<int>> elements_; // Array of elements, each consisting of a set of some global node IDs (e.g. three nodes for Tri3), with CCW node order
  Array<int> boundaryNodes_; // These are only the boundary nodes of the mesh; useful for applying BCs to e.g. all nodes with x_coord > 2.15, but without doing it on the inner nodes. This set is generally not equivalent to the given polygon (different IDs, nodes between polygon vertices)
};

class MeshTri3InConvexPolygon : Mesh {
  // Mesh generation algorithms
  void create_myAlgo_convexPolygonToTri3(Geo2D::ConvexPolygon& cPoly);
};

} // namespace MyFem::Mesh2D


//# TODO: We should also have a functionality where we can get only the boundary nodes. Because then we can from outside apply BCs to e.g. all nodes with x_coord > 2.15, but without doing it on the inner nodes. This set is generally not equivalent to the given polygon (different IDs, nodes between polygon vertices)
