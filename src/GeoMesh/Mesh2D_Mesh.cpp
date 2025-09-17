#include "Mesh2D_Mesh.hpp"

#include "Geo2D_Utils.hpp"

namespace MyFem::Mesh2D {

void MeshTri3InConvexPolygon::create_myAlgo_convexPolygonToTri3(Geo2D::ConvexPolygon& cPoly) {
  // boundary discretization
  // - we keep the boundary nodes as mesh nodes to make an exact mesh
  
  double minBoundarySpacing = 0.2;
  
  int nodeId = 0;
  auto& p = cPoly.pointsXY;
  FOR(i, p.size()-1) {
    nodeCoords_.push_back(p(0));
    boundaryNodes_.push_back(nodeId++);
    
    const double distPoly = Geo2D::Utils::dist(p(i), p(i+1));
    const int numSegments = std::ceil(distPoly/minBoundarySpacing);
    const double actualBoundarySpacing = distPoly/numSegments;
    FOR(j, numSegments-1) {
      nodeCoords_.push_back(Geo2D::Utils::lInterp(p(i), p(i+1), actualBoundarySpacing*(j+1)));
      boundaryNodes_.push_back(nodeId++);
    }
    
    nodeCoords_.push_back(p(0));
    boundaryNodes_.push_back(nodeId++);
  }
  // the first contiguous part of the nodeIDs are all boundary nodes, so we simply copy before putting the inner nodes
}

} // namespace MyFem::Mesh2D