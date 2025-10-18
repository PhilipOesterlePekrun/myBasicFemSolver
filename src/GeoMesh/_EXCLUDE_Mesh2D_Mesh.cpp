#include "Mesh2D_Mesh.hpp"

#include "Geo2D_Utils.hpp"

namespace MyFem::Mesh2D {

void MeshTri3InConvexPolygon::create_myAlgo_convexPolygonToTri3(Geo2D::ConvexPolygon& cPoly) {
  // boundary discretization
  // - we keep the boundary nodes as mesh nodes to make an exact mesh
  
  double minBoundarySpacing = 0.2;
  
  int nodeId = 0;
  auto& p = cPoly.pointsXY;
  Vectord actualBoundarySpacings;
  FOR(i, p.size()-1) {
    nodeCoords_.push_back(p(0));
    boundaryNodes_.push_back(nodeId++);
    
    const double distPoly = Geo2D::Utils::dist(p(i), p(i+1));
    const int numSegments = std::ceil(distPoly/minBoundarySpacing);
    const double actualBoundarySpacing = distPoly/numSegments;
    actualBoundarySpacings.push_back(actualBoundarySpacing);
    FOR(j, numSegments-1) {
      nodeCoords_.push_back(Geo2D::Utils::lInterp(p(i), p(i+1), actualBoundarySpacing*(j+1)));
      boundaryNodes_.push_back(nodeId++);
    }
    
    nodeCoords_.push_back(p(0));
    boundaryNodes_.push_back(nodeId++);
  }
  
  Array<int> virtualBoundaryNodes = boundaryNodes_; // the nodes on the "virtual boundary" within which meshing must still happen; starts same as actual boundary
  
  // condition for mesh complete/finished (no more meshing possible); assume true and get set false//#? or opposite?//#
  bool meshComplete = true;
  
  // the following will still also work for the initial actual (not really virtual) boundary, because at 180deg, we should be making a trisection. Basically, if angle<90deg we do no bisection, and if angle between 30+60n deg and 30+60(n+1) deg, we do n angle splits (bisection, trisection, whatever). The 30 is half so just make it an even chance to be smaller or greater than the perfect 60degs angles everywhere we optimally want
  while(1) {
    for(i=1; i<virtualBoundaryNodes.size(); ++i) {
      int nodeId = virtualBoundaryNodes(i);
      Array<Vectord> threePointsCoords;

      for(j=-1; j<3; ++j) {
        threePointsCoords.push_back(nodeCoords_(nodeId(i+j)));
      }
      double cornerAngle = Geo2D::Utils::lAngle(threePointsCoords);

      // very inefficient, but for now maybe the strat for inexact points on top of eachother is to do a full search and compare dist of all nodes. If there is a chance of the points ending up not super close but fairly close, we might need another parameter for this idk. Can be made more efficient later I'm sure
    }
  }
  
  // the first contiguous part of the nodeIDs are all boundary nodes, so we simply copy before putting the inner nodes
}

} // namespace MyFem::Mesh2D