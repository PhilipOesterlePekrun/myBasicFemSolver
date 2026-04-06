// // MAYBE THIS CLASS MAKES NO SENSE TO HAVE

#pragma once
#include "Global.hpp"

#include "Element_Tri3.hpp"

namespace MyFem {

namespace Element {

class Tri3Dyn : public Tri3 {
// // Values in time
  vector<Vectord> X_n;
  
  Vectord V_0;
  ///vector<Vectord> V_n;
  
// // Misc
  bool pseudoDynamic;
  
// // // ctor(s)
 public:
  Tri3Dyn(const vector<int>& globalNodeIds_, const Vectord& X_0_, const Vectord& V_0_, bool planeStressElsePlaneStrain_, std::function<Vectord(double, double)> youngPoisson_x_, std::function<double(double, double)> density_x_, bool pseudoDynamic_ = false)
  : Tri3(globalNodeIds_, X_0_, planeStressElsePlaneStrain_, youngPoisson_x_, density_x_), V_0(V_0_), pseudoDynamic(pseudoDynamic_) {}
  
  Tri3Dyn(const Vectord& globalX_0, const vector<int>& globalNodeIds_, const Vectord& V_0_, bool planeStressElsePlaneStrain_, std::function<Vectord(double, double)> youngPoisson_x_, std::function<double(double, double)> density_x_, bool pseudoDynamic_ = false)
  : Tri3(globalX_0, globalNodeIds_, planeStressElsePlaneStrain_, youngPoisson_x_, density_x_), V_0(V_0_), pseudoDynamic(pseudoDynamic_) {}
  
// // Dynamics
 private:
  Matrix2d integrandMmat_xi(double xi0, double xi1) {
    Matrix2d mat(ndof_, ndof_);
    double detJ = detMat2d(jacobian_xi(xi0, xi1));
    if(detJ<0) MyUtils::Db::warn("detJ is negative! It is "+std::to_string(detJ));
    Vectord N({1.0 - xi0 - xi1, xi0, xi1});

    // physical position of this quadrature point
    double x0 = 0.0;
    double x1 = 0.0;
    FOR(a, globalNodeIds_.size()) {
      x0 += N(a) * X_0_(ndofn_ * a + 0);
      x1 += N(a) * X_0_(ndofn_ * a + 1);
    }

    double rho = density_x(x0, x1);

    FOR(mLeft, ndof_) {
      int a = mLeft / ndofn_;   // node index
      int i = mLeft % ndofn_;   // component index

      FOR(mRight, ndof_) {
        int b = mRight / ndofn_;
        int j = mRight % ndofn_;

        mat(mLeft, mRight) = (i == j ? rho * N(a) * N(b) * detJ : 0.0);
      }
    }

    return mat;
  }
  
 public:
  // Get the mass matrix
  Matrix2d Mmat() {
    std::string outString = "Element ";
    FOR(i, globalNodeIds_.size()-1)
      outString += std::to_string(globalNodeIds_[i]) + "-";
    outString += std::to_string(globalNodeIds_[globalNodeIds_.size()-1]);
    outString += ": getting Mmat()\n";
    MyUtils::Db::pr(outString);
    
    Matrix2d result(ndof_, ndof_);
    
    for(int i=0;i<ndof_; ++i)
      for(int j=0;j<ndof_; ++j) {
        auto g = [&, i, j](double xi0)
        {
          auto f = [&, i, j, xi0](double xi1)
          {
            return integrandMmat_xi(xi0, xi1)(i, j);
          };
          
          // inner integration
          return MyUtils::NumIntegration::gaussianQuadrature(f, 1);
        };
        
        // outer integration
        result(i, j) = MyUtils::NumIntegration::gaussianQuadrature(g, 1);

      }
      
    return result;
  }
}; // class Tri3

} // namespace Element

} // namespace MyFem
