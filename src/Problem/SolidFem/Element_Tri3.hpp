#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include <LinAlg.hpp>
#include <Timer.hpp>

namespace MyFem {

namespace Element {
  
using namespace LinAlg; // mTODO: consider removing this and replacing it by specific usings. Maybe with a good macro I can reuse in other files too; In any case, except for structures like matrix vector array, I will try to be explicit about LinAlg::

class Tri3 {
// Characteristic of the element type
 public:
  static constexpr int nnode_ = 3;
  static constexpr int ndofn_ = 2;
  static constexpr int ndof_ = ndofn_*nnode_;
 private:
  static constexpr int dim_ = ndofn_;
  static constexpr double Xi_[nnode_][ndofn_] = {{0, 0}, {1, 0}, {0, 1}};
  /*
      |\
    ^ | \
xi1 | |  \
    | |___\
      --->
       xi0 (and first dof)
       
   N2    
    |\
    | \
    |  \
    |___\
   N0   N1
  
  */

 public:
  // This is the local (pos) to global (val) mapping
  Array<int> globalNodeIds_;
  Vectord X_0_;
  
  Array<int> getGlobalDofIds() {
    Array<int> arr = Array<int>();
    FOR(i, nnode_) {
      arr.push_back(2*globalNodeIds_(i));
      arr.push_back(2*globalNodeIds_(i) + 1);
    }
    return arr;
  }

// Node/dof mappings
 public:
  // node and node-local dof to element-global dof
  int nodeldof2dof(int n, int i) {
    return n*ndofn_ + i;
  }
  int nodeldof2dof(Array<int> ni) {
    return ni(0)*ndofn_ + ni(1);
  }
  // element-global dof to node and node-local dof
  Array<int> dof2nodeldof(int m) {
    return Array<int>({m/ndofn_, m%ndofn_});
  }

// Constitutive
 public:
  bool planeStressElsePlaneStrain_;
  
  // {E, \nu}
  Vectord YoungPoisson_x(double x0, double x1) const {
    return Vectord({200, 0.2});
  }
  
  // {\lambda, \mu}
  Vectord lameConsts_x(double x0, double x1) const {
    auto yp = YoungPoisson_x(x0, x1);
    if(planeStressElsePlaneStrain_) // plane stress
      return Vectord({yp(0)*yp(1) / (1.0 - yp(1)*yp(1)), yp(0) / (2*(1.0 + yp(1)))});
    else // plane strain
      return Vectord({yp(0)*yp(1) / ((1.0 + yp(1))*(1.0 - 2.0*yp(1))), yp(0) / (2.0*(1.0 + yp(1)))});
  }
 private:
  // St. Venant-Kirchoff
  Matrix2d stressFromStrain(const Matrix2d& strain, double x0, double x1) const {
    Matrix2d stress = Matrix2d(strain.nRows(), strain.nCols());
    const auto lambdaMu = lameConsts_x(x0, x1);
    double trE = strain.trace();
    for(int i=0; i<strain.nRows(); ++i)
      for(int j=0; j<strain.nCols(); ++j)
        stress(i, j) = lambdaMu(0)*diracDelta(i, j)*trE + 2.0*lambdaMu(1)*strain(i, j);
    return stress;
  }
  
  // St. Venant-Kirchoff tensor
  Matrix4d C_VK_x(double x0, double x1) {
    Matrix4d m(ndofn_, ndofn_, ndofn_, ndofn_);
    auto dD = LinAlg::diracDelta; // function alias
    
    // INEFFICIENT: C_VK inherently has a lot of symmetry and I'm not using that here, so a bit inefficient no doubt; but to use it, I really need symmetry-exploiting operators (if I can even use the operators rather than pure index notation)
    for(int i=0; i<ndofn_; ++i)
      for(int j=0; j<ndofn_; ++j)
        for(int k=0; k<ndofn_; ++k)
          for(int l=0; l<ndofn_; ++l)
            m(i,j,k,l) = lameConsts_x(x0, x1)(0) * dD(i, j) * dD(k, l) +
              lameConsts_x(x0, x1)(1) * (dD(i, k) * dD(j, l) + dD(i, l) * dD(j, k));
    return m;
  }
  inline Matrix4d C_VK_xi(double xi0, double xi1) {
    auto xVect = x_xi(xi0, xi1);
    return C_VK_x(xVect(0), xVect(1));
  }
  
// Kinematic
 public:
  enum strainMeasure {
    LINEAR,
    GREEN_LAGRANGE,
    EULER_ALMANSI,
    LOGARITHMIC
  };
  
  strainMeasure strainMeasure_ = LINEAR; // support linear first, then others
  
 private:
  // J_lj = \frac{\del x_j}{\del \xi_l}
  // i.e. gradR_x_wrtxi_xi
  Matrix2d jacobian_xi(double xi0, double xi1) {
    Matrix2d J(ndofn_, ndofn_);
    for(int l=0; l<ndofn_; ++l)
      for(int j=0; j<ndofn_; ++j) {
        double sumK = 0;
        for(int k=0; k<nnode_; ++k)
          sumK += gradL_shFct_wrtxi_xi(xi0, xi1)(k,l) * X_0_(nodeldof2dof(k, j));//#j or l in last l?
        J(l,j) = sumK;
      }
    return J;
  }
  // (J^{-1})_lj = \frac{\del \xi_l}{\del x_j}
  Matrix2d invJacobian_xi(double xi0, double xi1) {
    return LinAlg::invertMat2d(jacobian_xi(xi0, xi1));
  }

// ctor(s)
 public:
  Tri3(Array<int> nodes, Vectord X_0, bool planeStressElsePlaneStrain)
  : globalNodeIds_(nodes), X_0_(X_0), planeStressElsePlaneStrain_(planeStressElsePlaneStrain) {}
  
  Tri3(Vectord globalX_0, Array<int> nodes, bool planeStressElsePlaneStrain)
  : globalNodeIds_(nodes), planeStressElsePlaneStrain_(planeStressElsePlaneStrain) {
    X_0_ = Vectord();
    FOR(i, nodes.size()) {
      X_0_.push_back(globalX_0(nodes(i)*ndofn_));
      X_0_.push_back(globalX_0(nodes(i)*ndofn_+1));
    }
  }
  
// Shape functions
 private:
  ///static double lagrangePolynomialLin_xi(double xi, int i/*, int j for dim>1*/) {return (1.0/2) * (1.0 + Xi_[i]*xi);}
  
 private:
  static constexpr int shFct_xi_size_ = nnode_; // this and similar are ultimately unnecessary but I keep for now
  static inline Vectord shFct_xi(double xi0, double xi1) {
    return Vectord({1 - xi0 - xi1, xi0, xi1});
  }
  static constexpr int shFctMatrixForm_nRows_ = ndofn_;
  static constexpr int shFctMatrixForm_nCols_ = ndof_;
  /*///static inline Matrix2d shFctMatrixForm_xi(double xi0, double xi1) {
    Matrix2d mat(ndofn_, ndof_);
    auto shFct = shFct_xi(xi0, xi1);
    for(int i=0; i<ndofn_; ++i)
      for(int j=0; j<ndof_; ++j) {
        int nodeBlock = j/ndofn_; // floor
        mat(i, j) = diracDelta(i, j%ndofn_) * shFct(nodeBlock);
      }
    return mat;
  }*/
  // const because line shFcts
  static inline Matrix2d gradL_shFct_wrtxi_xi(double xi0 = -2, double xi = -2) {
    // \partial N_i / (\partial xi_j) (each row is grad(N_i) where i const)
    Matrix2d mat(3/*nnode_*/, 2/*ndofn_*/,
    {
      -1, -1, // grad shFct(0)
      1, 0, // grad shFct(1)
      0, 1 // grad shFct(2)
    }
    );
    return mat;
  }
  
  // row = node, col = node-local dof
  /*///inline Matrix2d X_0_asMat() const {
    Matrix2d mat(nnode_, ndofn_);
    for(int i=0; i<nnode_; ++i)
      for(int j=0; j<ndofn_; ++j)
        mat(i, j) = X_0_(i*ndofn_+j);
    return mat;
  }*/
  
  ///inline Matrix2d gradL_shFct_x(double x0, double x1) {
  ///  return mat2dTimesMat2d(invertMat2d(gradR_x_xi(x0, x1)), gradL_shFct_xi(x0, x1));
  ///}
  
  // \frac{\partial \mathrm N_k}{\partial x_j} = \sum_{l = 1}^{\text{ndofn}} \frac{\partial \mathrm N_k}{\partial \xi_l} (\bm J^{-1})_{lj} =: (\tilde{\bu N})_{kj}.
  // also \tilde{N}
  inline Matrix2d gradL_shFct_wrtx_xi(double xi0, double xi1) {
    Matrix2d mat(nnode_, ndofn_);
    auto invJ = invJacobian_xi(xi0, xi1);
    auto dNdxi = gradL_shFct_wrtxi_xi(xi0, xi1);
    for(int k=0; k<nnode_; ++k) {
      for(int j=0; j<ndofn_; ++j) {
        double tmpSum = 0;
        for(int l=0; l<ndofn_; ++l) {
          tmpSum += dNdxi(k,l) * invJ(l,j);
        mat(k,j) = tmpSum;
        }
      }
    }
    return mat;
  }
  
  // B operator, a 2x2x6 matrix B_{ijm}
  inline Matrix3d BOp_xi(double xi0, double xi1) {
    StandardTimer timer("BOp_xi()");
    timer.start();
    Matrix3d mat(ndofn_, ndofn_, ndof_);
    auto Ntilde = gradL_shFct_wrtx_xi(xi0, xi1);
    for(int i=0; i<ndofn_; ++i) {
      for(int j=0; j<ndofn_; ++j) {
        for(int m=0; m<ndof_; ++m) {
          auto vectm = dof2nodeldof(m);
          mat(i,j,m) = 0.5 * (Ntilde(vectm(0),j)*diracDelta(vectm(1),i) + Ntilde(vectm(0),i)*diracDelta(vectm(1),j));
        }
      }
    }
    timer.stop();
    return mat;
  }
  
  // Q and q will be my in general "some quantity"
  ///inline double q_xi(const Vectord& Q_xi_, double xi) {return vectDotProduct(shFct_xi(xi), Q_xi_);}
  
  inline Vectord x_xi(double xi0, double xi1) {
    Vectord x(ndofn_);
    FOR(i, ndofn_) {
      double tmpSum = 0;
      FOR(k, nnode_)
        tmpSum += shFct_xi(xi0, xi1)(k) * X_0_(nodeldof2dof(k,i));
      x(i) = tmpSum;
    }
    return x;
  }
  
  ///inline static double deriv(double (*f)(double), double xix, double h) {return (f(xix+h)-f(xix-h))/2.0;}
  
  
  static double integrateTrapz(double (*f)(double), double l, double r, double n) {
    double h = (r - l) / n;
    double result = 0;
    for(int i=0; i<n-1; ++i) {
      result += (f(i*h) + f((i+1)*h))/2.0 * h;
      //db::pr("result ="+std::to_string(result));
    }
    return result;
  }
  template<typename Fct>
  static double integrateTrapzTemplated(Fct f, double l, double r, double n) {
    double h = (r - l) / n;
    double result = 0;
    for(int i=0; i<n-1; ++i) {
      result += (f(i*h) + f((i+1)*h))/2.0 * h;
    }
    return result;
  }
  template<typename Fct>
  static double integrateGaussianQuadratureTemplated(Fct f, double l, double r, double n) {
    double result = 0;
    // TODO
    return result;
  }
  
  //\del E : S | For engineering/infinitesimal/linear strains: E = (\Grad u + (\Grad u)^T)/2 == makeMatSymmetric(\Grad u)
  //\del (1+F) : S(E)
  //(\del du/dx)*Cvk*(du/dx) | u = ND, N is a horizontally-long matrix, D is a vertical vector
  //\del D^T dN^T/dx*Cvk* dN/dx D | dN/dx = dN/dxi * dxi/dx = (+- 1/2) * 1/(dx/dxi); x = NX -> dx/dxi = dN/xi * X = (+- 1/2)*X --> dN/dx = (+- 1/2) * 1/((+- 1/2)*X)
  //           ^now a vertical vector                                                                               "  ^=:J   "                          "  ^=J    "
  // Integrating over xi: dx = J*dxi
  // --> \del D^T int_xi(dN^T/dx*Cvk* dN/dx * J dxi) * D
  // = \del D^T int_xi((deriv_shFct_xi^T*{J^-1}) * (deriv_shFct_xi*{J^-1}) * J dxi) * D * Cvk
  // = \del D^T int_xi(deriv_shFct_xi^T * deriv_shFct_xi dxi) * Jinv * D * Cvk // Jinv const so can be taken out of integral
  //                                    ^dyadic product because vertical vector times horizontal vector
  // Also, for Cvk(x) to Cvk(xi), we need x(xi). We have x(xi), it is just N(xi)*X
  Matrix2d integrandKmat_xi(double xi0, double xi1) {
    Matrix2d mat(ndof_, ndof_);
    double detJ = detMat2d(jacobian_xi(xi0, xi1));
    if(detJ<0) db::pr("detJ is negative! It is "+std::to_string(detJ));
    //db::pr("detJ="+std::to_string(detJ));
    auto CVK = C_VK_xi(xi0, xi1);
    auto B = BOp_xi(xi0, xi1);
    // TODO: inefficient as fuck
    FOR(mLeft, ndof_) {
      FOR(mRight, ndof_) {
        double tmpSum = 0;
        FOR(j, ndofn_) {
          FOR(i, ndofn_) {
            FOR(k, ndofn_) {
              FOR(l, ndofn_) {
                tmpSum += B(j,i,mLeft) * CVK(i,j,k,l) * B(k,l,mRight);
              }
            }
          }
        }
        mat(mLeft, mRight) = tmpSum * detJ;
      }
    }
    return mat;
  }
 public: //tmp
  Matrix2d Kmat() {
    std::string outString = "Element ";
    FOR(i, globalNodeIds_.size()-1)
      outString += std::to_string(globalNodeIds_(i)) + "-";
    outString += std::to_string(globalNodeIds_(globalNodeIds_.size()-1));
    outString += ": getting Kmat()\n";
    std::cout<<outString;
    
    Matrix2d result(ndof_, ndof_);
    
    // integration counts in xi0 and xi1 direction
    int n0 = 5;
    int n1 = 5;
    
    for(int i=0;i<ndof_; ++i)
      for(int j=0;j<ndof_; ++j) {
        ///auto integrandKmat_xi_ij = [this, i, j](double xi0, double xi1) {return integrandKmat_xi(xi0, xi1)(i, j);};
        
        auto g = [&, i, j, n1](double xi0)
        {
          auto f = [&, i, j, xi0](double xi1)
          {
            return integrandKmat_xi(xi0, xi1)(i, j);
          };
          
          // inner integration
          return integrateTrapzTemplated(f, 0.0, 1.0 - xi0, n1);
        };
        
        // outer integration
        result(i, j) = integrateTrapzTemplated(g, 0.0, 1.0, n0);

      }
    return result;
  }
  
  void test() {
    std::string outString = "--------------- Element ";
    FOR(i, globalNodeIds_.size()-1)
      outString += std::to_string(globalNodeIds_(i)) + "-";
    outString += std::to_string(globalNodeIds_(globalNodeIds_.size()-1));
    outString += " test() ---------------\n";
    
    std::cout<<outString;
    
    std::cout<<Kmat().toString(8)<<"\n";
    //gradL_shFct_wrtx_xi(-1, -1).print();
    
    //x_xi(0, 0.5).print();
  }
}; // class Tri3

} // namespace Element

} // namespace MyFem
