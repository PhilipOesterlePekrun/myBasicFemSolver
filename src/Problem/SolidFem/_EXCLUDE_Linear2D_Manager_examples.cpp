#include "Linear2D_Manager.hpp"

#include <SFML/Graphics.hpp> // TODO: delete

namespace MyFem::Problem {

void Linear2D::example_beam(double lx, double ly, int nx, int ny, int maxIter, double tol, const double density) {
  STATUS("Running example_beam()");
  printInfo();
  MyUtils::Timers::ScopedTimer timer("example_beam()");
  
  MyUtils::Timers::StandardTimer timer2("Example meshing");
  timer2.start();
  
  X_0_ = Vectord();
  
  FOR(j, ny)
    FOR(i, nx) {
      double xNom = i*lx/(nx-1);
      double yNom = j*ly/(ny-1);
      double x = xNom;
      double y = yNom;
      /*{
        double tf = 2.5*yNom;
        y = yNom+xNom*(tf-yNom)/lx; // asymmetric linear tapering
      }
      {
        MyUtils::Db::pr("yNom="+std::to_string(yNom));
        double tExtremum = 0.5*yNom;
        double tf = 0.9*yNom;
        double c = yNom;
        double b = NumMethods::solveScalarQuadraticEq(1, 4/lx*(yNom-tExtremum), -4/(lx*lx)*(tf-yNom)*(yNom-tExtremum))(1);
        MyUtils::Db::pr("b="+std::to_string(b));
        double a = (tf-yNom)/(lx*lx) - b/lx;
        double unitParabola = (a*(xNom*xNom)+b*xNom+c)/yNom;
        //y = yNom*unitParabola; // asymmetric parabola
      }*/
      X_0_.push_back(x);
      X_0_.push_back(y);
    }
    
  auto eleNodes = std::vector<std::vector<int>>();
  FOR(i, nx-1)
    FOR(j, ny-1) {
      vector<int> firstTri = {j*nx+ i, j*nx+ i+1, (j+1)*nx+ i};
      vector<int> secondTri = {j*nx+ i+1, (j+1)*nx+ i+1, (j+1)*nx+ i};
      
      //if(!(i>0&&i<nx-2 && j==1)) {
        eleNodes.push_back(firstTri);
        eleNodes.push_back(secondTri);
      }
    //}
    
  auto youngPoisson_x = [](double x0, double x1) {
    return Vectord(vector<double>{1, 0.2});
  };
  auto density_x = [=](double x0, double x1) {
    return double(density);
  };
    
  FOR(e, eleNodes.size()) {
    elements_.push_back(new Element::Tri3(
      X_0_,
      eleNodes[e],
      false,
      youngPoisson_x,
      density_x
    )
    );
  }
  
  timer2.stop();
  STATUS("Built elements");
  
  /*FOR(i, elements_.size()) {
    std::cout<<"================\nEle"<<i<<"\n";
    elements_(i)->test();
  }*/
  
  
  assembleAll(-9.81);
  //std::cout<<"K_ after assembly:\n";
  //std::cout<<K_.toString(8)<<"n";
  
  std::cout<<"rhs_:\n"<<rhsDirich_.toString();
  
  ///int numDofs = get_ndof();
  
  // //Neumann BC
  auto neumannIds = std::vector<size_t>();
  auto neumannVect = Vectord();
  
  FOR(j, ny) {
    // from 0 to 1
    double yFrac = (double)j/(ny-1);
    // right side
    //neumannIds.push_back(ndofn_*(nx*(j+1)-1)+ 0);
    //neumannVect.push_back(4*0.1*yFrac*(1-yFrac)); // second coeff is maximum
    //neumannIds.push_back(ndofn_*(nx*(j+1)-1)+ 1);
    //neumannVect.push_back(1*0.1*yFrac*(1-yFrac));
  }
  
  FOR(i, nx) {
    if(i==0||i==nx-1) continue;
    // from 0 to 1
    double yFrac = (double)i/(ny-1);
    // upper side
    //neumannIds.push_back(ndofn_*nx*(i+1)-1)+ 1);
    //neumannVect.push_back(4*0.05*yFrac*(1-yFrac));
  }
  
  
  // // Dirichlet BC
  auto dirichIds = std::vector<size_t>();
  auto dirichVect = Vectord();
  FOR(j, ny) {
    // left side
    dirichIds.push_back(ndofn_*nx*j+ 0);
    dirichVect.push_back(0.0);
    dirichIds.push_back(ndofn_*nx*j+ 1);
    dirichVect.push_back(0.0);
    
    
    /*dirichIds.push_back(ndofn_*nx*j+ 3);
    dirichVect.push_back(0.0);
    dirichIds.push_back(ndofn_*nx*j+ 4);
    dirichVect.push_back(0.0);*/
    
    
    // right side
    dirichIds.push_back(ndofn_*(nx*(j+1)-1)+ 0);
    dirichVect.push_back(0.0);
    dirichIds.push_back(ndofn_*(nx*(j+1)-1)+ 1);
    dirichVect.push_back(0.0);
  }
  
  FOR(i, nx) {
    // lower side
    /*dirichIds.push_back(ndofn_*i+ 0);
    dirichVect.push_back(0.0);
    dirichIds.push_back(ndofn_*i+ 1);
    dirichVect.push_back(0.0);*/
  }
  
  //dirichIds.push_back(ndofn_*(nx*(ny)-1)+ 0);
  //dirichVect.push_back(-2.0);
  
  /*dirichIds.push_back(ndofn_*(int)std::floor(nx/2.0)+ 1);
  dirichVect.push_back(-0.3);
  dirichIds.push_back(ndofn_*((int)std::floor(nx/2.0)-1)+ 1);
  dirichVect.push_back(-0.3);*/
  
  
  
  /*applyDirichlet(
    std::vector<size_t>{{
      ndofn_*0+ 0, ndofn_*0+ 1,
      ndofn_*nx*(ny-1)+ 0, ndofn_*nx*(ny-1)+ 1,
      ndofn_*(nx-1)+ 0, ndofn_*(nx-1)+ 1,
      ndofn_*(nx*ny-1)+ 0, ndofn_*nx*ny-1+ 1}},
    Vectord{{
      0.0, 0.0,
      0.0, 0.0,
      0.2, -0.1,
    0.2, -0.1}});
  */
  
  // remove duplicates
  /*for(size_t i=dirichIds.size()-1; i>0; --i) {
    if(StdVectorUtils::find(dirichIds, dirichIds[i]).size()>1) {
      StdVectorUtils::deleteIndices(dirichIds, {i});
      dirichVect.deleteIndices({i});
    }
  }*/
  
  applyBCs(neumannIds, neumannVect, dirichIds, dirichVect);
      
  
  //std::cout<<"K_ after applyDirichlet():\n";
  //K_.print();
  
  //std::cout<<"rhs_ after applyDirichlet():\n";
  //rhs_.print();
  
  {
    MyUtils::Timers::ScopedTimer timerLinSolve("GaussSeidel()");
    U_t_[1] = MyUtils::NumMethods::LinSolvers::GaussSeidel(KRed_, rhsDirich_, maxIter, tol, 1);
  }
  
  
  //solutionVect_.print(8);
  /*
  std::cout<<"\nComputation finished.\n\n";
  
  MyUtils::Db::pr("globalDofIds_:\n");
  StdVectorUtils::print(globalDofIds_);
  
  MyUtils::Db::pr("solutionVect_:\n");
  solutionVect_.print(8);
  MyUtils::Db::pr("dirichletDofIds_:\n");
  StdVectorUtils::print(dirichletDofIds_);
  MyUtils::Db::pr("dirichletVect_:\n");
  dirichletVect_.print(8);
  MyUtils::Db::pr("fullSolution():\n");
  fullSolution().print(8);
  
  MyUtils::Db::pr("\nX_0_:\n");
  X_0_.print(8);
  MyUtils::Db::pr("getX_t():\n");
  getX_t().print(8);
  */
  
  printInfo();
}

void Linear2D::example_beam_dyn(double lx, double ly, int nx, int ny, int maxIter, double tol, const double density) {
  STATUS("Running example_beam()");
  printInfo();
  MyUtils::Timers::ScopedTimer timer("example_beam_dyn()");
  
  MyUtils::Timers::StandardTimer timer2("Example meshing");
  timer2.start();
  
  X_0_ = Vectord(); // Initial displacement field (the mesh)
  auto V_0 = Vectord(); // Initial velocity field
  
  FOR(j, ny)
    FOR(i, nx) {
      double xNom = i*lx/(nx-1);
      double yNom = j*ly/(ny-1);
      double x = xNom;
      double y = yNom;
      X_0_.push_back(x);
      X_0_.push_back(y);
      
      V_0.push_back(0);
      V_0.push_back(0);
    }
    
  auto eleNodes = std::vector<std::vector<int>>();
  FOR(i, nx-1)
    FOR(j, ny-1) {
      vector<int> firstTri = {j*nx+ i, j*nx+ i+1, (j+1)*nx+ i};
      vector<int> secondTri = {j*nx+ i+1, (j+1)*nx+ i+1, (j+1)*nx+ i};
      
      //if(!(i>0&&i<nx-2 && j==1)) {
        eleNodes.push_back(firstTri);
        eleNodes.push_back(secondTri);
      }
    //}
    
  auto youngPoisson_x = [](double x0, double x1) {
    return Vectord(vector<double>{1, 0.2});
  };
  auto density_x = [=](double x0, double x1) {
    return double(density);
  };
    
  FOR(e, eleNodes.size()) {
    elements_.push_back(new Element::Tri3(
      X_0_,
      eleNodes[e],
      false,
      youngPoisson_x,
      density_x
    )
    );
  }
  
  timer2.stop();
  STATUS("Built elements");
  
  
  ///assembleAll(0);//-9.81);
  setNodeDofInfo();
  auto KFull = assembleKfull();
  auto MFull = assembleMfull();
  auto FGravity = assembleFGravity(0);
  if(KFull_.size() != MFull_.size() || KFull_.size() != rhsFull_.size()) THROW("KFull_.size() != MFull_.size() || KFull_.size() != rhsFull_.size()");
  
  // // Dirichlet BC
  auto dirichIds = std::vector<size_t>();
  auto dirichVect = Vectord();
  FOR(j, ny) {
    // left side
    dirichIds.push_back(ndofn_*nx*j+ 0);
    dirichVect.push_back(0.0);
    dirichIds.push_back(ndofn_*nx*j+ 1);
    dirichVect.push_back(0.0);
    
    
    // right side
    //dirichIds.push_back(ndofn_*(nx*(j+1)-1)+ 0);
    //dirichVect.push_back(0.0);
    //dirichIds.push_back(ndofn_*(nx*(j+1)-1)+ 1);
    //dirichVect.push_back(0.0);
  }
  //applyDirichlet(dirichIds, dirichVect);
  
  setFreeAndDirichDofIds(dirichIds, dirichVect);
  auto KRed = reducedMatrix(KFull, freeDofIds_);
  auto MRed = reducedMatrix(MFull, freeDofIds_);
  auto X_0Red = reducedVector(X_0_, freeDofIds_);
  
  int sysSizeReduced = KRed_.size();
  
  
  // All of these will be in reduced form (reduced size)
  Vectord U_last;
  Vectord V_last;
  Vectord A_last;
  
  // Small deformation dynamic problem time integration; no change in dirichlet currently
  // Newmark method
  double beta = 1.0/4;
  double gamma = 1.0/2;
  
  double finalT_ = 10;
  double deltaT_ = 0.05;
  int timeSteps = get_timeSteps();
  
  double deltaT_2 = deltaT_*deltaT_;
  
  for(int n = 0; n<=timeSteps; ++n) {
  
    // //Neumann BC
    auto neumannIds = std::vector<size_t>();
    auto neumannVect = Vectord();
    FOR(j, ny) {
      // from 0 to 1
      double yFrac = (double)j/(ny-1);
      // right side
      //neumannIds.push_back(ndofn_*(nx*(j+1)-1)+ 0);
      //neumannVect.push_back(4*0.1*yFrac*(1-yFrac)); // second coeff is maximum
      neumannIds.push_back(ndofn_*(nx*(j+1)-1)+ 1);
      neumannVect.push_back(deltaT_*timeSteps*1*0.1*yFrac*(1-yFrac));
    }
    auto rhs_curr = FGravity;
    applyNeumannToRhsFull(rhs_curr, neumannIds, neumannVect);
    applyDirichletToRhs(rhs_curr, dirichIds, dirichVect);
    
    if(n==0) {
      auto U_last = Vectord(sysSizeReduced, 0.0);
      auto V_last = reducedVector(V_0, freeDofIds_);  
      // Initial acceleration field; get this by solving M * A_0 = F_0 - K * X_0 // (all reduced)
      Vectord A_last = MyUtils::NumMethods::LinSolvers::GaussSeidel(MRed, vectdPlusVectd(rhs_curr, scaleVectd(-1, mat2dTimesVectd(KRed, X_0Red))), maxIter, tol);
      
      U_t_[n] = U_last;
      compute_X_t_single(dirichIds, dirichVect, n);
      
      continue;
    }
    
    
    Vectord U_predict = vectdPlusVectd(vectdPlusVectd(U_last, scaleVectd(deltaT_, V_last)), scaleVectd(deltaT_2*(1.0/2-beta), A_last));
    Vectord V_predict = vectdPlusVectd(V_last, scaleVectd(deltaT_*(1.0-gamma), A_last));
    
    // Solve for A_curr: (M + beta * deltaT^2 * K) A_curr = rhs_curr - K * U_predict
    MyUtils::Timers::StandardTimer timerLinSolve("GaussSeidel()");timerLinSolve.start();
    Vectord A_curr = MyUtils::NumMethods::LinSolvers::GaussSeidel(mat2dPlusMat2d(MRed, scaleMat2d(beta*deltaT_2, KRed)), vectdPlusVectd(rhs_curr, scaleVectd(-1, mat2dTimesVectd(KRed, U_predict))), maxIter, tol);
    timerLinSolve.stop();
    
    // U_curr for next iter
    U_t_[n] = U_last = vectdPlusVectd(U_predict, scaleVectd(beta*deltaT_2, A_curr));
    // V_curr for next iter
    V_last = vectdPlusVectd(V_predict, scaleVectd(gamma*deltaT_, A_curr));
    
    compute_X_t_single(dirichIds, dirichVect, n);
  }
  
  printInfo();
}

} // namespace MyFem::Problem
