#include "Linear2D_Manager.hpp"

#include <SFML/Graphics.hpp> // TODO: delete

namespace MyFem::Problem {

Vectord Linear2D::compute_X_t_single(vector<size_t> dirichletDofIds, Vectord dirichletVals, size_t n) {
  if(n>=X_t_.size()) X_t_.resize(n+1);
  X_t_[n] = get_X_0();
  
  if(n>0) {
    Vectord U_t_full = Vectord(globalDofIds_.size());
    
    std::vector<size_t> solutionDofIds(globalDofIds_);
    StdVectorUtils::deleteIndices(solutionDofIds, dirichletDofIds);
    //std::cout<<"solutionDofIds.print();line11\n";
    //StdVectorUtils::print(solutionDofIds);
    FOR(i, solutionDofIds.size()) {
      U_t_full(solutionDofIds[i]) = U_t_[n](i);
    }
    
    FOR(i, dirichletDofIds.size()) {
      U_t_full(dirichletDofIds[i]) = dirichletVals(i);
    }
    
    FOR(i, X_t_[n].size()) {
      X_t_[n](i) += U_t_full(i);
    }
  }
  
  return X_t_[n];
}

Matrix2d Linear2D::assembleKfull() const {
  MyUtils::Timers::ScopedTimer timer("assembleKfull()");
  
  Matrix2d assembledK(nnode_*ndofn_, nnode_*ndofn_); // TODO: we use a dyn matrix because we want to get rid of the first loop and just have one loop later, but I have to see how I can do that
  
  STATUS("assembleKfull(): starting assembly...");
  MyUtils::Db::LoadingBar lb(elements_.size(), 100);
  for(int e=0; e<elements_.size(); ++e) {
    const auto& eleGlobalDofIds = elements_[e]->getGlobalDofIds();
    const auto locK = elements_[e]->Kmat();
    for(int i=0; i<locK.nRows(); ++i)
      for(int j=0; j<locK.nCols(); ++j)
        assembledK(eleGlobalDofIds[i],eleGlobalDofIds[j]) += locK(i, j);
    lb();
  }
  
  return assembledK;
  
  //Db::pr("assembledK:");
  //assembledK.print();
}

Matrix2d Linear2D::assembleMfull() const {
  MyUtils::Timers::ScopedTimer timer("assembleMfull()");
  
  Matrix2d assembledM(nnode_*ndofn_, nnode_*ndofn_);
  
  STATUS("assembleMfull(): starting assembly...");
  MyUtils::Db::LoadingBar lb(elements_.size(), 100);
  for(int e=0; e<elements_.size(); ++e) {
    const auto& eleGlobalDofIds = elements_[e]->getGlobalDofIds();
    const auto locK = elements_[e]->Mmat();
    for(int i=0; i<locK.nRows(); ++i)
      for(int j=0; j<locK.nCols(); ++j)
        assembledM(eleGlobalDofIds[i],eleGlobalDofIds[j]) += locK(i, j);
    lb();
  }
  
  return assembledM;
}

Vectord Linear2D::assembleFGravity(double gravityAccel) const {
  if(gravityAccel==0)
    return Vectord(get_ndof(), 0.0);
  else {
    MyUtils::Timers::ScopedTimer timer("assembleFGravity()");

    Vectord assembledF(get_ndof());

    STATUS("assembleFGravity(): starting assembly...");
    MyUtils::Db::LoadingBar lb(elements_.size(), 100);
    for(int e = 0; e < elements_.size(); ++e) {
      const auto& eleGlobalDofIds = elements_[e]->getGlobalDofIds();
      const auto locF = elements_[e]->getFGravity(gravityAccel);

      for(int i = 0; i < locF.size(); ++i)
        assembledF(eleGlobalDofIds[i]) += locF(i);

      lb();
    }

    return assembledF;
  }
}

void Linear2D::setNodeDofInfo() {
  MyUtils::Timers::ScopedTimer timer("setNodeDofInfo()");
  
  ndofn_ = elements_[0]->ndofn_;
  
  STATUS("setNodeDofInfo(): getting dof ids");
  
  int globalNodeCount = 0;
  for(int e=0; e<elements_.size(); ++e) {
    auto eleGlobalNodeIds = elements_[e]->globalNodeIds_;
    for(int locId=0; locId<eleGlobalNodeIds.size(); ++locId)
      if(StdVectorUtils::find(globalDofIds_, ndofn_*eleGlobalNodeIds[locId]).size()==0)
        if(eleGlobalNodeIds[locId]+1 > globalNodeCount)
          globalNodeCount = eleGlobalNodeIds[locId]+1;
  }
  FOR(i, globalNodeCount) {
        globalDofIds_.push_back(ndofn_*i);
        globalDofIds_.push_back(ndofn_*i + 1);
  }
  
  nnode_ = globalNodeCount;
}

void Linear2D::assembleAll(double gravityAccel) {
  MyUtils::Timers::ScopedTimer timer("assembleFull()");
  
  setNodeDofInfo();
  KFull_ = assembleKfull();
  MFull_ = assembleMfull();
  rhsFull_ = assembleFGravity(gravityAccel);
  
  if(KFull_.size() != MFull_.size() || KFull_.size() != rhsFull_.size()) THROW("assembleAll(): KFull_.size() != MFull_.size() || KFull_.size() != rhsFull_.size()");
}

// static
void Linear2D::removeDuplicates(vector<size_t>& ids, Vectord& vals) {
  for(size_t i = ids.size(); i-- > 0; ) {
    if(StdVectorUtils::find(ids, ids[i]).size()>1) {
      StdVectorUtils::deleteIndices(ids, {i});
      vals.deleteIndices({i});
    }
  }
}

void Linear2D::applyNeumannToRhsFull(Vectord& rhsFull, std::vector<size_t>& ids, Vectord& vals) const {
  STATUS("applyNeumannToRhsFull()");
  MyUtils::Timers::ScopedTimer timer("applyNeumannToRhsFull()");
  
  removeDuplicates(ids, vals);
  
  FOR(i, ids.size())
    rhsFull(ids[i]) += vals(i);
}

void Linear2D::setFreeAndDirichDofIds(std::vector<size_t>& dirichDofIds, Vectord& dirichVals) {
  STATUS("setFreeAndDirichDofIds()");
  MyUtils::Timers::ScopedTimer timer("setFreeAndDirichDofIds()");
  
  removeDuplicates(dirichDofIds, dirichVals);
  
  dirichDofIds_ = dirichDofIds;
  
  freeDofIds_ = globalDofIds_;
  StdVectorUtils::deleteIndices(freeDofIds_, dirichDofIds); // or technically deleteIndices(...find(...)) but it will be the same
}

// static
Matrix2d Linear2D::reducedMatrix(const Matrix2d& A, const vector<size_t>& freeDofIds) {
  STATUS("reducedMatrix()");
  
  int n_reduced = freeDofIds.size();
  Matrix2d A_reduced(n_reduced, n_reduced);
  
  FOR(v, n_reduced) {
    int idS = freeDofIds[v];
    
    FOR(v2, n_reduced) {
      int id2 = freeDofIds[v2];
      A_reduced(v, v2) = A(idS, id2);
    }
  }
  
  return A_reduced;
}

// static
Vectord Linear2D::reducedVector(const Vectord& x, const vector<size_t>& freeDofIds) {
  STATUS("reducedVector()");
  
  int n_reduced = freeDofIds.size();
  Vectord x_reduced(n_reduced, n_reduced);
  
  FOR(v, n_reduced) {
    x_reduced(v) = x(freeDofIds[v]);
  }
  
  return x_reduced;
}

void Linear2D::applyDirichletToRhs(Vectord& rhsFull, const Matrix2d& KFull, vector<size_t>& dirichletDofIds, Vectord& dirichletVals) const {
  STATUS("applyDirichletToRhs()");
  
  int n_reduced = freeDofIds_.size();
  Vectord rhs_reduced(n_reduced);
  
  removeDuplicates(dirichletDofIds, dirichletVals);
  
  FOR(v, n_reduced) {
    int idS = freeDofIds_[v];
    
    double rhs_i = rhsFull(idS);
    
    FOR(i, dirichletDofIds.size()) {
      int dId = dirichletDofIds[i];
      rhs_i -= KFull(idS, dId) * dirichletVals(i);
    }
    rhs_reduced(v) = rhs_i;
  }
  
  rhsFull = rhs_reduced;
}

/*void Linear2D::applyDirichlet(vector<size_t>& ids, Vectord& vals){
  setFreeAndDirichDofIds(ids, vals);
  
  rhsDirich_ = rhsFull_;
  applyDirichletToRhs(rhsDirich_, ids, vals);
  
  KRed_ = reducedMatrix(KFull_, freeDofIds_);
  MRed_ = reducedMatrix(MFull_, freeDofIds_);
  ///CRed_ = reducedMatrix(CFull_, freeDofIds_);
}*/

void Linear2D::readMeshTxt(std::string inputFilePath) {
  std::vector<std::string> fileData;
  MyUtils::IO::readFileLines(inputFilePath, fileData);
  MyUtils::Db::pr("11111");
  
  for(int i = 0; i<fileData.size(); i++) {
    std::string current = fileData.at(i);
    int* cFI = MyUtils::Strings::checkForIn("//", current, 1);
    if(cFI[0]==0) {
      fileData.erase(fileData.begin() + i);
      i--;
    }
    
    free(cFI);
  }
  for(int i = 0; i<fileData.size(); i++) {
    std::string current = fileData.at(i);
    // The level is how many tabs to the right does the line text start (one tab = 2 spaces; actual '\t' characters are not supported)
    int whiteSpaceEndPos = MyUtils::Strings::getEndOfWhitespace(current);
    int numSpacesPerTab = 2; // But we can also make this an argument if we need to
    int level = whiteSpaceEndPos / numSpacesPerTab;
    MyUtils::Db::pr("line"+std::to_string(i)+",level="+std::to_string(level));
    if(current.length()==whiteSpaceEndPos) { // i.e. there is only white space in this line
      fileData.erase(fileData.begin() + i);
      i--;
    }
    //else if(MyUtils::Strings::keepInterval(current,))
  }
  MyUtils::IO::writeFileLines("test", fileData);
  
  //delete(fileData);
}

std::string Linear2D::infoString() {
  using std::string;
  using std::to_string;
  
  string str = "-- Problem Information --\n";
  str+="ndof="+to_string(get_ndof())+"\n";
  str+="nnode="+to_string(get_nnode())+"\n";
  str+="ndofn="+to_string(get_ndofn())+"\n";
  str+="nele="+to_string(elements_.size())+"\n";
  //str+="\n";
  //str+="num dirichlet dofs="+to_string(vals.size())+"\n";
  //str+="num neumann dofs="+to_string(neumannVect_.size())+"\n";
  
  
  return str;
}
void Linear2D::printInfo() {
  std::cout<<infoString()<<"\n";
}

void Linear2D::example_beam_dyn(double lx, double ly, int nx, int ny, int maxIter, double tol, const double density) {
  STATUS("Running example_beam_dyn()");
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
    return Vectord(vector<double>{10, 0.2});
  };
  auto density_x = [=](double x0, double x1) {
    return double(1);
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
  printInfo();
  auto KFull = assembleKfull();
  auto MFull = assembleMfull();
  auto FGravity = assembleFGravity(0);
  if(KFull.nRows() != MFull.nRows() || KFull.nRows() != FGravity.size()) THROW("KFull_.size() != MFull_.size() || KFull_.size() != rhsFull_.size()");
  
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
  
  int sysSizeReduced = KRed.nRows();
  
  
  // All of these will be in reduced form (reduced size)
  Vectord U_last;
  Vectord V_last;
  Vectord A_last;
  
  // Small deformation dynamic problem time integration; no change in dirichlet currently
  // Newmark method
  double beta = 1.0/4;
  double gamma = 1.0/2;
  
  finalT_ = 5;
  deltaT_ = 0.1;
  int timeSteps = get_timeSteps();
  
  double deltaT_2_ = deltaT_*deltaT_;
  
  for(int n = 0; n<=timeSteps; ++n) {
    double currT = n*deltaT_;
    STATUS("n="+std::to_string(n)+"; t="+std::to_string(currT))
  
    // //Neumann BC
    auto neumannIds = std::vector<size_t>();
    auto neumannVect = Vectord();
    FOR(j, ny) {
      // from 0 to 1
      double yFrac = (double)j/(ny-1);
      // right side
      
      // horizontal
      /*neumannIds.push_back(ndofn_*(nx*(j+1)-1)+ 0);
      if(currT>1.5&&currT<1.6) {
        neumannVect.push_back(0.1);//*1*0.5*yFrac(1-yFrac));
      }
      else {
        neumannVect.push_back(0);//*1*0.5*yFrac(1-yFrac));
      }*/
      
      // vertical
      neumannIds.push_back(ndofn_*(nx*(j+1)-1)+ 1);
      if(currT>0.4&&currT<0.7) {
        neumannVect.push_back(0.1*deltaT_);//*currT);//*1*0.5*yFrac(1-yFrac));
      }
      else {
        neumannVect.push_back(0);//*1*0.5*yFrac(1-yFrac));
      }
    }
    auto rhs_curr = FGravity;
    applyNeumannToRhsFull(rhs_curr, neumannIds, neumannVect);
    applyDirichletToRhs(rhs_curr, KFull, dirichIds, dirichVect);
    
    if(n==0) {
      U_last = Vectord(sysSizeReduced, 0.0);
      V_last = reducedVector(V_0, freeDofIds_);  
      // Initial acceleration field; get this by solving M * A_0 = F_0 - K * X_0 // (all reduced)
      A_last = MyUtils::NumMethods::LinSolvers::GaussSeidel(MRed, vectdPlusVectd(rhs_curr, scaleVectd(-1, mat2dTimesVectd(KRed, U_last))), maxIter, tol);
      
      U_t_.push_back(U_last);
      compute_X_t_single(dirichIds, dirichVect, n);
      
      continue;
    }
    
    
    Vectord U_predict = vectdPlusVectd(vectdPlusVectd(U_last, scaleVectd(deltaT_, V_last)), scaleVectd(deltaT_2_*(1.0/2-beta), A_last));
    Vectord V_predict = vectdPlusVectd(V_last, scaleVectd(deltaT_*(1.0-gamma), A_last));
  
    // Solve for A_curr: (M + beta * deltaT^2 * K) A_curr = rhs_curr - K * U_predict
    MyUtils::Timers::StandardTimer timerLinSolve("GaussSeidel()");timerLinSolve.start();
    Vectord A_curr = MyUtils::NumMethods::LinSolvers::GaussSeidel(mat2dPlusMat2d(MRed, scaleMat2d(beta*deltaT_2_, KRed)), vectdPlusVectd(rhs_curr, scaleVectd(-1, mat2dTimesVectd(KRed, U_predict))), maxIter, tol);
    timerLinSolve.stop();
    
    
    // U_curr for next iter
    U_last = vectdPlusVectd(U_predict, scaleVectd(beta*deltaT_2_, A_curr));
    U_t_.push_back(U_last);
    // V_curr for next iter
    V_last = vectdPlusVectd(V_predict, scaleVectd(gamma*deltaT_, A_curr));
    
    A_last = A_curr;
    
    compute_X_t_single(dirichIds, dirichVect, n);
  }
  
  printInfo();
}

} // namespace MyFem::Problem
