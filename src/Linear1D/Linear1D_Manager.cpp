#include "Linear1D_Manager.hpp"

namespace Problem {
  
void Linear1D::runNoInputExample() {
  //double nodeX_0[] = {0, 0.7, 1, 1.5, 1.9};
  
  dynArrayd nodeX_0(20);
  
  for(int i=0; i<nodeX_0.size(); ++i)
    nodeX_0[i] = 0.1*i;
  
  /*Element::line2* elements[] = {
    new Element::line2({0, 1}, {nodeX_0[0], nodeX_0[1]}),
    new Element::line2({1, 2}, {nodeX_0[1], nodeX_0[2]})
  };*/
  db::pr("line 15");
  for(int i=0; i<nodeX_0.size()-1; ++i)
    elements_.push_back(new Element::line2({i, i+1}, {nodeX_0[i], nodeX_0[i+1]}));
  
  
  assembleK();
  
  rhs_ = dynArrayd(K_.size()); // zero vect for now
  
  applyDirichlet(nodeX_0.size()-1, 0.1);
  db::pr("rhs_ after1");
  LinAlg::Dynamic::printVectd(rhs_);
  applyDirichlet(0, 0.0);
  
  db::pr("rhs_ after2");
  LinAlg::Dynamic::printVectd(rhs_);
  
  db::pr("K_ after apply dirichlet");
  LinAlg::Dynamic::printMatd(K_);
  
  solveSystem_Jacobi();
  
  /*for(int i=0; i<2; ++i)
    for(int j=0; j<2; ++j) {
      db::pr("line 19, i="+std::to_string(i)+", j="+std::to_string(j));
      std::cout<<"elements[0]->Kmat()[i][j] = "<<elements[0]->Kmat()[i][j]<<"\n";
    }
    
  //Element_line2::test();
  */
  free(elements_[0]);
  free(elements_[1]);
}

dynMatrixd Linear1D::assembleK() {
  int globalNodeCount = 0;
  for(int e=0; e<elements_.size(); ++e) {
    auto eleGlobalNodeIds = elements_[e]->globalNodeIds_;
    for(int locId=0; locId<eleGlobalNodeIds.size(); ++locId)
      if(eleGlobalNodeIds[locId]+1 > globalNodeCount)
        globalNodeCount = eleGlobalNodeIds[locId]+1;
  }
  dynMatrixd assembledK(globalNodeCount); // TODO we use a dyn matrix because we want to get rid of the first loop and just have one loop later, but I have to see how I can do that
  // TODO in general I need to make better matrix type
  db::pr("line41");
  for(int i = 0; i<globalNodeCount; ++i)
    assembledK[i] = dynArrayd(globalNodeCount);
  db::pr("line44");
  db::pr(std::to_string(globalNodeCount));
  for(int e=0; e<elements_.size(); ++e) {
    const auto& eleGlobalNodeIds = elements_[e]->globalNodeIds_;
    db::pr("line47");
    const auto locK = elements_[e]->Kmat();
    for(int i=0; i<locK.size(); ++i)
      for(int j=0; j<locK[0].size(); ++j)
        assembledK[eleGlobalNodeIds[i]][eleGlobalNodeIds[j]] += locK[i][j];
  }
  db::pr("assembledK =\n");
  LinAlg::Dynamic::printMatd(assembledK);
  
  K_ = assembledK;
  return assembledK;
}

void Linear1D::applyDirichlet(int globalNodeId/*, dofIndex or localDofindex of the node, needed for higher dim*/, double val) {
  // x_d==val means that any occurence of K_ij*x_d, ie when j==d, should be instead subtracted from the rhs. and if i==d, it doesnt matter anyays because we remove that row or practically remove it.
  // ok but actually we want the following approach: remove the rows of dirichlet nodes but then subtract from the rhs the upper right*x_d. ok actually maybe practically the same thing, which makes sense I guess
  int n = K_.size();
  // TODO very inefficient but whatever
  
  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      if(j==globalNodeId && K_[i][j]!=0){
        rhs_[i] -= K_[i][j]*val;
        db::pr("line89, ij="+std::to_string(i)+std::to_string(j));
        db::pr(std::to_string(K_[i][j]*val));
      }
    }
  }
  
  db::pr("rhs_");
  LinAlg::Dynamic::printVectd(rhs_);
  
  
  auto K2 = dynMatrixd(n-1, dynArrayd(n-1));
  auto rhs2 = dynArrayd(n-1);
  int i2 = 0;
  for(int i=0; i<n; ++i) {
    if(i!=globalNodeId) {
      int j2 = 0;
      for(int j=0; j<n; ++j) {
        if(j!=globalNodeId){
          K2[i2][j2] = K_[i][j];
          j2++;
        }
      }
      
      rhs2[i2] = rhs_[i];
      
      i2++;
    }
  }
  
  K_ = K2;
  rhs_ = rhs2;
}

dynArrayd Linear1D::solveSystem_gaussSeidel() {
  using namespace LinAlg::Dynamic;
  const int n = K_.size();
  auto L = K_;
  auto U = K_;
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j) {
      if(j>i)
        L[i][j] = 0;
      else//if(i>=j)
        U[i][j] = 0;
    }
    
  auto residual = [K = K_, rhs = rhs_](dynArrayd x){return vectAdd(matTimesVect(K, x), scaleVect(-1.0, rhs));};
  auto residualNorm = [&, residual](dynArrayd x){auto res = residual(x); return vectDotProduct(res, res);};
  dynArrayd x_0(n); // Initial guess, zeros here
  dynArrayd x_i = x_0;
  int iter = 0;
  constexpr int maxiter = 200;
  constexpr double maxResNorm = 0.01;
  while(iter < maxiter/* && residualNorm(x_i) < maxResNorm*/) {
    auto bMinusUx = vectAdd(rhs_, scaleVect(-1.0, matTimesVect(U, x_i)));
    //x_i =
  }
  return x_i;
}
dynArrayd Linear1D::solveSystem_Jacobi() {
  using namespace LinAlg::Dynamic;
  db::pr("Solve start");
  printMatd(K_);
  
  const int n = K_.size();
  auto L = K_;
  auto U = K_;
  auto D = K_;
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j) {
      if(j>=i)
        L[i][j] = 0;
      if(i>=j)
        U[i][j] = 0;
      if(i!=j)
        D[i][j] = 0;
    }
    
  auto residual = [K = K_, rhs = rhs_](dynArrayd x){return vectAdd(matTimesVect(K, x), scaleVect(-1.0, rhs));};
  auto residualNorm = [&, residual](dynArrayd x){auto res = residual(x); return vectDotProduct(res, res);};
  dynArrayd x_0(n); // Initial guess, zeros here
  dynArrayd x_i = x_0;
  int iter = 0;
  constexpr int maxiter = 50;
  constexpr double maxResNorm = 0.01;
  
  db::pr("K_ L U D");
  printMatd(K_);
  printMatd(L);
  printMatd(U);
  printMatd(D);
  
  db::pr("matAdd(L, U)");
  printMatd(matAdd(L, U));
  
  auto invD = D;
  for(int i=0; i<n; ++i)
    invD[i][i] = 1.0/D[i][i];
  db::pr("line 112");
  while(iter < maxiter/* && residualNorm(x_i) < maxResNorm*/) {
    x_i = matTimesVect(invD,
      vectAdd(rhs_,
        scaleVect(-1.0,
          matTimesVect(matAdd(L, U),
            x_i))));
    iter++;
    printVectd(x_i);
    std::cout<<residualNorm(x_i)<<std::endl;
  }
  
  return x_i;
  
}

void Linear1D::readMeshTxt(std::string inputFilePath) {
  std::vector<std::string>* fileData = new std::vector<std::string>();
  Utils::IO::readFileLines(inputFilePath, fileData, 5000);
  db::pr("11111");
  
  for(int i = 0; i<fileData->size(); i++) {
    std::string current = fileData->at(i);
    int* cFI = Utils::Strings::checkForIn("//", current, 1);
    if(cFI[0]==0) {
      fileData->erase(fileData->begin() + i);
      i--;
    }
    
    free(cFI);
  }
  for(int i = 0; i<fileData->size(); i++) {
    std::string current = fileData->at(i);
    // The level is how many tabs to the right does the line text start (one tab = 2 spaces; actual '\t' characters are not supported)
    int whiteSpaceEndPos = Utils::Strings::getEndOfWhitespace(current);
    int numSpacesPerTab = 2; // But we can also make this an argument if we need to
    int level = whiteSpaceEndPos / numSpacesPerTab;
    db::pr("line"+std::to_string(i)+",level="+std::to_string(level));
    if(current.length()==whiteSpaceEndPos) { // i.e. there is only white space in this line
      fileData->erase(fileData->begin() + i);
      i--;
    }
    //else if(Utils::Strings::keepInterval(current,))
  }
  Utils::IO::writeFileLines("test", fileData);
  
  //delete(fileData);
}

}