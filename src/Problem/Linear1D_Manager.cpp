#include "Linear1D_Manager.hpp"

namespace Problem {
  
void Linear1D::runNoInputExample() {
  //double nodeX_0[] = {0, 0.7, 1, 1.5, 1.9};
  
  Vectord nodeX_0(20);
  
  for(int i=0; i<nodeX_0.size(); ++i)
    nodeX_0(i) = 0.1*i;
  
  /*Element::line2* elements[] = {
    new Element::line2({0, 1}, {nodeX_0[0], nodeX_0[1]}),
    new Element::line2({1, 2}, {nodeX_0[1], nodeX_0[2]})
  };*/
  db::pr("line 15");
  for(int i=0; i<nodeX_0.size()-1; ++i)
    elements_.push_back(new Element::Line2(Array<int>({i, i+1}), Vectord({nodeX_0(i), nodeX_0(i+1)})));
  
  
  assembleK();
  
  rhs_ = Vectord(K_.nRows()); // zero vect for now
  
  applyDirichlet(nodeX_0.size()-1, 0.1);
  db::pr("rhs_ after1");
  rhs_.print();
  applyDirichlet(0, 0.0);
  
  db::pr("rhs_ after2");
  rhs_.print();
  
  solveSystem_Jacobi();
  
  /*for(int i=0; i<2; ++i)
    for(int j=0; j<2; ++j) {
      db::pr("line 19, i="+std::to_string(i)+", j="+std::to_string(j));
      std::cout<<"elements[0]->Kmat()(i, j) = "<<elements[0]->Kmat()(i, j)<<"\n";
    }
    
  //Element_line2::test();
  */
}

Matrix2d Linear1D::assembleK() {
  int globalNodeCount = 0;
  for(int e=0; e<elements_.size(); ++e) {
    auto eleGlobalNodeIds = elements_(e)->globalNodeIds_;
    for(int locId=0; locId<eleGlobalNodeIds.size(); ++locId) {
      db::pr("line54: "+std::to_string(eleGlobalNodeIds(locId)+1));
      if(eleGlobalNodeIds(locId)+1 > globalNodeCount)
        globalNodeCount = eleGlobalNodeIds(locId)+1;
    }
  }
  Matrix2d assembledK(globalNodeCount, globalNodeCount); // TODO we use a dyn matrix because we want to get rid of the first loop and just have one loop later, but I have to see how I can do that
  db::pr("line44");
  db::pr(std::to_string(globalNodeCount));
  for(int e=0; e<elements_.size(); ++e) {
    const auto& eleGlobalNodeIds = elements_(e)->globalNodeIds_;
    const auto locK = elements_(e)->Kmat();
    for(int i=0; i<locK.nRows(); ++i)
      for(int j=0; j<locK.nCols(); ++j)
        assembledK(eleGlobalNodeIds(i),eleGlobalNodeIds(j)) += locK(i, j);
  }
  
  K_ = assembledK;
  return assembledK;
}

void Linear1D::applyDirichlet(int globalNodeId/*, dofIndex or localDofindex of the node, needed for higher dim*/, double val) {
  // x_d==val means that any occurence of K_ij*x_d, ie when j==d, should be instead subtracted from the rhs. and if i==d, it doesnt matter anyays because we remove that row or practically remove it.
  // ok but actually we want the following approach: remove the rows of dirichlet nodes but then subtract from the rhs the upper right*x_d. ok actually maybe practically the same thing, which makes sense I guess
  int n = K_.nRows();
  // TODO very inefficient but whatever
  
  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      if(j==globalNodeId && K_(i, j)!=0){
        rhs_(i) -= K_(i, j)*val;
      }
    }
  }
  
  
  auto K2 = Matrix2d(n-1, n-1);
  auto rhs2 = Vectord(n-1);
  int i2 = 0;
  for(int i=0; i<n; ++i) {
    if(i!=globalNodeId) {
      int j2 = 0;
      for(int j=0; j<n; ++j) {
        if(j!=globalNodeId){
          K2(i2, j2) = K_(i, j);
          j2++;
        }
      }
      
      rhs2(i2) = rhs_(i);
      
      i2++;
    }
  }
  
  K_ = K2;
  rhs_ = rhs2;
}

Vectord Linear1D::solveSystem_Jacobi() {
  db::pr("Solve start");
  
  const int n = K_.nRows();
  auto L = K_;
  auto U = K_;
  auto D = K_;
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j) {
      if(j>=i)
        L(i, j) = 0;
      if(i>=j)
        U(i, j) = 0;
      if(i!=j)
        D(i, j) = 0;
    }
    
  auto residual = [K = K_, rhs = rhs_](Vectord x){return addVects(mat2TimesVect(K, x), scaleVect(-1.0, rhs));};
  auto residualNorm = [&, residual](Vectord x){auto res = residual(x); return vectDotProduct(res, res);};
  Vectord x_0(n); // Initial guess, zeros here
  Vectord x_i = x_0;
  int iter = 0;
  constexpr int maxiter = 50;
  constexpr double maxResNorm = 0.01;
  
  db::pr("K_ L U D");
  K_.print();
  L.print();
  U.print();
  D.print();
  
  auto invD = D;
  for(int i=0; i<n; ++i)
    invD(i, i) = 1.0/D(i, i);
  db::pr("line 112");
  invD.print();
  while(iter < maxiter/* && residualNorm(x_i) < maxResNorm*/) {
    x_i = mat2TimesVect(invD,
      addVects(rhs_,
        scaleVect(-1.0,
          mat2TimesVect(addMat2s(L, U),
            x_i))));
    iter++;
    std::cout<<"Iter "<<iter<<", resNorm="<<residualNorm(x_i)<<std::endl;
    residual(x_i).print(8);
    x_i.print(8);
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