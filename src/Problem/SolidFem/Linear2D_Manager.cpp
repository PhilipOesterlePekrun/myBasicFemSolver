#include "Linear2D_Manager.hpp"

#include <SFML/Graphics.hpp> // TODO: delete

namespace MyFem::Problem {
  
Vectord Linear2D::fullSolution() {
  Vectord v = Vectord(globalDofIds_.size());
  
  std::vector<size_t> solutionDofIds(globalDofIds_);
  StdVectorUtils::deleteIndices(solutionDofIds, dirichletDofIds_);
  //std::cout<<"solutionDofIds.print();line11\n";
  //StdVectorUtils::print(solutionDofIds);
  FOR(i, solutionDofIds.size()) {
    v(solutionDofIds[i]) = solutionVect_(i);
  }
  
  FOR(i, dirichletDofIds_.size()) {
    v(dirichletDofIds_[i]) = dirichletVect_(i);
  }
  
  return v;
}

Vectord Linear2D::getX_t() {
  Vectord X_0 = Vectord(getX_0());
  auto d = fullSolution();
  FOR(i, X_0.size()) {
    X_0(i) += d(i);
  }
  return X_0;
}

void Linear2D::assembleK() {
  MyUtils::Timers::ScopedTimer timer("assembleK()");
  
  Matrix2d assembledK(nnode_*ndofn_, nnode_*ndofn_); // TODO: we use a dyn matrix because we want to get rid of the first loop and just have one loop later, but I have to see how I can do that
  
  STATUS("assembleK(): starting assembly...");
  MyUtils::Db::LoadingBar lb(elements_.size(), 100);
  for(int e=0; e<elements_.size(); ++e) {
    const auto& eleGlobalDofIds = elements_[e]->getGlobalDofIds();
    const auto locK = elements_[e]->Kmat();
    for(int i=0; i<locK.nRows(); ++i)
      for(int j=0; j<locK.nCols(); ++j)
        assembledK(eleGlobalDofIds[i],eleGlobalDofIds[j]) += locK(i, j);
    lb();
  }
  
  K_ = assembledK;
  
  //Db::pr("assembledK:");
  //assembledK.print();
}

void Linear2D::assembleM() {
  MyUtils::Timers::ScopedTimer timer("assembleM()");
  
  Matrix2d assembledM(nnode_*ndofn_, nnode_*ndofn_);
  
  STATUS("assembleM(): starting assembly...");
  MyUtils::Db::LoadingBar lb(elements_.size(), 100);
  for(int e=0; e<elements_.size(); ++e) {
    const auto& eleGlobalDofIds = elements_[e]->getGlobalDofIds();
    const auto locK = elements_[e]->Mmat();
    for(int i=0; i<locK.nRows(); ++i)
      for(int j=0; j<locK.nCols(); ++j)
        assembledM(eleGlobalDofIds[i],eleGlobalDofIds[j]) += locK(i, j);
    lb();
  }
  
  M_ = assembledM;
}

void Linear2D::assembleFGravity(double gravityAccel) {
  if(gravityAccel==0)
    rhs_ = Vectord(nnode_ * ndofn_, 0.0);
  else {
    MyUtils::Timers::ScopedTimer timer("assembleFGravity()");

    Vectord assembledF(nnode_ * ndofn_);

    STATUS("assembleFGravity(): starting assembly...");
    MyUtils::Db::LoadingBar lb(elements_.size(), 100);
    for(int e = 0; e < elements_.size(); ++e) {
      const auto& eleGlobalDofIds = elements_[e]->getGlobalDofIds();
      const auto locF = elements_[e]->getFGravity(gravityAccel);

      for(int i = 0; i < locF.size(); ++i)
        assembledF(eleGlobalDofIds[i]) += locF(i);

      lb();
    }

    rhs_ = assembledF;
  }
}

void Linear2D::setNodeDofInfo() {
  MyUtils::Timers::ScopedTimer timer("setNodeDofInfo()");
  
  ndofn_ = elements_[0]->ndofn_;
  
  STATUS("setNodeDofInfo(): getting dof ids");
  
  int globalNodeCount = 0;
  for(int e=0; e<elements_.size(); ++e) {
    auto eleGlobalNodeIds = elements_[e]->globalNodeIds;
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

void Linear2D::assembleFull(double gravityAccel) {
  MyUtils::Timers::ScopedTimer timer("assembleFull()");
  
  setNodeDofInfo();
  assembleK();
  assembleM();
  assembleFGravity(gravityAccel);
  
  if(K_.size() != M_.size()) THROW("assembleFull(): K_.size() != M_.size()");
}

void Linear2D::applyNeumann(const std::vector<size_t>& ids, const Vectord& vals) {
  STATUS("applyNeumann()");
  MyUtils::Timers::ScopedTimer timer("applyNeumann()");
  
  neumannDofIds_ = ids;
  neumannVect_ = vals;
  
  // remove duplicates; if there are duplicates, the later one is ignored (no override)
  for(size_t i = ids.size(); i-- > 0; ) {
    if(StdVectorUtils::find(neumannDofIds_, neumannDofIds_[i]).size()>1) {
      StdVectorUtils::deleteIndices(neumannDofIds_, {i});
      neumannVect_.deleteIndices({i});
    }
  }
  
  FOR(i, neumannDofIds_.size())
    rhs_(neumannDofIds_[i]) += neumannVect_(i);
}
void Linear2D::applyDirichlet(const std::vector<size_t>& ids, const Vectord& vals) {
  STATUS("applyDirichlet()");
  MyUtils::Timers::ScopedTimer timer("applyDirichlet()");
  
  dirichletDofIds_ = ids;
  dirichletVect_ = vals;
  
  // remove duplicates
  for(size_t i = ids.size(); i-- > 0; ) {
    if(StdVectorUtils::find(dirichletDofIds_, dirichletDofIds_[i]).size()>1) {
      StdVectorUtils::deleteIndices(dirichletDofIds_, {i});
      dirichletVect_.deleteIndices({i});
    }
  }
  
  int n = K_.nRows();
  
  vector<size_t> solutionDofIds = globalDofIds_;
  StdVectorUtils::deleteIndices(solutionDofIds, dirichletDofIds_); // or technically deleteIndices(...find(...)) but it will be the same
  
  int n_reduced = solutionDofIds.size();
  Matrix2d K_reduced(n_reduced, n_reduced);
  Vectord rhs_reduced(n_reduced);
  
  FOR(v, n_reduced) {
    int idS = solutionDofIds[v];
    
    double rhs_i = rhs_(idS);
    
    FOR(v2, n_reduced) {
      int id2 = solutionDofIds[v2];
      K_reduced(v, v2) = K_(idS, id2);
    }
    
    FOR(i, dirichletDofIds_.size()) {
      int dId = dirichletDofIds_[i];
      rhs_i -= K_(idS, dId) * dirichletVect_(i);
    }
    rhs_reduced(v) = rhs_i;
  }
  
  rhs_ = rhs_reduced;
  K_ = K_reduced;
}
void Linear2D::applyBCs(const std::vector<size_t>& neumannIds, const Vectord& neumannVals, const std::vector<size_t>& dirichletIds, const Vectord& dirichletVals) {
  applyNeumann(neumannIds, neumannVals);
  applyDirichlet(dirichletIds, dirichletVals);
}

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
  str+="num dirichlet dofs="+to_string(dirichletVect_.size())+"\n";
  str+="num neumann dofs="+to_string(neumannVect_.size())+"\n";
  
  
  return str;
}
void Linear2D::printInfo() {
  std::cout<<infoString()<<"\n";
}

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
  
  
  assembleFull(-9.81);
  //std::cout<<"K_ after assembly:\n";
  //std::cout<<K_.toString(8)<<"n";
  
  std::cout<<"rhs_:\n"<<rhs_.toString();
  
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
    solutionVect_ = MyUtils::NumMethods::LinSolvers::GaussSeidel(K_, rhs_, maxIter, tol, 1);
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

} // namespace MyFem::Problem
