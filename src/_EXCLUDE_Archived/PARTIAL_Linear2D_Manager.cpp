void Linear2D::example_torus(double x0, double y0, double ri, double ro, int nc, int nr) {
  MyUtils::Timers::ScopedTimer timer("example_torus()");
  
  MyUtils::Timers::StandardTimer timer2("Example meshing");
  timer2.start();
  
  std::vector<int> dirichLeftIds = std::vector<int>();
  std::vector<int> dirichRightIds = std::vector<int>();
  int dofId = 0;
  X_0_ = Vectord();
  FOR(ic, nc) {
    double angleRad = ic*(2.0*MyUtils::Constants::pi/(nc-1));
    double sinAngle = sin(angleRad);
    double cosAngle = cos(angleRad);
    FOR(ir, nr) {
      double r = ri + (ro-ri)*ir/(nr-1);
      double x = x0+cosAngle*r;
      double y = y0+sinAngle*r;
      X_0_.push_back(x);
      X_0_.push_back(y);
      
      if(x<x0-0.95 && ir==nr-1) {
        dirichLeftIds.push_back(dofId);
        dirichLeftIds.push_back(dofId+1);
      }
      else if(x>x0+1.0 && ir==nr-1) {
        dirichRightIds.push_back(dofId);
        dirichRightIds.push_back(dofId+1);
      }
      
      dofId+=2;
    }
  }
    
  auto eleNodes = std::vector<std::vector<int>>();
  FOR(ir, nr-1)
    FOR(ic, nc-2) {
      eleNodes.push_back(std::vector<int>{{ic*nr+ ir, ic*nr+ ir+1, (ic+1)*nr+ ir}});
      eleNodes.push_back(std::vector<int>{{ic*nr+ ir+1, (ic+1)*nr+ ir+1, (ic+1)*nr+ ir}});
    }
  FOR(ir, nr-1) {
    int ic = nc-2;
    eleNodes.push_back(std::vector<int>{{ic*nr+ ir, ic*nr+ ir+1, ir}});
    eleNodes.push_back(std::vector<int>{{ic*nr+ ir+1, ir+1, ir}});
    
    MyUtils::Db::pr("eleNodes.print();");
    StdVectorUtils::print<int>({ic*nr+ ir, ic*nr+ ir+1, (0)*nr+ ir});
    StdVectorUtils::print<int>({ic*nr+ ir+1, (0)*nr+ ir+1, (0)*nr+ ir});
  }
  
  auto YoungPoisson_x = [](double x0, double x1) {
    return Vectord(vector<double>{1, 0.2});
  };
    
  FOR(e, eleNodes.size()) {
    elements_.push_back(new Element::Tri3(
      X_0_,
      eleNodes[e],
      false,
      YoungPoisson_x)
    );
  }
  
  timer2.stop();
  std::cout<<"Built elements\n\n";
  
  /*FOR(i, elements_.size()) {
    std::cout<<"================\nEle"<<i<<"\n";
    elements_(i)->test();
  }*/
  
  
  assembleK();
  //std::cout<<"K_ after assembly:\n";
  //std::cout<<K_.toString(8)<<"n";
  
  rhs_ = Vectord(K_.nRows()); // zero vect for now
  
  auto dirichIds = std::vector<size_t>();
  auto dirichVect = Vectord();
  for(int i = 0; i < dirichLeftIds.size(); i+=2) {
    // left side
    dirichIds.push_back(dirichLeftIds[i]);
    dirichVect.push_back(0.0);
    dirichIds.push_back(dirichLeftIds[i+1]);
    dirichVect.push_back(0.0);
  }
  for(int i = 0; i < dirichRightIds.size(); i+=2) {
    // left side
    dirichIds.push_back(dirichRightIds[i]);
    dirichVect.push_back(-0.1);
    //dirichIds.push_back(dirichRightIds(i+1));
    //dirichVect.push_back(0.5);
  }
  
  MyUtils::Db::pr("dirichRightIds.print();");
  StdVectorUtils::print(dirichRightIds);
  
  MyUtils::Db::pr("dirichLeftIds.print();");
  StdVectorUtils::print(dirichLeftIds);
  
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
  
  applyDirichlet(dirichIds, dirichVect);
      
  
  //std::cout<<"K_ after applyDirichlet():\n";
  //K_.print();
  
  std::cout<<"rhs_ after applyDirichlet():\n";
  rhs_.print();
  
  solutionVect_ = MyUtils::NumMethods::LinSolvers::GaussSeidel(K_, rhs_, 500, 1e-9);
  
  
  //solutionVect_.print(8);
  
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
}
