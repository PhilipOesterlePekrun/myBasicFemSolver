#include "Linear2D_Manager.hpp"

#include <Methods/Polynomial.hpp>

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

// deprecated
void Linear2D::visualize() {// TODO: delete
  Db::warn("visualize() is deprecated");
  
  sf::RenderWindow rw;
  rw.create(sf::VideoMode({1600, 1000}), "tmp");
	rw.setPosition(sf::Vector2i(200,200));
	rw.setFramerateLimit(20);

  while (rw.isOpen()) {
    while (auto event = rw.pollEvent()) {
      if(event->is<sf::Event::Closed>())
        rw.close();
    }

    rw.clear(sf::Color::White);
    
    auto vect = [](double x, double y) {return sf::Vector2f(200+x, 1000-(200+y));};
    
    auto fullSol = fullSolution();
    //auto fullSol = solutionVect_;
    
    sf::ConvexShape quad;
    quad.setFillColor(sf::Color(0,0,0));
    quad.setPointCount(4);
    int scaleBy = 400;
    FOR(i, nnode_) {
      sf::CircleShape pt(10);
      pt.setFillColor(sf::Color(0,240,0));
      pt.setPosition(vect(scaleBy*(X_0_(2*i)), scaleBy*(X_0_(2*i+1))));
      rw.draw(pt);
    }
    FOR(i, nnode_) {
      sf::CircleShape pt(10);
      pt.setFillColor(sf::Color(0,0,0));
      pt.setPosition(vect(scaleBy*(X_0_(2*i)+fullSol(2*i)), scaleBy*(X_0_(2*i+1)+fullSol(2*i+1))));
      rw.draw(pt);
      //quad.setPoint(i, vect(scaleBy*(X_0_(2*i)+fullSol(2*i)), scaleBy*(X_0_(2*i+1)+fullSol(2*i+1))));
    }
    //rw.draw(quad);
    
    rw.display();
  }

}

void Linear2D::runNoInputExample_SingleEle() {
  X_0_ = Vectord({
    0.0, 0.0, // x, y
    2.0, 0.0,
    0.0, 1.0
  });
  
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{0, 1, 2}},
    true
  ));
  
  
  FOR(i, elements_.size()) {
    std::cout<<"================\nEle"<<i<<"\n";
    elements_[i]->test();
  }
  
  
  assembleK();
  std::cout<<K_.toString(8)<<"n";
  
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  
  rhs_ = Vectord(K_.nRows()); // zero vect for now
  
  
  applyDirichlet(2*2+0, 0.0);
  applyDirichlet(2*2+1, 0.0);
  MyUtils::Db::pr("rhs_ after2");
  rhs_.print();
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  applyDirichlet(2*1+0, 0.2);
  
  applyDirichlet(2*0+1, 0.0);
  MyUtils::Db::pr("rhs_ after3");
  rhs_.print();
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  
  K_.print();
  
  solveSystem_Jacobi(500, 0.00001);
  
  MyUtils::Db::pr("line117");
  solutionVect_.print();
  dirichletVect_.print();
  StdVectorUtils::print(dirichletDofIds_);
  
  fullSolution().print();
  //solutionDofIds_.print();
  
  
  
  //visualize();
  /**/
  /*for(int i=0; i<2; ++i)
    for(int j=0; j<2; ++j) {
      MyUtils::Db::pr("line 19, i="+std::to_string(i)+", j="+std::to_string(j));
      std::cout<<"elements[0]->Kmat()(i, j) = "<<elements[0]->Kmat()(i, j)<<"\n";
    }
    
  //Element_line2::test();
  */
}

void Linear2D::runNoInputExample() {
  X_0_ = Vectord({
    0.0, 0.0, // x, y
    2.0, 0.0,
    0.0, 1.0,
    2.0, 1.0
  });
  
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{0, 1, 2}},
    true
  ));
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{1, 3, 2}},
    true
  ));
  
  
  FOR(i, elements_.size()) {
    std::cout<<"================\nEle"<<i<<"\n";
    elements_[i]->test();
  }
  
  
  assembleK();
  std::cout<<K_.toString(8)<<"n";
  
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  
  rhs_ = Vectord(K_.nRows()); // zero vect for now
  
  applyDirichlet(2*3+1, -0.4);
  applyDirichlet(2*3+0, 0.4);
  MyUtils::Db::pr("rhs_ after1");
  rhs_.print();
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  
  
  applyDirichlet(2*2+0, 0.0);
  MyUtils::Db::pr("rhs_ after2");
  rhs_.print();
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  
  applyDirichlet(2*0+1, 0.0);
  applyDirichlet(2*0+0, 0.0);
  MyUtils::Db::pr("rhs_ after3");
  rhs_.print();
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  
  K_.print();
  
  solveSystem_Jacobi(500, 0.00001);
  
  MyUtils::Db::pr("line117");
  solutionVect_.print();
  dirichletVect_.print();
  StdVectorUtils::print(dirichletDofIds_);
  
  fullSolution().print();
  //solutionDofIds_.print();
  
  
  
  visualize();
  /**/
  /*for(int i=0; i<2; ++i)
    for(int j=0; j<2; ++j) {
      MyUtils::Db::pr("line 19, i="+std::to_string(i)+", j="+std::to_string(j));
      std::cout<<"elements[0]->Kmat()(i, j) = "<<elements[0]->Kmat()(i, j)<<"\n";
    }
    
  //Element_line2::test();
  */
}

void Linear2D::example_beam(double lx, double ly, int nx, int ny) {
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
      
      //if(X_0_[2*firstTri[0]]<0.1&&)
      if(!(i>0&&i<nx-2 && j==1)) {
        eleNodes.push_back(firstTri);
        eleNodes.push_back(secondTri);
      }
    }
    
  FOR(e, eleNodes.size()) {
    elements_.push_back(new Element::Tri3(
      X_0_,
      eleNodes[e],
      false)
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
  
  int numDofs = K_.nRows();
  
  rhs_ = Vectord(numDofs); // zeros
  
  // //Neumann BC
  auto neumannIds = std::vector<size_t>();
  auto neumannVect = Vectord();
  
  FOR(j, ny) {
    // from 0 to 1
    double yFrac = (double)j/(ny-1);
    // right side
    neumannIds.push_back(ndofn_*(nx*(j+1)-1)+ 0);
    neumannVect.push_back(4*0.05*yFrac*(1-yFrac)); // second coeff is maximum
    neumannIds.push_back(ndofn_*(nx*(j+1)-1)+ 1);
    neumannVect.push_back(4*0.05*yFrac*(1-yFrac));
  }
  
  FOR(i, nx) {
    if(i==0||i==nx-1) continue;
    // from 0 to 1
    double yFrac = (double)i/(ny-1);
    // upper side
    //neumannIds.push_back(ndofn_*nx*(i+1)-1)+ 1);
    //neumannVect.push_back(4*0.05*yFrac*(1-yFrac));
  }
  
  // K x = F
  FOR(i, neumannIds.size()) // if there are duplicates, later overrides earlier
    rhs_(neumannIds[i]) = neumannVect(i);
    
  std::cout<<"rhs_ after Neumann():\n";
  rhs_.print();
  
  
  // // Dirichlet BC
  auto dirichIds = std::vector<size_t>();
  auto dirichVect = Vectord();
  FOR(j, ny) {
    // left side
    dirichIds.push_back(ndofn_*nx*j+ 0);
    dirichVect.push_back(0.0);
    dirichIds.push_back(ndofn_*nx*j+ 1);
    dirichVect.push_back(0.0);
    
    
    /*dirichIds.push_back(ndofn_*nx*i+ 3);
    dirichVect.push_back(0.0);
    dirichIds.push_back(ndofn_*nx*i+ 4);
    dirichVect.push_back(0.0);*/
    
    
    // right side
    //dirichIds.push_back(ndofn_*(nx*(i+1)-1)+ 0);
    //dirichVect.push_back(1.0);
    //dirichIds.push_back(ndofn_*(nx*(i+1)-1)+ 1);
    //dirichVect.push_back(0.0);
  }
  
  FOR(i, nx) {
    // lower side
    dirichIds.push_back(ndofn_*i+ 0);
    dirichVect.push_back(0.0);
    dirichIds.push_back(ndofn_*i+ 1);
    dirichVect.push_back(0.0);
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
  for(size_t i=dirichIds.size()-1; i>0; --i) {
    if(StdVectorUtils::find(dirichIds, dirichIds[i]).size()>1) {
      std::cout<<"i="<<i<<"\n";
      StdVectorUtils::deleteIndices(dirichIds, {i});
      dirichVect.deleteIndices({i});
    }
  }
  
  applyDirichlet(dirichIds, dirichVect);
      
  
  //std::cout<<"K_ after applyDirichlet():\n";
  //K_.print();
  
  std::cout<<"rhs_ after applyDirichlet():\n";
  rhs_.print();
  
  solveSystem_GaussSeidel(500, 0.0000005);
  
  
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
    
  FOR(e, eleNodes.size()) {
    elements_.push_back(new Element::Tri3(
      X_0_,
      eleNodes[e],
      false)
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
  
  solveSystem_GaussSeidel(500, 1e-9);
  
  
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

void Linear2D::runNoInputExample1() {
  X_0_ = Vectord({
    0.0, 0.0, // x, y
    1.0, 0.0,
    2.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    2.0, 1.0,
  });
  
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{0, 1, 3}},
    true
  ));
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{1, 4, 3}},
    true
  ));
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{1, 2, 4}},
    true
  ));
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{2, 5, 4}},
    true
  ));
  
  
  /*FOR(i, elements_.size()) {
    std::cout<<"================\nEle"<<i<<"\n";
    elements_(i)->test();
  }*/
  
  
  assembleK();
  std::cout<<"K_ after assembly:\n";
  std::cout<<K_.toString(8)<<"n";
  
  rhs_ = Vectord(K_.nRows()); // zero vect for now
  
  applyDirichlet(
    std::vector<size_t>{{2*0+0, 2*2+1, 2*3+0, 2*3+1,
      2*2+1, 5*2+1}},
    Vectord{{0.0, 0.0, 0.0, 0.0,
      -0.2, -0.2}});
  
      
  
  std::cout<<"K_ after applyDirichlet():\n";
  K_.print();
  
  std::cout<<"rhs_ after applyDirichlet():\n";
  rhs_.print();
  
  solveSystem_GaussSeidel(500, 0.0000005);
  
  
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

void Linear2D::runNoInputExample2() {
  X_0_ = Vectord({
    0.0, 0.0, // x, y
    1.0, 0.0,
    2.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    2.0, 1.0,
  });
  
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{0, 1, 3}},
    true
  ));
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{1, 2, 4}},
    true
  ));
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{1, 4, 3}},
    true
  ));
  elements_.push_back(new Element::Tri3(
    X_0_,
    std::vector<int>{{5, 4, 2}},
    true
  ));
  
  
  //elements_(0)->test();
  
  
  assembleK();
  K_.print();
  
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  
  rhs_ = Vectord(K_.nRows()); // zero vect for now
  
  MyUtils::Db::pr("rhs_ after1");
  rhs_.print();
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  applyDirichlet(2*5+1, 0.05);
  
  applyDirichlet(2*3+0, 0.0);
  
  applyDirichlet(2*0+1, 0.0);
  applyDirichlet(2*0+0, 0.0);
  MyUtils::Db::pr("rhs_ after2");
  rhs_.print();
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  
  
  MyUtils::Db::pr("rhs_ after3");
  rhs_.print();
  std::cout<<"globalDofs_\n";
  //solutionDofIds_.print();
  
  K_.print(8);
  
  std::cout<<"line193\n";
  StdVectorUtils::print(dirichletDofIds_);
  dirichletVect_.print();
  //solutionDofIds_.print(); // with removed ids after applying dirichlet
  solutionVect_.print();
  
  solveSystem_Jacobi(500, 0.0001);
  
  fullSolution().print(10);
  
  visualize();
  
  /*for(int i=0; i<2; ++i)
    for(int j=0; j<2; ++j) {
      MyUtils::Db::pr("line 19, i="+std::to_string(i)+", j="+std::to_string(j));
      std::cout<<"elements[0]->Kmat()(i, j) = "<<elements[0]->Kmat()(i, j)<<"\n";
    }
    
  //Element_line2::test();
  */
}

Matrix2d Linear2D::assembleK() {
  MyUtils::Timers::ScopedTimer timer("assembleK()");
  ndofn_ = elements_[0]->ndofn_;
  
  int globalNodeCount = 0;
  for(int e=0; e<elements_.size(); ++e) {
    auto eleGlobalNodeIds = elements_[e]->globalNodeIds_;
    //std::cout<<"line698 globalDofIds_:\n";
    StdVectorUtils::print(globalDofIds_);
    for(int locId=0; locId<eleGlobalNodeIds.size(); ++locId)
      if(StdVectorUtils::find(globalDofIds_, ndofn_*eleGlobalNodeIds[locId]).size()==0)
        if(eleGlobalNodeIds[locId]+1 > globalNodeCount)
          globalNodeCount = eleGlobalNodeIds[locId]+1;
  }
  FOR(i, globalNodeCount) {
        globalDofIds_.push_back(ndofn_*i);
        globalDofIds_.push_back(ndofn_*i + 1);
  }
  
  //std::cout<<"line710 globalDofIds_:\n";
    //StdVectorUtils::print(globalDofIds_);
  
  nnode_ = globalNodeCount;
  Matrix2d assembledK(globalNodeCount*ndofn_, globalNodeCount*ndofn_); // TODO: we use a dyn matrix because we want to get rid of the first loop and just have one loop later, but I have to see how I can do that
  
  std::cout<<"Starting Assembly\n\n";
  for(int e=0; e<elements_.size(); ++e) {
    const auto& eleGlobalDofIds = elements_[e]->getGlobalDofIds();
    const auto locK = elements_[e]->Kmat();
    for(int i=0; i<locK.nRows(); ++i)
      for(int j=0; j<locK.nCols(); ++j)
        assembledK(eleGlobalDofIds[i],eleGlobalDofIds[j]) += locK(i, j);
  }
  
  K_ = assembledK;
  
  //Db::pr("assembledK:");
  //assembledK.print();
  
  return assembledK;
}

// old
void Linear2D::applyDirichlet(int dofId/*, dofIndex or localDofindex of the node, needed for higher dim*/, double val) {
  // x_d==val means that any occurence of K_ij*x_d, ie when j==d, should be instead subtracted from the rhs. and if i==d, it doesnt matter anyays because we remove that row or practically remove it.
  // ok but actually we want the following approach: remove the rows of dirichlet nodes but then subtract from the rhs the upper right*x_d. ok actually maybe practically the same thing, which makes sense I guess
  int n = K_.nRows();
  // TODO very inefficient but whatever
  dirichletDofIds_.push_back(dofId);
  dirichletVect_.push_back(val);
  /*FOR(i, solutionDofIds_.size()) {
    if(solutionDofIds_(i) == globalDofId)
      solutionDofIds_.deleteIndices(std::vector<size_t>{{(size_t)i}});
  }*/
  
  for(int i=0; i<n; ++i) {
    for(int j=0; j<n; ++j) {
      if(j==dofId && K_(i, j)!=0){
        rhs_(i) -= K_(i, j)*val;
      }
    }
  }
  
  // TODOx: fix applyDirichlet
  auto K2 = Matrix2d(n-1, n-1);
  auto rhs2 = Vectord(n-1);
  int i2 = 0;
  for(int i=0; i<n-1; ++i) {
    if(i!=dofId) {
      int j2 = 0;
      for(int j=0; j<n-1; ++j) {
        if(j!=dofId){
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

void Linear2D::applyDirichlet(const std::vector<size_t>& ids, const Vectord& vals) {
  MyUtils::Timers::ScopedTimer timer("applyDirichlet()");
  int n = K_.nRows();
  
  vector<size_t> solutionDofIds = globalDofIds_;
  StdVectorUtils::deleteIndices(solutionDofIds, ids); // or technically deleteIndices(...find(...)) but it will be the same
  
  std::cout<<"applyDirichlet(); print globalDofIds_ and solutionDofIds:\n";
  StdVectorUtils::print(globalDofIds_);
  StdVectorUtils::print(solutionDofIds);
  
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
    
    FOR(i, ids.size()) {
      int dId = ids[i];
      rhs_i -= K_(idS, dId) * vals(i);
    }
    rhs_reduced(v) = rhs_i;
  }
  
  rhs_ = rhs_reduced;
  K_ = K_reduced;
  
  dirichletDofIds_ = ids;
  dirichletVect_ = vals;
}

Vectord Linear2D::solveSystem_Jacobi(int maxiter, double maxRelResNorm) {
  MyUtils::Db::pr("Solve start");
  
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
    
  auto residual = [&](Vectord& x){return vect2dPlusVect2d(mat2dTimesVectd(K_, x), scaleVect2d(-1.0, rhs_));};
  auto residualNorm = [&](Vectord& x){auto res = residual(x); return vect2dDotVect2d(res, res);};
  Vectord x_0(n); // Initial guess, zeros here
  auto relativeResidualNorm = [&](Vectord& x){return residualNorm(x)/residualNorm(x_0);};
  Vectord x_i = x_0;
  int iter = 0;
  
  /*MyUtils::Db::pr("K_ L U D");
  K_.print();
  L.print();
  U.print();
  D.print();*/
  
  auto invD = D;
  for(int i=0; i<n; ++i)
    if(invD(i, i) != 0)
      invD(i, i) = 1.0/D(i, i);
  invD.print();

  while(iter < maxiter && [&](Vectord& x_iIn){if(maxRelResNorm==-1.0) return true; else return relativeResidualNorm(x_iIn) > maxRelResNorm;}(x_i)) {
    x_i = mat2dTimesVectd(invD,
      vect2dPlusVect2d(rhs_,
        scaleVect2d(-1.0,
          mat2dTimesVectd(mat2dPlusMat2d(L, U),
            x_i))));
    iter++;
    std::cout<<"Iter "<<iter<<"\n";
    std::cout<<"relResNorm="<<relativeResidualNorm(x_i)<<"\n";
    std::cout<<"residual = ";
    residual(x_i).print(8);
    std::cout<<"x_i = ";
    x_i.print(8);
    std::cout<<"\n";
  }
  
  solutionVect_ = x_i; // we do both return and set member variable why not
  return x_i;
}

Vectord Linear2D::solveSystem_GaussSeidel(int maxiter, double maxRelResNorm) {
  MyUtils::Timers::ScopedTimer timer("solveSystem_GaussSeidel()");
  MyUtils::Db::pr("Solve start - Gauss Seidel");
  
  const int n = K_.nRows();
  auto L = K_;
  auto U = K_;
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j) {
      if(j>i)
        L(i, j) = 0;
      else// if(i>=j)
        U(i, j) = 0;
    }
    
  auto residual = [&](Vectord& x){return vect2dPlusVect2d(mat2dTimesVectd(K_, x), scaleVect2d(-1.0, rhs_));};
  auto residualNorm = [&](Vectord& x){auto res = residual(x); return vect2dDotVect2d(res, res);};
  Vectord x_0(n); // Initial guess, zeros here
  auto relativeResidualNorm = [&](Vectord& x){return residualNorm(x)/residualNorm(x_0);};
  Vectord x_i = x_0;
  int iter = 0;
  
  /*MyUtils::Db::pr("K_ L U");
  K_.print();
  L.print();
  U.print();*/

  while(iter < maxiter && [&](Vectord& x_iIn){if(maxRelResNorm==-1.0) return true; else return relativeResidualNorm(x_iIn) > maxRelResNorm;}(x_i)) {
    x_i = LinAlg::solveLxb(L,
      vect2dPlusVect2d(rhs_,
        scaleVect2d(-1.0,
          mat2dTimesVectd(U, x_i))));
    iter++;
    std::cout<<"Iter "<<iter<<"\t";
    std::cout<<"relResNorm="<<relativeResidualNorm(x_i)<<"\n";
    //std::cout<<"residual = ";
    //residual(x_i).print(8);
    //std::cout<<"x_i = ";
    //x_i.print(8);
    std::cout<<"\n";
  }
  
  solutionVect_ = x_i; // we do both return and set member variable why not
  return x_i;
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

}