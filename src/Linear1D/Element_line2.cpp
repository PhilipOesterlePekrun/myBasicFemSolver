#include "Element_line2.hpp"

namespace Element {
  
Element_line2::Element_line2(arrayi<nnodee_> nodes, arrayd<ndof_> X_0)
  : nodes_(nodes), X_0_(X_0) {
    X_t_ = new arrayd<ndof_>(X_0_);
    D_t_xi = new arrayd<ndof_>({0, 0});
    D_t_ = new arrayd<ndof_>({0, 0});
  }

} // namespace Element