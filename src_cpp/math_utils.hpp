#ifndef MATH_UTILS_TEMPLATE_H
#define MATH_UTILS_TEMPLATE_H

namespace cbasis {

  // -- to be removed
  template<class F>
  class array3 {
  public:
    typedef F Field;

    // ---- member field ----
  private:
    F xs_[3];

  public:
    // ---- constructors ----
    array3(F x, F y, F z) {
      this->xs_[0] = x;
      this->xs_[1] = y;
      this->xs_[2] = z;
    }
    template<class F2>
    array3(const array3<F2>& o) {
      for(int i = 0; i < 3; i++) 
	this->xs_[i] = o[i];
    }

    // ---- accessors ----
    F operator [](int i) const { return xs_[i]; }
    F x() const { return xs_[0]; }
    F y() const { return xs_[1]; }
    F z() const { return xs_[2]; }
    void set_x(F x) { xs_[0] = x; }
    void set_y(F y) { xs_[1] = y; }
    void set_z(F z) { xs_[2] = z; }
  };

  template<class F> F ConjugateIfPossible(F x);

}

#endif
