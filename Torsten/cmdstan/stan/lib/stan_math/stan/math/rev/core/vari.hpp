#ifndef STAN_MATH_REV_CORE_VARI_HPP
#define STAN_MATH_REV_CORE_VARI_HPP

#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/prim/meta.hpp>
#include <ostream>
#include <type_traits>

namespace stan {
namespace math {

// forward decleration of vari_value
template <typename T, typename = void>
class vari_value;

// forward declaration of var
template <typename T, typename = void>
class var_value;

/**
 * Abstract base class that all `vari_value` and it's derived classes inherit.
 *
 * The chain() method applies the chain rule. Concrete extensions
 * of this class will represent base variables or the result
 * of operations such as addition or subtraction. These extended
 * classes will store operand variables and propagate derivative
 * information via an implementation of chain().
 */
class vari_base {
 public:
  /**
   * Apply the chain rule to this variable based on the variables
   * on which it depends.
   */
  virtual void chain() = 0;
  virtual void set_zero_adjoint() = 0;

  /**
   * Allocate memory from the underlying memory pool.  This memory is
   * is managed as a whole externally.
   *
   * Warning: Classes should not be allocated with this operator
   * if they have non-trivial destructors.
   *
   * @param nbytes Number of bytes to allocate.
   * @return Pointer to allocated bytes.
   */
  static inline void* operator new(size_t nbytes) noexcept {
    return ChainableStack::instance_->memalloc_.alloc(nbytes);
  }

  /**
   * Delete a pointer from the underlying memory pool.
   *
   * This no-op implementation enables a subclass to throw
   * exceptions in its constructor.  An exception thrown in the
   * constructor of a subclass will result in an error being
   * raised, which is in turn caught and calls delete().
   *
   * See the discussion of "plugging the memory leak" in:
   *   http://www.parashift.com/c++-faq/memory-pools.html
   */
  static inline void operator delete(
      void* /* ignore arg */) noexcept { /* no op */
  }
};

/**
 * The variable implementation for floating point types.
 *
 * This class is complete (not abstract) and may be used for
 * constants.
 *
 * A variable implementation is constructed with a constant
 * value. It also stores the adjoint for storing the partial
 * derivative with respect to the root of the derivative tree.
 *
 */
template <typename T>
class vari_value<T, require_t<std::is_floating_point<T>>> : public vari_base {
 public:
  using value_type = std::decay_t<T>;
  /**
   * The value of this variable.
   */
  const value_type val_;
  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  value_type adj_{0.0};

  /**
   * Construct a variable implementation from a value.  The
   * adjoint is initialized to zero.
   *
   * All constructed variables are added to the stack.  Variables
   * should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.  During
   * derivative propagation, the chain() method of each variable
   * will be called in the reverse order of construction.
   *
   * @tparam S a floating point type.
   * @param x Value of the constructed variable.
   */
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  vari_value(S x) noexcept : val_(x) {  // NOLINT
    ChainableStack::instance_->var_stack_.push_back(this);
  }

  /**
   * Construct a variable implementation from a value.  The
   *  adjoint is initialized to zero and if `stacked` is `false` this vari
   *  will be not be put on the var_stack. Instead it will only be put on
   *  a stack to keep track of whether the adjoint needs to be set to zero.
   *  Variables should be constructed before variables on which they depend
   *  to insure proper partial derivative propagation.  During
   *  derivative propagation, the chain() method of each variable
   *  will be called in the reverse order of construction.
   *
   * @tparam S n floating point type.
   * @param x Value of the constructed variable.
   * @param stacked If false will put this this vari on the nochain stack so
   * that its `chain()` method is not called.
   */
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  vari_value(S x, bool stacked) noexcept : val_(x) {
    if (stacked) {
      ChainableStack::instance_->var_stack_.push_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.push_back(this);
    }
  }

  inline void chain() {}

  /**
   * Initialize the adjoint for this (dependent) variable to 1.
   * This operation is applied to the dependent variable before
   * propagating derivatives, setting the derivative of the
   * result with respect to itself to be 1.
   */
  inline void init_dependent() noexcept { adj_ = 1.0; }

  /**
   * Set the adjoint value of this variable to 0.  This is used to
   * reset adjoints before propagating derivatives again (for
   * example in a Jacobian calculation).
   */
  inline void set_zero_adjoint() noexcept final { adj_ = 0.0; }

  /**
   * Insertion operator for vari. Prints the current value and
   * the adjoint value.
   *
   * @param os [in, out] ostream to modify
   * @param v [in] vari object to print.
   *
   * @return The modified ostream.
   */
  friend std::ostream& operator<<(std::ostream& os, const vari_value<T>* v) {
    return os << v->val_ << ":" << v->adj_;
  }

 private:
  template <typename, typename>
  friend class var_value;
};

// For backwards compatability the default is double
using vari = vari_value<double>;

/**
 * A `vari_view` is used to read from a slice of a `vari_value` with an inner
 * eigen type. It can only accept expressions which do not allocate dynamic
 * memory.
 * @tparam T An eigen expression referencing memory allocated in a `vari_value`.
 */
template <typename T, typename = void>
class vari_view;

/**
 * This struct is follows the CRTP for methods common to `vari_view<>` and
 * `vari_value<Matrix>`.
 * @tparam Derived A `var_value<>` or `vari_view` with an inner type that has
 *  defined methods for subslices of the value and adjoint.
 */
template <typename Derived>
class vari_view_eigen {
 private:
  /**
   * Making the base constructor private while making the derived class a friend
   * help's catch if derived types inherit from another derived types
   * base class. See the fluentcpp article on CRTP for more information.
   */
  vari_view_eigen() = default;
  friend Derived;

  /**
   * Helper function to return a reference to the derived type
   */
  inline Derived& derived() { return static_cast<Derived&>(*this); }
  /**
   * Helper function to return a constant reference to the derived type
   */
  inline const Derived& derived() const {
    return static_cast<const Derived&>(*this);
  }

 public:
  /**
   * A block view of the underlying Eigen matrices.
   * @param start_row Starting row of block.
   * @param start_col Starting columns of block.
   * @param num_rows Number of rows to return.
   * @param num_cols Number of columns to return.
   */
  inline auto block(Eigen::Index start_row, Eigen::Index start_col,
                    Eigen::Index num_rows, Eigen::Index num_cols) const {
    using inner_type = decltype(
        derived().val_.block(start_row, start_col, num_rows, num_cols));
    return vari_view<inner_type>(
        derived().val_.block(start_row, start_col, num_rows, num_cols),
        derived().adj_.block(start_row, start_col, num_rows, num_cols));
  }
  inline auto block(Eigen::Index start_row, Eigen::Index start_col,
                    Eigen::Index num_rows, Eigen::Index num_cols) {
    using inner_type = decltype(
        derived().val_.block(start_row, start_col, num_rows, num_cols));
    return vari_view<inner_type>(
        derived().val_.block(start_row, start_col, num_rows, num_cols),
        derived().adj_.block(start_row, start_col, num_rows, num_cols));
  }

  /**
   * View of the head of Eigen vector types.
   * @param n Number of elements to return from top of vector.
   */
  inline auto head(Eigen::Index n) const {
    using inner_type = decltype(derived().val_.head(n));
    return vari_view<inner_type>(derived().val_.head(n),
                                 derived().adj_.head(n));
  }
  inline auto head(Eigen::Index n) {
    using inner_type = decltype(derived().val_.head(n));
    return vari_view<inner_type>(derived().val_.head(n),
                                 derived().adj_.head(n));
  }

  /**
   * View of the tail of the Eigen vector types.
   * @param n Number of elements to return from bottom of vector.
   */
  inline auto tail(Eigen::Index n) const {
    using inner_type = decltype(derived().val_.tail(n));
    return vari_view<inner_type>(derived().val_.tail(n),
                                 derived().adj_.tail(n));
  }
  inline auto tail(Eigen::Index n) {
    using inner_type = decltype(derived().val_.tail(n));
    return vari_view<inner_type>(derived().val_.tail(n),
                                 derived().adj_.tail(n));
  }

  /**
   * View block of N elements starting at position `i`
   * @param i Starting position of block.
   * @param n Number of elements in block
   */
  inline auto segment(Eigen::Index i, Eigen::Index n) const {
    using inner_type = decltype(derived().val_.segment(i, n));
    return vari_view<inner_type>(derived().val_.segment(i, n),
                                 derived().adj_.segment(i, n));
  }
  inline auto segment(Eigen::Index i, Eigen::Index n) {
    using inner_type = decltype(derived().val_.segment(i, n));
    return vari_view<inner_type>(derived().val_.segment(i, n),
                                 derived().adj_.segment(i, n));
  }

  /**
   * View row of eigen matrices.
   * @param i Row index to slice.
   */
  inline auto row(Eigen::Index i) const {
    using inner_type = decltype(derived().val_.row(i));
    return vari_view<inner_type>(derived().val_.row(i), derived().adj_.row(i));
  }
  inline auto row(Eigen::Index i) {
    using inner_type = decltype(derived().val_.row(i));
    return vari_view<inner_type>(derived().val_.row(i), derived().adj_.row(i));
  }

  /**
   * View column of eigen matrices
   * @param i Column index to slice
   */
  inline auto col(Eigen::Index i) const {
    using inner_type = decltype(derived().val_.col(i));
    return vari_view<inner_type>(derived().val_.col(i), derived().adj_.col(i));
  }
  inline auto col(Eigen::Index i) {
    using inner_type = decltype(derived().val_.col(i));
    return vari_view<inner_type>(derived().val_.col(i), derived().adj_.col(i));
  }

  /**
   * Get coefficient of eigen matrices
   * @param i Row index
   * @param j Column index
   */
  inline auto coeff(Eigen::Index i, Eigen::Index j) const {
    return vari_value<double>(derived().val_.coeffRef(i, j),
                              derived().adj_.coeffRef(i, j));
  }
  inline auto coeff(Eigen::Index i, Eigen::Index j) {
    return vari_value<double>(derived().val_.coeffRef(i, j),
                              derived().adj_.coeffRef(i, j));
  }

  /**
   * Get coefficient of eigen matrices
   * @param i Column index to slice
   */
  inline auto coeff(Eigen::Index i) const {
    return vari_value<double>(derived().val_.coeffRef(i),
                              derived().adj_.coeffRef(i));
  }
  inline auto coeff(Eigen::Index i) {
    return vari_value<double>(derived().val_.coeffRef(i),
                              derived().adj_.coeffRef(i));
  }

  /**
   * Get coefficient of eigen matrices
   * @param i Column index to slice
   */
  inline auto operator()(Eigen::Index i) const { return this->coeff(i); }
  inline auto operator()(Eigen::Index i) { return this->coeff(i); }

  /**
   * Get coefficient of eigen matrices
   * @param i Row index
   * @param j Column index
   */
  inline auto operator()(Eigen::Index i, Eigen::Index j) const {
    return this->coeff(i, j);
  }
  inline auto operator()(Eigen::Index i, Eigen::Index j) {
    return this->coeff(i, j);
  }

  /**
   * Return an expression that operates on the rows of the matrix `vari`
   */
  inline auto rowwise_reverse() const {
    using inner_type = decltype(derived().val_.rowwise().reverse());
    return vari_view<inner_type>(derived().val_.rowwise().reverse(),
                                 derived().adj_.rowwise().reverse());
  }
  inline auto rowwise_reverse() {
    using inner_type = decltype(derived().val_.rowwise().reverse());
    return vari_view<inner_type>(derived().val_.rowwise().reverse(),
                                 derived().adj_.rowwise().reverse());
  }

  /**
   * Return an expression that operates on the columns of the matrix `vari`
   */
  inline auto colwise_reverse() const {
    using inner_type = decltype(derived().val_.colwise().reverse());
    return vari_view<inner_type>(derived().val_.colwise().reverse(),
                                 derived().adj_.colwise().reverse());
  }
  inline auto colwise_reverse() {
    using inner_type = decltype(derived().val_.colwise().reverse());
    return vari_view<inner_type>(derived().val_.colwise().reverse(),
                                 derived().adj_.colwise().reverse());
  }

  /**
   * Return an expression an expression to reverse the order of the coefficients
   * inside of a `vari` matrix
   */
  inline auto reverse() const {
    using inner_type = decltype(derived().val_.reverse());
    return vari_view<inner_type>(derived().val_.reverse(),
                                 derived().adj_.reverse());
  }
  inline auto reverse() {
    using inner_type = decltype(derived().val_.reverse());
    return vari_view<inner_type>(derived().val_.reverse(),
                                 derived().adj_.reverse());
  }

  /**
   * Return the number of rows for this class's `val_` member
   */
  inline Eigen::Index rows() const { return derived().val_.rows(); }
  /**
   * Return the number of columns for this class's `val_` member
   */
  inline Eigen::Index cols() const { return derived().val_.cols(); }
  /**
   * Return the size of this class's `val_` member
   */
  inline Eigen::Index size() const { return derived().val_.size(); }
};

template <typename T>
class vari_view<T, require_not_plain_type_t<T>> final
    : public vari_base,
      public vari_view_eigen<vari_view<T, require_not_plain_type_t<T>>> {
 public:
  using PlainObject = plain_type_t<T>;
  using value_type = std::decay_t<T>;  // The underlying type for this class
  /**
   * Number of rows known at compile time
   */
  static constexpr int RowsAtCompileTime = PlainObject::RowsAtCompileTime;
  /**
   * Number of columns known at compile time
   */
  static constexpr int ColsAtCompileTime = PlainObject::ColsAtCompileTime;

  T val_;
  T adj_;
  template <typename S, typename K,
            require_convertible_t<S&, value_type>* = nullptr,
            require_convertible_t<K&, value_type>* = nullptr>
  vari_view(const S& val, const K& adj) noexcept : val_(val), adj_(adj) {}

  void set_zero_adjoint() {}
  void chain() {}
};

/**
 * The variable implementation for Eigen dense matrix types.
 *
 * This class is complete (not abstract) and may be used for
 * constants.
 *
 * A variable implementation is constructed with a constant
 * value. It also stores the adjoint for storing the partial
 * derivative with respect to the root of the derivative tree.
 *
 */
template <typename T>
class vari_value<T, require_all_t<is_plain_type<T>, is_eigen_dense_base<T>>>
    : public vari_base,
      public vari_view_eigen<vari_value<
          T, require_all_t<is_plain_type<T>, is_eigen_dense_base<T>>>> {
 public:
  /**
   * `PlainObject` represents a user constructible type such as Matrix or Array
   */
  using PlainObject = plain_type_t<T>;
  using value_type = PlainObject;  // The underlying type for this class
  using eigen_scalar = value_type_t<PlainObject>;  // A floating point type
  using vari_type = vari_value<T>;
  /**
   * Number of rows known at compile time
   */
  static constexpr int RowsAtCompileTime = PlainObject::RowsAtCompileTime;
  /**
   * Number of columns known at compile time
   */
  static constexpr int ColsAtCompileTime = PlainObject::ColsAtCompileTime;

  /**
   * The value of this variable.
   */
  arena_matrix<PlainObject> val_;

  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  arena_matrix<PlainObject> adj_;

  /**
   * Construct a dense Eigen variable implementation from a value. The
   * adjoint is initialized to zero.
   *
   * All constructed variables are added to the stack. Variables
   * should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.  During
   * derivative propagation, the chain() method of each variable
   * will be called in the reverse order of construction.
   *
   * @tparam S A dense Eigen type that is convertible to `value_type`
   * @param x Value of the constructed variable.
   */
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  explicit vari_value(const S& x) : val_(x), adj_(x.rows(), x.cols()) {
    adj_.setZero();
    ChainableStack::instance_->var_stack_.push_back(this);
  }

  /**
   * Construct a dense Eigen variable implementation from a value. The
   *  adjoint is initialized to zero and if `stacked` is `false` this vari
   *  will be not be put on the var_stack. Instead it will only be put on
   *  a stack to keep track of whether the adjoint needs to be set to zero.
   *  Variables should be constructed before variables on which they depend
   *  to insure proper partial derivative propagation.  During
   *  derivative propagation, the chain() method of each variable
   *  will be called in the reverse order of construction.
   *
   * @tparam S A dense Eigen type that is convertible to `value_type`
   * @param x Value of the constructed variable.
   * @param stacked If false will put this this vari on the nochain stack so
   * that its `chain()` method is not called.
   */
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  vari_value(const S& x, bool stacked) : val_(x), adj_(x.rows(), x.cols()) {
    adj_.setZero();
    if (stacked) {
      ChainableStack::instance_->var_stack_.push_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.push_back(this);
    }
  }

  virtual void chain() {}
  /**
   * Initialize the adjoint for this (dependent) variable to 1.
   * This operation is applied to the dependent variable before
   * propagating derivatives, setting the derivative of the
   * result with respect to itself to be 1.
   */
  inline void init_dependent() { adj_.setOnes(); }

  /**
   * Set the adjoint value of this variable to 0.  This is used to
   * reset adjoints before propagating derivatives again (for
   * example in a Jacobian calculation).
   */
  inline void set_zero_adjoint() final { adj_.setZero(); }

  /**
   * Insertion operator for vari. Prints the current value and
   * the adjoint value.
   *
   * @param os [in, out] ostream to modify
   * @param v [in] vari object to print.
   *
   * @return The modified ostream.
   */
  friend std::ostream& operator<<(std::ostream& os, const vari_value<T>* v) {
    return os << "val: \n" << v->val_ << " \nadj: \n" << v->adj_;
  }

 private:
  template <typename, typename>
  friend class var_value;
};

/**
 * The variable implementation for Eigen sparse matrix types.
 *
 * This class is complete (not abstract) and may be used for
 * constants.
 *
 * A variable implementation is constructed with a constant
 * value. It also stores the adjoint for storing the partial
 * derivative with respect to the root of the derivative tree.
 *
 */
template <typename T>
class vari_value<T, require_eigen_sparse_base_t<T>> : public vari_base,
                                                      chainable_alloc {
 public:
  using PlainObject = plain_type_t<T>;  // Base type of Eigen class
  using value_type = PlainObject;       // vari's adj_ and val_ member type
  /**
   * Rows at compile time
   */
  static constexpr int RowsAtCompileTime = T::RowsAtCompileTime;
  /**
   * Columns at compile time
   */
  static constexpr int ColsAtCompileTime = T::ColsAtCompileTime;

  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  PlainObject adj_;

  /**
   * The value of this variable.
   */
  const PlainObject val_;

  /**
   * Construct a variable implementation from a value. The
   * adjoint is initialized to zero.
   *
   * All constructed variables are added to the stack. For a sparse eigen matrix
   * this includes the nozero values as well the inner and outer indices.
   * Variables should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.  During
   * derivative propagation, the chain() method of each variable
   * will be called in the reverse order of construction.
   *
   * @tparam S A sparse Eigen type that is convertible to `value_type`
   * @param x Value of the constructed variable.
   */
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  explicit vari_value(S&& x)
      : adj_(x), val_(std::forward<S>(x)), chainable_alloc() {
    this->set_zero_adjoint();
    ChainableStack::instance_->var_stack_.push_back(this);
  }
  /**
   * Construct an sparse Eigen variable implementation from a value. The
   *  adjoint is initialized to zero and if `stacked` is `false` this vari
   *  will be not be put on the var_stack. Instead it will only be put on
   *  a stack to keep track of whether the adjoint needs to be set to zero.
   *
   * All constructed variables are added to a stack.  Variables
   *  should be constructed before variables on which they depend
   *  to insure proper partial derivative propagation.  During
   *  derivative propagation, the chain() method of each variable
   *  will be called in the reverse order of construction.
   *
   * @tparam S A sparse Eigen type that is convertible to `value_type`
   * @param x Value of the constructed variable.
   * @param stacked If false will put this this vari on the nochain stack so
   * that its `chain()` method is not called.
   */
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  vari_value(S&& x, bool stacked)
      : adj_(x), val_(std::forward<S>(x)), chainable_alloc() {
    this->set_zero_adjoint();
    if (stacked) {
      ChainableStack::instance_->var_stack_.push_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.push_back(this);
    }
  }

  /**
   * Return the number of rows for this class's `val_` member
   */
  Eigen::Index rows() const { return val_.rows(); }
  /**
   * Return the number of columns for this class's `val_` member
   */
  Eigen::Index cols() const { return val_.cols(); }
  /**
   * Return the size of this class's `val_` member
   */
  Eigen::Index size() const { return val_.size(); }

  void chain() {}
  /**
   * Initialize the adjoint for this (dependent) variable to 1.
   * This operation is applied to the dependent variable before
   * propagating derivatives, setting the derivative of the
   * result with respect to itself to be 1.
   */
  inline void init_dependent() {
    for (int k = 0; k < adj_.outerSize(); ++k) {
      for (typename PlainObject::InnerIterator it(adj_, k); it; ++it) {
        it.valueRef() = 1.0;
      }
    }
  }

  /**
   * Set the adjoint value of this variable to 0.  This is used to
   * reset adjoints before propagating derivatives again (for
   * example in a Jacobian calculation).
   */
  inline void set_zero_adjoint() noexcept final {
    for (int k = 0; k < adj_.outerSize(); ++k) {
      for (typename PlainObject::InnerIterator it(adj_, k); it; ++it) {
        it.valueRef() = 0.0;
      }
    }
  }

  /**
   * Insertion operator for vari. Prints the current value and
   * the adjoint value.
   *
   * @param os [in, out] ostream to modify
   * @param v [in] vari object to print.
   *
   * @return The modified ostream.
   */
  friend std::ostream& operator<<(std::ostream& os, const vari_value<T>* v) {
    return os << "val: \n" << v->val_ << " \nadj: \n" << v->adj_;
  }

 private:
  template <typename, typename>
  friend class var_value;
};

}  // namespace math
}  // namespace stan
#endif
