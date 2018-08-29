/*
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date:   28.12.2017
 */

#include "bits/MatrixArithmetic.hpp"

namespace anpi
{

  // -------------------------------------------
  // Implementation of Matrix::_Matrix_impl
  // -------------------------------------------

  template<typename T,class Alloc>
  Matrix<T,Alloc>::_Matrix_impl::_Matrix_impl()
    : allocator_type(), _data(), _rows(), _cols(), _dcols() { }

  template<typename T,class Alloc>
  Matrix<T,Alloc>::_Matrix_impl::
  _Matrix_impl(allocator_type const& _a) noexcept
    : allocator_type(_a), _data(), _rows(), _cols(), _dcols() { }
      
  template<typename T,class Alloc>
  Matrix<T,Alloc>::_Matrix_impl::
  _Matrix_impl(allocator_type&& _a) noexcept
    : allocator_type(std::move(_a)),
      _data(), _rows(), _cols(), _dcols() { }
  
  template<typename T,class Alloc>
  void Matrix<T,Alloc>::_Matrix_impl::
  _swap_data(_Matrix_impl& _x) noexcept {
    std::swap(_data,  _x._data);
    std::swap(_rows,  _x._rows);
    std::swap(_cols,  _x._cols);
    std::swap(_dcols, _x._dcols);
  }
     
  // ------------------------
  // Implementation of Matrix
  // ------------------------

  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix() : _impl() {}

  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(const allocator_type& _a) noexcept
    : _impl(_a) { }

  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(const size_t _rows,
                          const size_t _cols,
                          const value_type _initVal)
    : Matrix(_rows,_cols,DoNotInitialize) {
    fill(_initVal);
  }

  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(const size_t _rows,
                          const size_t _cols,
                          const value_type _initVal,
                          const allocator_type& _a)
    : Matrix(_rows,_cols,DoNotInitialize,_a) {
    fill(_initVal);
  }

  
  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(const size_t _rows,
                          const size_t _cols,
                          const InitializationType)
    : _impl() {
    _create_storage(_rows,_cols);
  }

  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(const size_t _rows,
                          const size_t _cols,
                          const InitializationType,
                          const allocator_type& _a)
    : _impl(_a) {
    _create_storage(_rows,_cols);
  }
  
  /*
   * Construct a matrix rows x cols and initialize all
   * elements with the memory content at the given pointer
   */
  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(const size_t _rows,
                          const size_t _cols,
                          const const_pointer _initMem)
    : Matrix(_rows,_cols,DoNotInitialize) {
    fill(_initMem);
  }

  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(const size_t _rows,
                          const size_t _cols,
                          const const_pointer _initMem,
                          const allocator_type& _a)
    : Matrix(_rows,_cols,DoNotInitialize,_a) {
    fill(_initMem);
  }
  
  template<typename T,class Alloc>
  Matrix<T,Alloc>::
  Matrix(std::initializer_list< std::initializer_list<value_type> > _lst)
    : Matrix(_lst.size(),
             (_lst.size()>0) ? _lst.begin()->size() : 0,
             DoNotInitialize) {
    fill(_lst);
  }

  template<typename T,class Alloc>
  Matrix<T,Alloc>::
  Matrix(std::initializer_list< std::initializer_list<value_type> > _lst,
         const allocator_type& _a)
    : Matrix(_lst.size(),
             (_lst.size()>0) ? _lst.begin()->size() : 0,
             DoNotInitialize, _a) {
    fill(_lst);
  }
  

  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(const Matrix<T,Alloc>& _other)
    : Matrix(_other.rows(),_other.cols(),DoNotInitialize) {
    
    fill(_other.data());
  }

  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(const Matrix<T,Alloc>& _other,
                          const allocator_type& _a)
    : Matrix(_other.rows(),_other.cols(),DoNotInitialize,_a) {
    
    fill(_other.data());
  }
  
  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(Matrix<T,Alloc>&& _other)
    : _impl(std::move(_other._get_allocator())) {
    this->_impl._swap_data(_other._impl);
  }

  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(Matrix<T,Alloc>&& _other,
                          const allocator_type& _a)
    : _impl(_a) {

    if (_other._get_allocator() == _a) {
      this->_impl._swap_data(_other._impl);
    } else {
      _create_storage(_other._impl._rows,_other._impl._cols);
    }
  }
  
  template<typename T,class Alloc>
  Matrix<T,Alloc>::Matrix(allocator_type&& _a) noexcept
    : _impl(std::move(_a)) { }

  
  template<typename T,class Alloc>
  Matrix<T,Alloc>::~Matrix() noexcept {
    _deallocate();
  }
  
  template<typename T,class Alloc>
  Matrix<T,Alloc>& Matrix<T,Alloc>::operator=(const Matrix<T,Alloc>& other) {
    allocate(other._impl._rows,other._impl._cols);
    fill(other.data());
    return *this;
  }

  template<typename T,class Alloc>
  Matrix<T,Alloc>& Matrix<T,Alloc>::operator=(Matrix<T,Alloc>&& other) {
    if (this->data() != other.data() ) { // alias detection first
      this->_impl._swap_data(other._impl);
    }
    other.clear();
    return *this;
  }
  
  template<typename T,class Alloc>
  bool Matrix<T,Alloc>::operator==(const Matrix<T,Alloc>& other) const {
    if (&other==this) return true; // alias detection

    // same size of matrices?
    if ((other.rows() != this->rows()) ||
        (other.cols() != this->cols())) return false;

    // check the content with pointers
    if (this->_impl._dcols == this->_impl._cols)
      return (memcmp(this->_impl._data,
                     other._impl._data,
                     this->_impl.tentries()*sizeof(T))==0);

    // we have to compare row by row, becase the padding may differ
    for (size_t i=0;i<this->_impl._rows;++i) {
      if (memcmp(this->operator[](i),
                 other[i],
                 this->_impl._cols*sizeof(T))!=0) {
        return false;
      }
    }
    return true;
  }

  template<typename T,class Alloc>
  bool Matrix<T,Alloc>::operator!=(const Matrix<T,Alloc>& other) const {
    if (&other==this) return false; // alias detection
    
    return !operator==(other);
  }

  template<typename T,class Alloc>
  void Matrix<T,Alloc>::swap(Matrix<T,Alloc>& other) {
    this->_impl._swap_data(other._impl);
  }
    
  template<typename T,class Alloc>
  void Matrix<T,Alloc>::allocate(const size_t r,
                                 const size_t c) {
    // only reserve iff the desired size is different to the current one
    if ( (r!=rows()) || (c!=cols()) ) {
      _deallocate();
      _create_storage(r,c);
    }
  }

  template<typename T,class Alloc>
  void Matrix<T,Alloc>::clear() {
    _deallocate();
  }

  
  template<typename T,class Alloc>
  void Matrix<T,Alloc>::_create_storage(size_t __rows,size_t __cols) {

    // ensure that the type T fits into the alignment
    static_assert( ( (_Matrix_impl::alignment <= sizeof(T) )
                     ? (sizeof(T) % _Matrix_impl::alignment  == 0)
                     : (_Matrix_impl::alignment % sizeof(T) == 0) ),
                   "Type of elements must fit an integral number of "
                   "times into an aligned block" );

    // report an error if this has not previously been initialized
    assert(this->_impl._data == nullptr);

    size_t n,dcols;

    if (_Matrix_impl::rowAlign) {
      // how many aligned "blocks" are required to hold __cols
      const size_t blocks = (__cols*sizeof(T) + (_Matrix_impl::alignment-1) ) /
                            _Matrix_impl::alignment;
      // dominant columns are determined by the # blocks per row
      dcols               = blocks*_Matrix_impl::alignment/sizeof(T);
      // total number of entries already padded
      n                   = __rows*dcols;
      
    } else { // do not align the rows, just the complete memory block

      // total number of blocks
      const size_t blocks
        = (__cols*__rows*sizeof(T)+(_Matrix_impl::alignment-1) ) /
          _Matrix_impl::alignment;
      // dominant columns is the same as columns in this case
      dcols = __cols;
      // the total number of entries of type T to be allocated 
      n     = blocks*_Matrix_impl::alignment/sizeof(T);
    } 
          
    // Call the allocator to reserve the required memory
    this->_impl._data
      = (n != 0)
      ? std::allocator_traits<allocator_type>::allocate(_impl, n) 
      : pointer();
    
    // Initialize the rest of the attributes
    this->_impl._rows = __rows;
    this->_impl._cols = __cols;
    this->_impl._dcols = dcols;
  }

  template<typename T,class Alloc>
  void Matrix<T,Alloc>::_deallocate() {
    if (this->_impl._data) {
      std::allocator_traits<allocator_type>::deallocate(this->_impl,
                                                        this->_impl._data,
                                                        this->_impl.tentries());
    }
    
    this->_impl._data  = 0;
    this->_impl._rows  = 0;
    this->_impl._cols  = 0;
    this->_impl._dcols = 0;
  }

  template<typename T,class Alloc>
  typename Matrix<T,Alloc>::allocator_type&
  Matrix<T,Alloc>::_get_allocator() noexcept {
    return *static_cast<allocator_type*>(&this->_impl);
  }
    
  template<typename T,class Alloc>
  const typename Matrix<T,Alloc>::allocator_type&
  Matrix<T,Alloc>::_get_allocator() const noexcept {
    return *static_cast<const allocator_type*>(&this->_impl);
  }
  
  template<typename T,class Alloc>
  void Matrix<T,Alloc>::fill(const T val) {
    T* end = this->_impl._data + ( this->_impl._rows * this->_impl._dcols );
    for (T* ptr = this->_impl._data;ptr!=end;++ptr) {
      *ptr = val;
    }

    // ACCELERATE ME!
  }

  template<typename T,class Alloc>
  void Matrix<T,Alloc>::fill(const T* mem) {
    std::memcpy(this->_impl._data,mem,sizeof(T)*this->_impl.tentries());
  }

  template<typename T,class Alloc>
  template<class OAlloc>
  void Matrix<T,Alloc>::fill(const Matrix<T,OAlloc>& _other) {

    // we can only copy this number of rows
    const size_t r=std::min(_other.rows(),this->rows());

    // we can only copy this number of columns
    const size_t c=std::min(_other.cols(),this->cols());

    // copy each row separately, ignoring the differences of sizes
    for (size_t i=0u;i<r;++i) {
      std::memcpy(this->operator[](i),_other[i],sizeof(value_type)*c);
    }
  }

  

  template<typename T,class Alloc>
  void Matrix<T,Alloc>::
  fill(const std::initializer_list< std::initializer_list<value_type> >& lst) {

    const size_t r = lst.size();
    const size_t c = (r>0u) ? lst.begin()->size() : 0u;
    
    assert(r==rows() && "Check number of rows");
    assert(c==cols() && "Check number of cols");

    pointer rowPtr = this->_impl._data;
    for (const auto& r : lst) {
      pointer ptr=rowPtr;
      
      for (const auto& c : r) {
        *ptr++ = c;
      }
      rowPtr += this->_impl._dcols;
    }
  }

  template<typename T,class Alloc>
  std::vector< typename Matrix<T,Alloc>::value_type >
  Matrix<T,Alloc>::column(const size_t col) const {
    
    std::vector< typename Matrix<T,Alloc>::value_type > vct(this->_impl._rows);

    const_pointer ptr = this->_impl._data + col;    
    for (auto it=vct.begin(); it!=vct.end(); ++it, ptr+=this->_impl._dcols) {
      *it = *ptr;
    }

    return vct;
  }
  
  template<typename T,class Alloc>
  Matrix<T,Alloc>& Matrix<T,Alloc>::operator+=(const Matrix<T,Alloc>& other) {

    ::anpi::aimpl::add(*this,other);
    
    return *this;
  }

  template<typename T,class Alloc>
  Matrix<T,Alloc>& Matrix<T,Alloc>::operator-=(const Matrix<T,Alloc>& other) {

    ::anpi::aimpl::subtract(*this,other);
      
    return *this;
  }

  template<typename T,class Alloc>
  Matrix<T,Alloc> operator+(const Matrix<T,Alloc>& a,
                            const Matrix<T,Alloc>& b) {

    assert( (a.rows()==b.rows()) && (a.cols()==b.cols()) );    

    Matrix<T,Alloc> c(a.rows(),a.cols(),anpi::DoNotInitialize);
    ::anpi::aimpl::add(a,b,c);
    return c;
  }

  template<typename T,class Alloc>
  Matrix<T,Alloc> operator-(const Matrix<T,Alloc>& a,
                            const Matrix<T,Alloc>& b) {
    
    assert( (a.rows()==b.rows()) && (a.cols()==b.cols()) );
    
    Matrix<T,Alloc> c(a.rows(),a.cols(),anpi::DoNotInitialize);
    ::anpi::aimpl::subtract(a,b,c);
    return c;
  }
  
} // namespace ANPI
