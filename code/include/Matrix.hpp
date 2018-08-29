/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 28.12.2017
 */

#ifndef ANPI_MATRIX_HPP
#define ANPI_MATRIX_HPP

#include <cstddef>
#include <cstring>
#include <cassert>
#include <memory>

#include <initializer_list>

#include <AnpiConfig.hpp>
#include <Allocator.hpp>

#include <typeinfo>

namespace anpi
{
  /**
   * Type to make explicit a desired unitialized matrix
   */
  enum InitializationType {
    DoNotInitialize
  };
  
  /**
   * Row-major matrix class.
   *
   * The allocator is used to reserve the memory.  If the allocator
   * has a static const attribute named Alignment, then each row is
   * padded to ensure that also each row is properly aligned.  In such
   * case, the dcols() method will return the total number of entries
   * (including padding), while cols() returns the effective number of
   * entries.  The value of cols() is always less or equal to dcols()
   *
   * To just align the memory buffer (without any alignment), you can
   * use anpi::aligned_allocator and for forcing the alignment of each
   * row you can use anpi::aligned_row_allocator, both defined in
   * <Allocator.hpp>.
   */
  template<typename T,class Alloc=anpi::aligned_row_allocator<T> >
  class Matrix {
  public:   
    /**
     * @name Standard types
     */
    //@{
    typedef typename std::allocator_traits<Alloc>::template
      rebind_alloc<T> allocator_type;

    typedef typename std::allocator_traits<allocator_type>::value_type
      value_type;

    typedef typename std::allocator_traits<allocator_type>::pointer
      pointer;
      
    typedef typename std::pointer_traits<pointer>::template
      rebind<const value_type> const_pointer;
    //@}
   
  protected:
    /**
     * Data is held in an instance of this.
     *
     * This class has limited functionality but ensures
     * exception safety.
     */
    struct _Matrix_impl : public allocator_type {
      /// All matrix data
      pointer _data;
      
      /// Number of rows
      size_t _rows;
      
      /// Effective number of columns
      size_t _cols;
      
      /// Dominant (real) number of columns
      size_t _dcols;

      /// Return the total number of entries (i.e. real buffer size)
      inline size_t tentries() const {return _rows*_dcols;}

      /// Alignment in use for the rows
      static constexpr size_t alignment =
        extract_alignment<allocator_type>::value;

      /// Allocator indicates we must align each row
      static constexpr bool rowAlign =
        extract_alignment<allocator_type>::row_aligned;
      
      /**
       * @name Default constructor
       */
      //@{
      _Matrix_impl();
      _Matrix_impl(allocator_type const& _a) noexcept;
      _Matrix_impl(allocator_type&& _a) noexcept;
      //@}
      
      void _swap_data(_Matrix_impl& _x) noexcept;
    };

    /// The instance doing all the allocation work
    _Matrix_impl _impl;


  public:

    /**
     * @name Constructors
     */
    //@{
    Matrix();
    Matrix(const allocator_type& _a) noexcept;
    explicit Matrix(const size_t _rows,
                    const size_t _cols,
                    const InitializationType _it);
    explicit Matrix(const size_t rows,
                    const size_t cols,
                    const InitializationType _it,
                    const allocator_type& _a);
    explicit Matrix(const size_t _rows,
                    const size_t _cols,
                    const value_type _initVal=value_type());
    explicit Matrix(const size_t _rows,
                    const size_t _cols,
                    const value_type _initVal,
                    const allocator_type& _a);
    explicit Matrix(const size_t _rows,
                    const size_t _cols,
                    const const_pointer _initMem);
    explicit Matrix(const size_t _rows,
                    const size_t _cols,
                    const const_pointer _initMem,
                    const allocator_type& _a);
    Matrix(const Matrix<T,Alloc>& _other);
    Matrix(const Matrix<T,Alloc>& _other,const allocator_type& _a);
    Matrix(allocator_type&& _a) noexcept;
    Matrix(Matrix<T,Alloc>&& _other);
    Matrix(Matrix<T,Alloc>&& _other,const allocator_type& _a);
    ~Matrix() noexcept;
    
    /**
     * Constructs a matrix from a std::initializer_list
     *
     * This allows to construct and initialize a 2x3 matrix in this way:
     *
     * \code
     * anpi::Matrix<int> m={ {1,2,3}, {4,5,6} };
     * \endcode
     */
    Matrix(std::initializer_list< std::initializer_list<value_type> > _lst);
    Matrix(std::initializer_list< std::initializer_list<value_type> > _lst,
           const allocator_type& _a);
    
    //@}

    
    /**
     * Deep copy another matrix
     */
    Matrix<T,Alloc>& operator=(const Matrix<T,Alloc>& other);

    /**
     * Move assignment operator
     */
    Matrix<T,Alloc>& operator=(Matrix<T,Alloc>&& other);

    /**
     * Compare two matrices for equality
     *
     * This is slow, as all componentes are elementwise compared
     */
    bool operator==(const Matrix<T,Alloc>& other) const;

    /**
     * Compare two matrices for equality
     *
     * This is slow, as all componentes are elementwise compared
     */
    bool operator!=(const Matrix<T,Alloc>& other) const;
    
    /// Return pointer to a given row
    inline T* operator[](const size_t row) {
      return this->_impl._data + row * this->_impl._dcols;
    }

    /// Return read-only pointer to a given row
    const T* operator[](const size_t row) const {
      return this->_impl._data + row * this->_impl._dcols;
    }

    /// Return reference to the element at the r row and c column
    T& operator()(const size_t row,const size_t col) {
      return *(this->_impl._data +
               (row*this->_impl._dcols + col));
    }

    /// Return const reference to the element at the r row and c column
    const T& operator()(const size_t row,const size_t col) const {
      return *(this->_impl._data +
               (row*this->_impl._dcols + col));
    }

    /**
     * Swap the contents of the other matrix with this one
     */
    void swap(Matrix<T,Alloc>& other);
    
    /**
     * Allocate memory for the given number of rows and cols
     */
    void allocate(const size_t row,const size_t col);

    /**
     * Reset this matrix to a default constructed empty state
     */
    void clear();
    
    /**
     * Fill all elements of the matrix with the given value
     */
    void fill(const T val);

    /**
     * Fill all elements of the matrix with the given memory block
     *
     * The user must ensure that the given memory block has enough elements
     * and with the same layout used in this matrix, including the 
     * row padding.
     */
    void fill(const T* mem);

    /**
     * Fill this matrix with the content of the other matrix.  Even if
     * the padding of both matrices differ, the content will be
     * appropriately copied.
     */
    template<class OAlloc>
    void fill(const Matrix<T,OAlloc>& _other);
    
    /**
     * Check if the matrix is empty (zero rows or columns)
     */
    inline bool empty() const {
      return ( (this->_impl._rows==0) ||
               (this->_impl._cols==0) );
    }
    
    /**
     * Number of rows
     */
    inline size_t rows() const {
      return this->_impl._rows;
    }

    /**
     * Number of matrix columns (without any padding)
     */
    inline size_t cols() const {
      return this->_impl._cols;
    }

    /**
     * Number of number of columns including padding
     */
    inline size_t dcols() const {
      return this->_impl._dcols;
    }

    
    /**
     * Total number of entries (rows x cols)
     */
    inline size_t entries() const {
      return this->_impl._rows * this->_impl._cols;
    }

    /**
     * Pointer to data block
     */
    inline T* data() { return this->_impl._data; }

    /**
     * Pointer to data block
     */
    inline const T* data() const { return this->_impl._data; }

    /**
     * Extract one particular column
     *
     * This method has to copy the column, and hence it is relatively slow
     */
    inline std::vector<value_type> column(const size_t col) const;
    
    /**
     * @name Arithmetic operators
     */
    //@{

    /// Sum this and another matrix, and leave the result in here
    Matrix& operator+=(const Matrix& other);

    /// Subtract another matrix to this one, and leave the result in here
    Matrix& operator-=(const Matrix& other);
    
    //@}

  private:

    // Call the memory deallocation 
    void _deallocate();

    
    /// Use the allocator to create the necessary storage
    void _create_storage(size_t _rows,size_t _cols);

    /// Read-writable reference to the allocator in use
    allocator_type& _get_allocator() noexcept;

    /// Read-only reference to the allocator in use
    const allocator_type& _get_allocator() const noexcept;    

    /**
     * Fill with an initializer list
     */
    void fill(const std::initializer_list<std::initializer_list<value_type> >&);
  }; // class Matrix


  // External arithmetic operators
  template<typename T,class Alloc>
  Matrix<T,Alloc> operator+(const Matrix<T,Alloc>& a,
                            const Matrix<T,Alloc>& b);

  template<typename T,class Alloc>
  Matrix<T,Alloc> operator-(const Matrix<T,Alloc>& a,
                            const Matrix<T,Alloc>& b);
  
} // namespace ANPI

// include the template implementations
#include "Matrix.tpp"

#endif
