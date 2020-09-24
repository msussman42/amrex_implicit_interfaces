/* adept_source.h - Source code for the Adept library

  Copyright (C) 2012-2015 The University of Reading
  Copyright (C) 2015-2017 European Centre for Medium-Range Weather Forecasts

  Licensed under the Apache License, Version 2.0 (the "License"); you
  may not use this file except in compliance with the License.  You
  may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
  implied.  See the License for the specific language governing
  permissions and limitations under the License.


  This file was created automatically by script ./create_adept_source_header 
  on Sun 28 Jan 21:05:14 GMT 2018

  It contains a concatenation of the source files from the Adept
  library. The idea is that a program may #include this file in one of
  its source files (typically the one containing the main function),
  and then the Adept library will be built into the executable without
  the need to link to an external library. All other source files
  should just #include <adept.h> or <adept_arrays.h>. The ability to
  use Adept in this way makes it easier to distribute an Adept package
  that is usable on non-Unix platforms that are unable to use the
  autoconf configure script to build external libraries.

  If HAVE_BLAS is defined below then matrix multiplication will be
  enabled; the BLAS library should be provided at the link stage
  although no header file is required.  If HAVE_LAPACK is defined
  below then linear algebra routines will be enabled (matrix inverse
  and solving linear systems of equations); again, the LAPACK library
  should be provided at the link stage although no header file is
  required.

*/

/* Feel free to delete this warning: */
#ifdef _MSC_FULL_VER 
#pragma message("warning: the adept_source.h header file has not been edited so BLAS matrix multiplication and LAPACK linear-algebra support have been disabled")
#else
#warning "The adept_source.h header file has not been edited so BLAS matrix multiplication and LAPACK linear-algebra support have been disabled"
#endif

/* Uncomment this if you are linking to the BLAS library (header file
   not required) to enable matrix multiplication */
//#define HAVE_BLAS 1

/* Uncomment this if you are linking to the LAPACK library (header
   file not required) */
//#define HAVE_LAPACK 1

/* Uncomment this if you have the cblas.h header from OpenBLAS */
//#define HAVE_OPENBLAS_CBLAS_HEADER

/*

  The individual source files now follow.

*/

#ifndef AdeptSource_H
#define AdeptSource_H 1




// =================================================================
// Contents of config_platform_independent.h
// =================================================================

/* config_platform_independent.h.  Generated from config_platform_independent.h.in by configure.  */
/* config_platform_independent.h.in. */

/* Name of package */
#define PACKAGE "adept"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "r.j.hogan@ecmwf.int"

/* Define to the full name of this package. */
#define PACKAGE_NAME "adept"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "adept 2.0.5"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "adept"

/* Define to the home page for this package. */
#define PACKAGE_URL "http://www.met.reading.ac.uk/clouds/adept/"

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.0.5"

/* Version number of package */
#define VERSION "2.0.5"



// =================================================================
// Contents of cpplapack.h
// =================================================================

/* cpplapack.h -- C++ interface to LAPACK

    Copyright (C) 2015-2016 European Centre for Medium-Range Weather Forecasts

    Author: Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.
*/

#ifndef AdeptCppLapack_H
#define AdeptCppLapack_H 1                       

#include <vector>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LAPACK

extern "C" {
  // External LAPACK Fortran functions
  void sgetrf_(const int* m, const int* n, float*  a, const int* lda, int* ipiv, int* info);
  void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
  void sgetri_(const int* n, float* a, const int* lda, const int* ipiv, 
	       float* work, const int* lwork, int* info);
  void dgetri_(const int* n, double* a, const int* lda, const int* ipiv, 
	       double* work, const int* lwork, int* info);
  void ssytrf_(const char* uplo, const int* n, float* a, const int* lda, int* ipiv,
	       float* work, const int* lwork, int* info);
  void dsytrf_(const char* uplo, const int* n, double* a, const int* lda, int* ipiv,
	       double* work, const int* lwork, int* info);
  void ssytri_(const char* uplo, const int* n, float* a, const int* lda, 
	       const int* ipiv, float* work, int* info);
  void dsytri_(const char* uplo, const int* n, double* a, const int* lda, 
	       const int* ipiv, double* work, int* info);
  void ssysv_(const char* uplo, const int* n, const int* nrhs, float* a, const int* lda, 
	      int* ipiv, float* b, const int* ldb, float* work, const int* lwork, int* info);
  void dsysv_(const char* uplo, const int* n, const int* nrhs, double* a, const int* lda, 
	      int* ipiv, double* b, const int* ldb, double* work, const int* lwork, int* info);
  void sgesv_(const int* n, const int* nrhs, float* a, const int* lda, 
	      int* ipiv, float* b, const int* ldb, int* info);
  void dgesv_(const int* n, const int* nrhs, double* a, const int* lda, 
	      int* ipiv, double* b, const int* ldb, int* info);
}

namespace adept {

  // Overloaded functions provide both single &
  // double precision versions, and prevents the huge lapacke.h having
  // to be included in all user code
  namespace internal {
    typedef int lapack_int;
    // Factorize a general matrix
    inline
    int cpplapack_getrf(int n, float* a,  int lda, int* ipiv) {
      int info;
      sgetrf_(&n, &n, a, &lda, ipiv, &info);
      return info;
    }
    inline
    int cpplapack_getrf(int n, double* a, int lda, int* ipiv) {
      int info;
      dgetrf_(&n, &n, a, &lda, ipiv, &info);
      return info;
    }

    // Invert a general matrix
    inline
    int cpplapack_getri(int n, float* a,  int lda, const int* ipiv) {
      int info;
      float work_query;
      int lwork = -1;
      // Find out how much work memory required
      sgetri_(&n, a, &lda, ipiv, &work_query, &lwork, &info);
      lwork = static_cast<int>(work_query);
      std::vector<float> work(static_cast<size_t>(lwork));
      // Do full calculation
      sgetri_(&n, a, &lda, ipiv, &work[0], &lwork, &info);
      return info;
    }
    inline
    int cpplapack_getri(int n, double* a,  int lda, const int* ipiv) {
      int info;
      double work_query;
      int lwork = -1;
      // Find out how much work memory required
      dgetri_(&n, a, &lda, ipiv, &work_query, &lwork, &info);
      lwork = static_cast<int>(work_query);
      std::vector<double> work(static_cast<size_t>(lwork));
      // Do full calculation
      dgetri_(&n, a, &lda, ipiv, &work[0], &lwork, &info);
      return info;
    }

    // Factorize a symmetric matrix
    inline
    int cpplapack_sytrf(char uplo, int n, float* a, int lda, int* ipiv) {
      int info;
      float work_query;
      int lwork = -1;
      // Find out how much work memory required
      ssytrf_(&uplo, &n, a, &lda, ipiv, &work_query, &lwork, &info);
      lwork = static_cast<int>(work_query);
      std::vector<float> work(static_cast<size_t>(lwork));
      // Do full calculation
      ssytrf_(&uplo, &n, a, &lda, ipiv, &work[0], &lwork, &info);
      return info;
    }
    inline
    int cpplapack_sytrf(char uplo, int n, double* a, int lda, int* ipiv) {
      int info;
      double work_query;
      int lwork = -1;
      // Find out how much work memory required
      dsytrf_(&uplo, &n, a, &lda, ipiv, &work_query, &lwork, &info);
      lwork = static_cast<int>(work_query);
      std::vector<double> work(static_cast<size_t>(lwork));
      // Do full calculation
      dsytrf_(&uplo, &n, a, &lda, ipiv, &work[0], &lwork, &info);
      return info;
    }

    // Invert a symmetric matrix
    inline
    int cpplapack_sytri(char uplo, int n, float* a, int lda, const int* ipiv) {
      int info;
      std::vector<float> work(n);
      ssytri_(&uplo, &n, a, &lda, ipiv, &work[0], &info);
      return info;
    }
    inline
    int cpplapack_sytri(char uplo, int n, double* a, int lda, const int* ipiv) {
      int info;
      std::vector<double> work(n);
      dsytri_(&uplo, &n, a, &lda, ipiv, &work[0], &info);
      return info;
    }

    // Solve system of linear equations with general matrix
    inline
    int cpplapack_gesv(int n, int nrhs, float* a, int lda,
		       int* ipiv, float* b, int ldb) {
      int info;
      sgesv_(&n, &nrhs, a, &lda, ipiv, b, &lda, &info);
      return info;
    }
    inline
    int cpplapack_gesv(int n, int nrhs, double* a, int lda,
		       int* ipiv, double* b, int ldb) {
      int info;
      dgesv_(&n, &nrhs, a, &lda, ipiv, b, &lda, &info);
      return info;
    }

    // Solve system of linear equations with symmetric matrix
    inline
    int cpplapack_sysv(char uplo, int n, int nrhs, float* a, int lda, int* ipiv,
		       float* b, int ldb) {
      int info;
      float work_query;
      int lwork = -1;
      // Find out how much work memory required
      ssysv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &work_query, &lwork, &info);
      lwork = static_cast<int>(work_query);
      std::vector<float> work(static_cast<size_t>(lwork));
      // Do full calculation
      ssysv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &work[0], &lwork, &info);
      return info;
    }
    inline
    int cpplapack_sysv(char uplo, int n, int nrhs, double* a, int lda, int* ipiv,
		       double* b, int ldb) {
      int info;
      double work_query;
      int lwork = -1;
      // Find out how much work memory required
      dsysv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &work_query, &lwork, &info);
      lwork = static_cast<int>(work_query);
      std::vector<double> work(static_cast<size_t>(lwork));
      // Do full calculation
      dsysv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &work[0], &lwork, &info);
      return info;
    }

  }
}

#endif

#endif


// =================================================================
// Contents of Array.cpp
// =================================================================

/* Array.cpp -- Functions and global variables controlling array behaviour

    Copyright (C) 2015-2016 European Centre for Medium-Range Weather Forecasts

    Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.
*/


#include <adept/Array.h>

namespace adept {
  namespace internal {
    bool array_row_major_order = true;
    //    bool array_print_curly_brackets = true;

    // Variables describing how arrays are written to a stream
    ArrayPrintStyle array_print_style = PRINT_STYLE_CURLY;
    std::string vector_separator = ", ";
    std::string vector_print_before = "{";
    std::string vector_print_after = "}";
    std::string array_opening_bracket = "{";
    std::string array_closing_bracket = "}";
    std::string array_contiguous_separator = ", ";
    std::string array_non_contiguous_separator = ",\n";
    std::string array_print_before = "\n{";
    std::string array_print_after = "}";
    std::string array_print_empty_before = "(empty rank-";
    std::string array_print_empty_after = " array)";
    bool array_print_indent = true;
    bool array_print_empty_rank = true;
  }

  void set_array_print_style(ArrayPrintStyle ps) {
    using namespace internal;
    switch (ps) {
    case PRINT_STYLE_PLAIN:
       vector_separator = " ";
       vector_print_before = "";
       vector_print_after = "";
       array_opening_bracket = "";
       array_closing_bracket = "";
       array_contiguous_separator = " ";
       array_non_contiguous_separator = "\n";
       array_print_before = "";
       array_print_after = "";
       array_print_empty_before = "(empty rank-";
       array_print_empty_after = " array)";
       array_print_indent = false;
       array_print_empty_rank = true;
       break;
    case PRINT_STYLE_CSV:
       vector_separator = ", ";
       vector_print_before = "";
       vector_print_after = "";
       array_opening_bracket = "";
       array_closing_bracket = "";
       array_contiguous_separator = ", ";
       array_non_contiguous_separator = "\n";
       array_print_before = "";
       array_print_after = "";
       array_print_empty_before = "empty";
       array_print_empty_after = "";
       array_print_indent = false;
       array_print_empty_rank = false;
       break;
    case PRINT_STYLE_MATLAB:
       vector_separator = " ";
       vector_print_before = "[";
       vector_print_after = "]";
       array_opening_bracket = "";
       array_closing_bracket = "";
       array_contiguous_separator = " ";
       array_non_contiguous_separator = ";\n";
       array_print_before = "[";
       array_print_after = "]";
       array_print_empty_before = "[";
       array_print_empty_after = "]";
       array_print_indent = true;
       array_print_empty_rank = false;
       break;
    case PRINT_STYLE_CURLY:
       vector_separator = ", ";
       vector_print_before = "{";
       vector_print_after = "}";
       array_opening_bracket = "{";
       array_closing_bracket = "}";
       array_contiguous_separator = ", ";
       array_non_contiguous_separator = ",\n";
       array_print_before = "\n{";
       array_print_after = "}";
       array_print_empty_before = "(empty rank-";
       array_print_empty_after = " array)";
       array_print_indent = true;
       array_print_empty_rank = true;
       break;
    default:
      throw invalid_operation("Array print style not understood");
    }
    array_print_style = ps;
  }

}


// =================================================================
// Contents of Stack.cpp
// =================================================================

/* Stack.cpp -- Stack for storing automatic differentiation information

     Copyright (C) 2012-2014 University of Reading
    Copyright (C) 2015 European Centre for Medium-Range Weather Forecasts

    Author: Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.

*/


#include <iostream>
#include <cstring> // For memcpy


#ifdef _OPENMP
#include <omp.h>
#endif

#include <adept/Stack.h>


namespace adept {

  using namespace internal;

  // Global pointers to the current thread, the second of which is
  // thread safe. The first is only used if ADEPT_STACK_THREAD_UNSAFE
  // is defined.
  ADEPT_THREAD_LOCAL Stack* _stack_current_thread = 0;
  Stack* _stack_current_thread_unsafe = 0;

  // MEMBER FUNCTIONS OF THE STACK CLASS

  // Destructor: frees dynamically allocated memory (if any)
  Stack::~Stack() {
    // If this is the currently active stack then set to NULL as
    // "this" is shortly to become invalid
    if (is_thread_unsafe_) {
      if (_stack_current_thread_unsafe == this) {
	_stack_current_thread_unsafe = 0; 
      }
    }
    else if (_stack_current_thread == this) {
      _stack_current_thread = 0; 
    }
#ifndef ADEPT_STACK_STORAGE_STL
    if (gradient_) {
      delete[] gradient_;
    }
#endif
  }
  
  // Make this stack "active" by copying its "this" pointer to a
  // global variable; this makes it the stack that aReal objects
  // subsequently interact with when being created and participating
  // in mathematical expressions
  void
  Stack::activate()
  {
    // Check that we don't already have an active stack in this thread
    if ((is_thread_unsafe_ && _stack_current_thread_unsafe 
	 && _stack_current_thread_unsafe != this)
	|| ((!is_thread_unsafe_) && _stack_current_thread
	    && _stack_current_thread != this)) {
      throw(stack_already_active());
    }
    else {
      if (!is_thread_unsafe_) {
	_stack_current_thread = this;
      }
      else {
	_stack_current_thread_unsafe = this;
      }
    }    
  }

  
  // Set the maximum number of threads to be used in Jacobian
  // calculations, if possible. A value of 1 indicates that OpenMP
  // will not be used, while a value of 0 indicates that the number
  // will match the number of available processors. Returns the
  // maximum that will be used, which will be 1 if the Adept library
  // was compiled without OpenMP support. Note that a value of 1 will
  // disable the use of OpenMP with Adept, so Adept will then use no
  // OpenMP directives or function calls. Note that if in your program
  // you use OpenMP with each thread performing automatic
  // differentiaion with its own independent Adept stack, then
  // typically only one OpenMP thread is available for each Jacobian
  // calculation, regardless of whether you call this function.
  int
  Stack::set_max_jacobian_threads(int n)
  {
#ifdef _OPENMP
    if (have_openmp_) {
      if (n == 1) {
	openmp_manually_disabled_ = true;
	return 1;
      }
      else if (n < 1) {
	openmp_manually_disabled_ = false;
	omp_set_num_threads(omp_get_num_procs());
	return omp_get_max_threads();
      }
      else {
	openmp_manually_disabled_ = false;
	omp_set_num_threads(n);
	return omp_get_max_threads();
      }
    }
#endif
    return 1;
  }


  // Return maximum number of OpenMP threads to be used in Jacobian
  // calculation
  int 
  Stack::max_jacobian_threads() const
  {
#ifdef _OPENMP
    if (have_openmp_) {
      if (openmp_manually_disabled_) {
	return 1;
      }
      else {
	return omp_get_max_threads();
      }
    }
#endif
    return 1;
  }


  // Perform to adjoint computation (reverse mode). It is assumed that
  // some gradients have been assigned already, otherwise the function
  // returns with an error.
  void
  Stack::compute_adjoint()
  {
    if (gradients_are_initialized()) {
      // Loop backwards through the derivative statements
      for (uIndex ist = n_statements_-1; ist > 0; ist--) {
	const Statement& statement = statement_[ist];
	// We copy the RHS gradient (LHS in the original derivative
	// statement but swapped in the adjoint equivalent) to "a" in
	// case it appears on the LHS in any of the following statements
	Real a = gradient_[statement.index];
	gradient_[statement.index] = 0.0;
	// By only looping if a is non-zero we gain a significant speed-up
	if (a != 0.0) {
	  // Loop over operations
	  for (uIndex i = statement_[ist-1].end_plus_one;
	       i < statement.end_plus_one; i++) {
	    gradient_[index_[i]] += multiplier_[i]*a;
	  }
	}
      }
    }  
    else {
      throw(gradients_not_initialized());
    }  
  }


  // Perform tangent linear computation (forward mode). It is assumed
  // that some gradients have been assigned already, otherwise the
  // function returns with an error.
  void
  Stack::compute_tangent_linear()
  {
    if (gradients_are_initialized()) {
      // Loop forward through the statements
      for (uIndex ist = 1; ist < n_statements_; ist++) {
	const Statement& statement = statement_[ist];
	// We copy the LHS to "a" in case it appears on the RHS in any
	// of the following statements
	Real a = 0.0;
	for (uIndex i = statement_[ist-1].end_plus_one;
	     i < statement.end_plus_one; i++) {
	  a += multiplier_[i]*gradient_[index_[i]];
	}
	gradient_[statement.index] = a;
      }
    }
    else {
      throw(gradients_not_initialized());
    }
  }



  // Register n gradients
  uIndex
  Stack::do_register_gradients(const uIndex& n) {
    n_gradients_registered_ += n;
    if (!gap_list_.empty()) {
      uIndex return_val;
      // Insert in a gap, if there is one big enough
      for (GapListIterator it = gap_list_.begin();
	   it != gap_list_.end(); it++) {
	uIndex len = it->end + 1 - it->start;
	if (len > n) {
	  // Gap a bit larger than needed: reduce its size
	  return_val = it->start;
	  it->start += n;
	  return return_val;
	}
	else if (len == n) {
	  // Gap exactly the size needed: fill it and remove from list
	  return_val = it->start;
	  if (most_recent_gap_ == it) {
	    gap_list_.erase(it);
	    most_recent_gap_ = gap_list_.end();
	  }
	  else {
	    gap_list_.erase(it);
	  }
	  return return_val;
	}
      }
    }
    // No suitable gap found; instead add to end of gradient vector
    i_gradient_ += n;
    if (i_gradient_ > max_gradient_) {
      max_gradient_ = i_gradient_;
    }
    return i_gradient_ - n;
  }
  

  // If an aReal object is deleted, its gradient_index is
  // unregistered from the stack.  If this is at the top of the stack
  // then this is easy and is done inline; this is the usual case
  // since C++ trys to deallocate automatic objects in the reverse
  // order to that in which they were allocated.  If it is not at the
  // top of the stack then a non-inline function is called to ensure
  // that the gap list is adjusted correctly.
  void
  Stack::unregister_gradient_not_top(const uIndex& gradient_index)
  {
    enum {
      ADDED_AT_BASE,
      ADDED_AT_TOP,
      NEW_GAP,
      NOT_FOUND
    } status = NOT_FOUND;
    // First try to find if the unregistered element is at the
    // start or end of an existing gap
    if (!gap_list_.empty() && most_recent_gap_ != gap_list_.end()) {
      // We have a "most recent" gap - check whether the gradient
      // to be unregistered is here
      Gap& current_gap = *most_recent_gap_;
      if (gradient_index == current_gap.start - 1) {
	current_gap.start--;
	status = ADDED_AT_BASE;
      }
      else if (gradient_index == current_gap.end + 1) {
	current_gap.end++;
	status = ADDED_AT_TOP;
      }
      // Should we check for erroneous removal from middle of gap?
    }
    if (status == NOT_FOUND) {
      // Search other gaps
      for (GapListIterator it = gap_list_.begin();
	   it != gap_list_.end(); it++) {
	if (gradient_index <= it->end + 1) {
	  // Gradient to unregister is either within the gap
	  // referenced by iterator "it", or it is between "it"
	  // and the previous gap in the list
	  if (gradient_index == it->start - 1) {
	    status = ADDED_AT_BASE;
	    it->start--;
	    most_recent_gap_ = it;
	  }
	  else if (gradient_index == it->end + 1) {
	    status = ADDED_AT_TOP;
	    it->end++;
	    most_recent_gap_ = it;
	  }
	  else {
	    // Insert a new gap of width 1; note that list::insert
	    // inserts *before* the specified location
	    most_recent_gap_
	      = gap_list_.insert(it, Gap(gradient_index));
	    status = NEW_GAP;
	  }
	  break;
	}
      }
      if (status == NOT_FOUND) {
	gap_list_.push_back(Gap(gradient_index));
	most_recent_gap_ = gap_list_.end();
	most_recent_gap_--;
      }
    }
    // Finally check if gaps have merged
    if (status == ADDED_AT_BASE
	&& most_recent_gap_ != gap_list_.begin()) {
      // Check whether the gap has merged with the next one
      GapListIterator it = most_recent_gap_;
      it--;
      if (it->end == most_recent_gap_->start - 1) {
	// Merge two gaps
	most_recent_gap_->start = it->start;
	gap_list_.erase(it);
      }
    }
    else if (status == ADDED_AT_TOP) {
      GapListIterator it = most_recent_gap_;
      it++;
      if (it != gap_list_.end()
	  && it->start == most_recent_gap_->end + 1) {
	// Merge two gaps
	most_recent_gap_->end = it->end;
	gap_list_.erase(it);
      }
    }
  }	


  // Unregister n gradients starting at gradient_index
  void
  Stack::unregister_gradients(const uIndex& gradient_index,
			      const uIndex& n)
  {
    n_gradients_registered_ -= n;
    if (gradient_index+n == i_gradient_) {
      // Gradient to be unregistered is at the top of the stack
      i_gradient_ -= n;
      if (!gap_list_.empty()) {
	Gap& last_gap = gap_list_.back();
	if (i_gradient_ == last_gap.end+1) {
	  // We have unregistered the elements between the "gap" of
	  // unregistered element and the top of the stack, so can set
	  // the variables indicating the presence of the gap to zero
	  i_gradient_ = last_gap.start;
	  GapListIterator it = gap_list_.end();
	  it--;
	  if (most_recent_gap_ == it) {
	    most_recent_gap_ = gap_list_.end();
	  }
	  gap_list_.pop_back();
	}
      }
    }
    else { // Gradients to be unregistered not at top of stack.
      enum {
	ADDED_AT_BASE,
	ADDED_AT_TOP,
	NEW_GAP,
	NOT_FOUND
      } status = NOT_FOUND;
      // First try to find if the unregistered element is at the start
      // or end of an existing gap
      if (!gap_list_.empty() && most_recent_gap_ != gap_list_.end()) {
	// We have a "most recent" gap - check whether the gradient
	// to be unregistered is here
	Gap& current_gap = *most_recent_gap_;
	if (gradient_index == current_gap.start - n) {
	  current_gap.start -= n;
	  status = ADDED_AT_BASE;
	}
	else if (gradient_index == current_gap.end + 1) {
	  current_gap.end += n;
	  status = ADDED_AT_TOP;
	}
	/*
	else if (gradient_index > current_gap.start - n
		 && gradient_index < current_gap.end + 1) {
	  std::cout << "** Attempt to find " << gradient_index << " in gaps ";
	  print_gaps();
	  std::cout << "\n";
	  throw invalid_operation("Gap list corruption");
	}
	*/
	// Should we check for erroneous removal from middle of gap?
      }
      if (status == NOT_FOUND) {
	// Search other gaps
	for (GapListIterator it = gap_list_.begin();
	     it != gap_list_.end(); it++) {
	  if (gradient_index <= it->end + 1) {
	    // Gradient to unregister is either within the gap
	    // referenced by iterator "it", or it is between "it" and
	    // the previous gap in the list
	    if (gradient_index == it->start - n) {
	      status = ADDED_AT_BASE;
	      it->start -= n;
	      most_recent_gap_ = it;
	    }
	    else if (gradient_index == it->end + 1) {
	      status = ADDED_AT_TOP;
	      it->end += n;
	      most_recent_gap_ = it;
	    }
	    /*
	    else if (gradient_index > it->start - n) {
	      std::cout << "*** Attempt to find " << gradient_index << " in gaps ";
	      print_gaps();
	      std::cout << "\n";
	      throw invalid_operation("Gap list corruption");
	    }
	    */
	    else {
	      // Insert a new gap; note that list::insert inserts
	      // *before* the specified location
	      most_recent_gap_
		= gap_list_.insert(it, Gap(gradient_index,
					   gradient_index+n-1));
	      status = NEW_GAP;
	    }
	    break;
	  }
	}
	if (status == NOT_FOUND) {
	  gap_list_.push_back(Gap(gradient_index,
				  gradient_index+n-1));
	  most_recent_gap_ = gap_list_.end();
	  most_recent_gap_--;
	}
      }
      // Finally check if gaps have merged
      if (status == ADDED_AT_BASE
	  && most_recent_gap_ != gap_list_.begin()) {
	// Check whether the gap has merged with the next one
	GapListIterator it = most_recent_gap_;
	it--;
	if (it->end == most_recent_gap_->start - 1) {
	  // Merge two gaps
	  most_recent_gap_->start = it->start;
	  gap_list_.erase(it);
	}
      }
      else if (status == ADDED_AT_TOP) {
	GapListIterator it = most_recent_gap_;

	it++;
	if (it != gap_list_.end()
	    && it->start == most_recent_gap_->end + 1) {
	  // Merge two gaps
	  most_recent_gap_->end = it->end;
	  gap_list_.erase(it);
	}
      }
    }
  }
  
  
  // Print each derivative statement to the specified stream (standard
  // output if omitted)
  void
  Stack::print_statements(std::ostream& os) const
  {
    for (uIndex ist = 1; ist < n_statements_; ist++) {
      const Statement& statement = statement_[ist];
      os << ist
		<< ": d[" << statement.index
		<< "] = ";
      
      if (statement_[ist-1].end_plus_one == statement_[ist].end_plus_one) {
	os << "0\n";
      }
      else {    
	for (uIndex i = statement_[ist-1].end_plus_one;
	     i < statement.end_plus_one; i++) {
	  os << " + " << multiplier_[i] << "*d[" << index_[i] << "]";
	}
	os << "\n";
      }
    }
  }
  
  // Print the current gradient list to the specified stream (standard
  // output if omitted)
  bool
  Stack::print_gradients(std::ostream& os) const
  {
    if (gradients_are_initialized()) {
      for (uIndex i = 0; i < max_gradient_; i++) {
	if (i%10 == 0) {
	  if (i != 0) {
	    os << "\n";
	  }
	  os << i << ":";
	}
	os << " " << gradient_[i];
      }
      os << "\n";
      return true;
    }
    else {
      os << "No gradients initialized\n";
      return false;
    }
  }

  // Print the list of gaps in the gradient list to the specified
  // stream (standard output if omitted)
  void
  Stack::print_gaps(std::ostream& os) const
  {
    for (std::list<Gap>::const_iterator it = gap_list_.begin();
	 it != gap_list_.end(); it++) {
      os << it->start << "-" << it->end << " ";
    }
  }


#ifndef ADEPT_STACK_STORAGE_STL
  // Initialize the vector of gradients ready for the adjoint
  // calculation
  void
  Stack::initialize_gradients()
  {
    if (max_gradient_ > 0) {
      if (n_allocated_gradients_ < max_gradient_) {
	if (gradient_) {
	  delete[] gradient_;
	}
	gradient_ = new Real[max_gradient_];
	n_allocated_gradients_ = max_gradient_;
      }
      for (uIndex i = 0; i < max_gradient_; i++) {
	gradient_[i] = 0.0;
      }
    }
    gradients_initialized_ = true;
  }
#else
  void
  Stack::initialize_gradients()
  {
    gradient_.resize(max_gradient_+10, 0.0);
      gradients_initialized_ = true;
  }
#endif

  // Report information about the stack to the specified stream, or
  // standard output if omitted; note that this is synonymous with
  // sending the Stack object to a stream using the "<<" operator.
  void
  Stack::print_status(std::ostream& os) const
  {
    os << "Automatic Differentiation Stack (address " << this << "):\n";
    if ((!is_thread_unsafe_) && _stack_current_thread == this) {
      os << "   Currently attached - thread safe\n";
    }
    else if (is_thread_unsafe_ && _stack_current_thread_unsafe == this) {
      os << "   Currently attached - thread unsafe\n";
    }
    else {
      os << "   Currently detached\n";
    }
    os << "   Recording status:\n";
    if (is_recording_) {
      os << "      Recording is ON\n";  
    }
    else {
      os << "      Recording is PAUSED\n";
    }
    // Account for the null statement at the start by subtracting one
    os << "      " << n_statements()-1 << " statements (" 
       << n_allocated_statements() << " allocated)";
    os << " and " << n_operations() << " operations (" 
       << n_allocated_operations() << " allocated)\n";
    os << "      " << n_gradients_registered() << " gradients currently registered ";
    os << "and a total of " << max_gradients() << " needed (current index "
       << i_gradient() << ")\n";
    if (gap_list_.empty()) {
      os << "      Gradient list has no gaps\n";
    }
    else {
      os << "      Gradient list has " << gap_list_.size() << " gaps (";
      print_gaps(os);
      os << ")\n";
    }
    os << "   Computation status:\n";
    if (gradients_are_initialized()) {
      os << "      " << max_gradients() << " gradients assigned (" 
	 << n_allocated_gradients() << " allocated)\n";
    }
    else {
      os << "      0 gradients assigned (" << n_allocated_gradients()
	 << " allocated)\n";
    }
    os << "      Jacobian size: " << n_dependents() << "x" << n_independents() << "\n";
    if (n_dependents() <= 10 && n_independents() <= 10) {
      os << "      Independent indices:";
      for (std::size_t i = 0; i < independent_index_.size(); ++i) {
	os << " " << independent_index_[i];
      }
      os << "\n      Dependent indices:  ";
      for (std::size_t i = 0; i < dependent_index_.size(); ++i) {
	os << " " << dependent_index_[i];
      }
      os << "\n";
    }

#ifdef _OPENMP
    if (have_openmp_) {
      if (openmp_manually_disabled_) {
	os << "      Parallel Jacobian calculation manually disabled\n";
      }
      else {
	os << "      Parallel Jacobian calculation can use up to "
	   << omp_get_max_threads() << " threads\n";
	os << "      Each thread treats " << ADEPT_MULTIPASS_SIZE 
	   << " (in)dependent variables\n";
      }
    }
    else {
#endif
      os << "      Parallel Jacobian calculation not available\n";
#ifdef _OPENMP
    }
#endif
  }
} // End namespace adept



// =================================================================
// Contents of StackStorageOrig.cpp
// =================================================================

/* StackStorageOrig.cpp -- Original storage of stacks using STL containers

    Copyright (C) 2014-2015 University of Reading

    Author: Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.

   The Stack class inherits from a class providing the storage (and
   interface to the storage) for the derivative statements that are
   accumulated during the execution of an algorithm.  The derivative
   statements are held in two stacks described by Hogan (2014): the
   "statement stack" and the "operation stack".

   This file provides one of the original storage engine, which used
   std::vector to hold the two stacks. Note that these stacks are
   contiguous in memory, which is not ideal for very large algorithms.

*/

#include <cstring>

#include <adept/StackStorageOrig.h>

namespace adept {
  namespace internal {

    StackStorageOrig::~StackStorageOrig() {
      if (statement_) {
	delete[] statement_;
      }
      if (multiplier_) {
	delete[] multiplier_;
      }
      if (index_) {
	delete[] index_;
      }
    }


    // Double the size of the operation stack, or grow it even more if
    // the requested minimum number of extra entries (min) is greater
    // than this would allow
    void
    StackStorageOrig::grow_operation_stack(uIndex min)
    {
      uIndex new_size = 2*n_allocated_operations_;
      if (min > 0 && new_size < n_allocated_operations_+min) {
	new_size += min;
      }
      Real* new_multiplier = new Real[new_size];
      uIndex* new_index = new uIndex[new_size];
      
      std::memcpy(new_multiplier, multiplier_, n_operations_*sizeof(Real));
      std::memcpy(new_index, index_, n_operations_*sizeof(uIndex));
      
      delete[] multiplier_;
      delete[] index_;
      
      multiplier_ = new_multiplier;
      index_ = new_index;
      
      n_allocated_operations_ = new_size;
    }
    
    // ... likewise for the statement stack
    void
    StackStorageOrig::grow_statement_stack(uIndex min)
    {
      uIndex new_size = 2*n_allocated_statements_;
      if (min > 0 && new_size < n_allocated_statements_+min) {
	new_size += min;
      }
      Statement* new_statement = new Statement[new_size];
      std::memcpy(new_statement, statement_,
		  n_statements_*sizeof(Statement));
      delete[] statement_;
      
      statement_ = new_statement;
      
      n_allocated_statements_ = new_size;
    }

  }
}


// =================================================================
// Contents of Storage.cpp
// =================================================================

/* Storage.cpp -- Global variables recording use of Storage objects

    Copyright (C) 2015 European Centre for Medium-Range Weather Forecasts

    Author: Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.

*/

#include <adept/Storage.h>

namespace adept {
  namespace internal {
    Index n_storage_objects_created_;
    Index n_storage_objects_deleted_;
  }
}


// =================================================================
// Contents of cppblas.cpp
// =================================================================

/* cppblas.cpp -- C++ interface to BLAS functions

    Copyright (C) 2015-2016 European Centre for Medium-Range Weather Forecasts

    Author: Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.

   This file provides a C++ interface to selected Level-2 and -3 BLAS
   functions in which the precision of the arguments (float versus
   double) is inferred via overloading

*/

#include <adept/exception.h>
#include <adept/cppblas.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_BLAS

extern "C" {
  void sgemm_(const char* TransA, const char* TransB, const int* M,
	      const int* N, const int* K, const float* alpha,
	      const float* A, const int* lda, const float* B, const int* ldb,
	      const float* beta, const float* C, const int* ldc);
  void dgemm_(const char* TransA, const char* TransB, const int* M,
	      const int* N, const int* K, const double* alpha,
	      const double* A, const int* lda, const double* B, const int* ldb,
	      const double* beta, const double* C, const int* ldc);
  void sgemv_(const char* TransA, const int* M, const int* N, const float* alpha,
	      const float* A, const int* lda, const float* X, const int* incX,
	      const float* beta, const float* Y, const int* incY);
  void dgemv_(const char* TransA, const int* M, const int* N, const double* alpha,
	      const double* A, const int* lda, const double* X, const int* incX,
	      const double* beta, const double* Y, const int* incY);
  void ssymm_(const char* side, const char* uplo, const int* M, const int* N,
	      const float* alpha, const float* A, const int* lda, const float* B,
	      const int* ldb, const float* beta, float* C, const int* ldc);
  void dsymm_(const char* side, const char* uplo, const int* M, const int* N,
	      const double* alpha, const double* A, const int* lda, const double* B,
	      const int* ldb, const double* beta, double* C, const int* ldc);
  void ssymv_(const char* uplo, const int* N, const float* alpha, const float* A, 
	      const int* lda, const float* X, const int* incX, const float* beta, 
	      const float* Y, const int* incY);
  void dsymv_(const char* uplo, const int* N, const double* alpha, const double* A, 
	      const int* lda, const double* X, const int* incX, const double* beta, 
	      const double* Y, const int* incY);
  void sgbmv_(const char* TransA, const int* M, const int* N, const int* kl, 
	      const int* ku, const float* alpha, const float* A, const int* lda,
	      const float* X, const int* incX, const float* beta, 
	      const float* Y, const int* incY);
  void dgbmv_(const char* TransA, const int* M, const int* N, const int* kl, 
	      const int* ku, const double* alpha, const double* A, const int* lda,
	      const double* X, const int* incX, const double* beta, 
	      const double* Y, const int* incY);
};

namespace adept {

  namespace internal {
    
    // Matrix-matrix multiplication for general dense matrices
#define ADEPT_DEFINE_GEMM(T, FUNC, FUNC_COMPLEX)		\
    void cppblas_gemm(BLAS_ORDER Order,				\
		      BLAS_TRANSPOSE TransA,			\
		      BLAS_TRANSPOSE TransB,			\
		      int M, int N,				\
		      int K, T alpha, const T *A,		\
		      int lda, const T *B, int ldb,		\
		      T beta, T *C, int ldc) {			\
      if (Order == BlasColMajor) {				\
        FUNC(&TransA, &TransB, &M, &N, &K, &alpha, A, &lda,	\
	     B, &ldb, &beta, C, &ldc);				\
      }								\
      else {							\
        FUNC(&TransB, &TransA, &N, &M, &K, &alpha, B, &ldb,	\
	     A, &lda, &beta, C, &ldc);				\
      }								\
    }
    ADEPT_DEFINE_GEMM(double, dgemm_, zgemm_);
    ADEPT_DEFINE_GEMM(float,  sgemm_, cgemm_);
#undef ADEPT_DEFINE_GEMM
    
    // Matrix-vector multiplication for a general dense matrix
#define ADEPT_DEFINE_GEMV(T, FUNC, FUNC_COMPLEX)		\
    void cppblas_gemv(const BLAS_ORDER Order,			\
		      const BLAS_TRANSPOSE TransA,		\
		      const int M, const int N,			\
		      const T alpha, const T *A, const int lda,	\
		      const T *X, const int incX, const T beta,	\
		      T *Y, const int incY) {			\
      if (Order == BlasColMajor) {				\
        FUNC(&TransA, &M, &N, &alpha, A, &lda, X, &incX, 	\
	     &beta, Y, &incY);					\
      }								\
      else {							\
        BLAS_TRANSPOSE TransNew					\
	  = TransA == BlasTrans ? BlasNoTrans : BlasTrans;	\
        FUNC(&TransNew, &N, &M, &alpha, A, &lda, X, &incX, 	\
	     &beta, Y, &incY);					\
      }								\
    }
    ADEPT_DEFINE_GEMV(double, dgemv_, zgemv_);
    ADEPT_DEFINE_GEMV(float,  sgemv_, cgemv_);
#undef ADEPT_DEFINE_GEMV
    
    // Matrix-matrix multiplication where matrix A is symmetric
    // FIX! CHECK ROW MAJOR VERSION IS RIGHT			
#define ADEPT_DEFINE_SYMM(T, FUNC, FUNC_COMPLEX)			\
    void cppblas_symm(const BLAS_ORDER Order,				\
		      const BLAS_SIDE Side,				\
		      const BLAS_UPLO Uplo,				\
		      const int M, const int N,				\
		      const T alpha, const T *A, const int lda,		\
		      const T *B, const int ldb, const T beta,		\
		      T *C, const int ldc) {				\
      if (Order == BlasColMajor) {					\
        FUNC(&Side, &Uplo, &M, &N, &alpha, A, &lda,			\
	     B, &ldb, &beta, C, &ldc);					\
      }									\
      else {								\
	BLAS_SIDE SideNew = Side == BlasLeft  ? BlasRight : BlasLeft;	\
	BLAS_UPLO UploNew = Uplo == BlasUpper ? BlasLower : BlasUpper;  \
        FUNC(&SideNew, &UploNew, &N, &M, &alpha, A, &lda,		\
	     B, &ldb, &beta, C, &ldc);					\
      }									\
    }
    ADEPT_DEFINE_SYMM(double, dsymm_, zsymm_);
    ADEPT_DEFINE_SYMM(float,  ssymm_, csymm_);
#undef ADEPT_DEFINE_SYMM
    
    // Matrix-vector multiplication where the matrix is symmetric
#define ADEPT_DEFINE_SYMV(T, FUNC, FUNC_COMPLEX)			\
    void cppblas_symv(const BLAS_ORDER Order,				\
		      const BLAS_UPLO Uplo,				\
		      const int N, const T alpha, const T *A,		\
		      const int lda, const T *X, const int incX,	\
		      const T beta, T *Y, const int incY) {		\
      if (Order == BlasColMajor) {					\
        FUNC(&Uplo, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);	\
      }									\
      else {								\
        BLAS_UPLO UploNew = Uplo == BlasUpper ? BlasLower : BlasUpper;  \
        FUNC(&UploNew, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);	\
      }									\
    }
    ADEPT_DEFINE_SYMV(double, dsymv_, zsymv_);
    ADEPT_DEFINE_SYMV(float,  ssymv_, csymv_);
#undef ADEPT_DEFINE_SYMV
    
    // Matrix-vector multiplication for a general band matrix
#define ADEPT_DEFINE_GBMV(T, FUNC, FUNC_COMPLEX)		\
    void cppblas_gbmv(const BLAS_ORDER Order,			\
		      const BLAS_TRANSPOSE TransA,		\
		      const int M, const int N,			\
		      const int KL, const int KU, const T alpha,\
		      const T *A, const int lda, const T *X,	\
		      const int incX, const T beta, T *Y,	\
		      const int incY) {				\
      if (Order == BlasColMajor) {				\
        FUNC(&TransA, &M, &N, &KL, &KU, &alpha, A, &lda,	\
	     X, &incX, &beta, Y, &incY);			\
      }								\
      else {							\
	BLAS_TRANSPOSE TransNew					\
	  = TransA == BlasTrans ? BlasNoTrans : BlasTrans;	\
	FUNC(&TransNew, &N, &M, &KU, &KL, &alpha, A, &lda,	\
	     X, &incX, &beta, Y, &incY);			\
      }								\
    }
    ADEPT_DEFINE_GBMV(double, dgbmv_, zgbmv_);
    ADEPT_DEFINE_GBMV(float,  sgbmv_, cgbmv_);
#undef ADEPT_DEFINE_GBMV
  
  } // End namespace internal
  
} // End namespace adept
  

#else // Don't have BLAS


namespace adept {

  namespace internal {
    
    // Matrix-matrix multiplication for general dense matrices
#define ADEPT_DEFINE_GEMM(T, FUNC, FUNC_COMPLEX)		\
    void cppblas_gemm(BLAS_ORDER Order,				\
		      BLAS_TRANSPOSE TransA,			\
		      BLAS_TRANSPOSE TransB,			\
		      int M, int N,				\
		      int K, T alpha, const T *A,		\
		      int lda, const T *B, int ldb,		\
		      T beta, T *C, int ldc) {			\
      throw feature_not_available("Cannot perform matrix-matrix multiplication because compiled without BLAS"); \
    }
    ADEPT_DEFINE_GEMM(double, dgemm_, zgemm_);
    ADEPT_DEFINE_GEMM(float,  sgemm_, cgemm_);
#undef ADEPT_DEFINE_GEMM
    
    // Matrix-vector multiplication for a general dense matrix
#define ADEPT_DEFINE_GEMV(T, FUNC, FUNC_COMPLEX)		\
    void cppblas_gemv(const BLAS_ORDER Order,			\
		      const BLAS_TRANSPOSE TransA,		\
		      const int M, const int N,			\
		      const T alpha, const T *A, const int lda,	\
		      const T *X, const int incX, const T beta,	\
		      T *Y, const int incY) {			\
      throw feature_not_available("Cannot perform matrix-vector multiplication because compiled without BLAS"); \
    }
    ADEPT_DEFINE_GEMV(double, dgemv_, zgemv_);
    ADEPT_DEFINE_GEMV(float,  sgemv_, cgemv_);
#undef ADEPT_DEFINE_GEMV
    
    // Matrix-matrix multiplication where matrix A is symmetric
    // FIX! CHECK ROW MAJOR VERSION IS RIGHT			
#define ADEPT_DEFINE_SYMM(T, FUNC, FUNC_COMPLEX)			\
    void cppblas_symm(const BLAS_ORDER Order,				\
		      const BLAS_SIDE Side,				\
		      const BLAS_UPLO Uplo,				\
		      const int M, const int N,				\
		      const T alpha, const T *A, const int lda,		\
		      const T *B, const int ldb, const T beta,		\
		      T *C, const int ldc) {				\
      throw feature_not_available("Cannot perform symmetric matrix-matrix multiplication because compiled without BLAS"); \
    }
    ADEPT_DEFINE_SYMM(double, dsymm_, zsymm_);
    ADEPT_DEFINE_SYMM(float,  ssymm_, csymm_);
#undef ADEPT_DEFINE_SYMM
    
    // Matrix-vector multiplication where the matrix is symmetric
#define ADEPT_DEFINE_SYMV(T, FUNC, FUNC_COMPLEX)			\
    void cppblas_symv(const BLAS_ORDER Order,				\
		      const BLAS_UPLO Uplo,				\
		      const int N, const T alpha, const T *A,		\
		      const int lda, const T *X, const int incX,	\
		      const T beta, T *Y, const int incY) {		\
      throw feature_not_available("Cannot perform symmetric matrix-vector multiplication because compiled without BLAS"); \
    }
    ADEPT_DEFINE_SYMV(double, dsymv_, zsymv_);
    ADEPT_DEFINE_SYMV(float,  ssymv_, csymv_);
#undef ADEPT_DEFINE_SYMV
    
    // Matrix-vector multiplication for a general band matrix
#define ADEPT_DEFINE_GBMV(T, FUNC, FUNC_COMPLEX)		\
    void cppblas_gbmv(const BLAS_ORDER Order,			\
		      const BLAS_TRANSPOSE TransA,		\
		      const int M, const int N,			\
		      const int KL, const int KU, const T alpha,\
		      const T *A, const int lda, const T *X,	\
		      const int incX, const T beta, T *Y,	\
		      const int incY) {				\
      throw feature_not_available("Cannot perform band matrix-vector multiplication because compiled without BLAS"); \
    }
    ADEPT_DEFINE_GBMV(double, dgbmv_, zgbmv_);
    ADEPT_DEFINE_GBMV(float,  sgbmv_, cgbmv_);
#undef ADEPT_DEFINE_GBMV

  }
}

#endif


// =================================================================
// Contents of index.cpp
// =================================================================

/* index.cpp -- Definitions of "end" and "__" for array indexing

    Copyright (C) 2015 European Centre for Medium-Range Weather Forecasts

    Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.
*/

#include <adept/RangeIndex.h>

namespace adept {

  ::adept::internal::EndIndex end;
  ::adept::internal::AllIndex __;

}


// =================================================================
// Contents of inv.cpp
// =================================================================

/* inv.cpp -- Invert matrices

    Copyright (C) 2015-2016 European Centre for Medium-Range Weather Forecasts

    Author: Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.
*/
                             
#include <vector>

#include <adept/Array.h>
#include <adept/SpecialMatrix.h>

#ifndef AdeptSource_H
#include "cpplapack.h"
#endif

#ifdef HAVE_LAPACK

namespace adept {

  // -------------------------------------------------------------------
  // Invert general square matrix A
  // -------------------------------------------------------------------
  template <typename Type>
  Array<2,Type,false> 
  inv(const Array<2,Type,false>& A) {
    using internal::cpplapack_getrf;
    using internal::cpplapack_getri;

    if (A.dimension(0) != A.dimension(1)) {
      throw invalid_operation("Only square matrices can be inverted"
			      ADEPT_EXCEPTION_LOCATION);
    }

    Array<2,Type,false> A_;

    // LAPACKE is more efficient with column-major input
    A_.resize_column_major(A.dimensions());
    A_ = A;

    std::vector<lapack_int> ipiv(A_.dimension(0));

    //    lapack_int status = LAPACKE_dgetrf(LAPACK_COL_MAJOR, A_.dimension(0), A_.dimension(1),
    //				       A_.data(), A_.offset(1), &ipiv[0]);

    lapack_int status = cpplapack_getrf(A_.dimension(0),
					A_.data(), A_.offset(1), &ipiv[0]);
    if (status != 0) {
      std::stringstream s;
      s << "Failed to factorize matrix: LAPACK ?getrf returned code " << status;
      throw(matrix_ill_conditioned(s.str() ADEPT_EXCEPTION_LOCATION));
    }

    //    status = LAPACKE_dgetri(LAPACK_COL_MAJOR, A_.dimension(0),
    //			    A_.data(), A_.offset(1), &ipiv[0]);
    status = cpplapack_getri(A_.dimension(0),
			     A_.data(), A_.offset(1), &ipiv[0]);

    if (status != 0) {
      std::stringstream s;
      s << "Failed to invert matrix: LAPACK ?getri returned code " << status;
      throw(matrix_ill_conditioned(s.str() ADEPT_EXCEPTION_LOCATION));
    }
    return A_;
  }



  // -------------------------------------------------------------------
  // Invert symmetric matrix A
  // -------------------------------------------------------------------
  template <typename Type, SymmMatrixOrientation Orient>
  SpecialMatrix<Type,SymmEngine<Orient>,false> 
  inv(const SpecialMatrix<Type,SymmEngine<Orient>,false>& A) {
    using internal::cpplapack_sytrf;
    using internal::cpplapack_sytri;

    SpecialMatrix<Type,SymmEngine<Orient>,false> A_;

    A_.resize(A.dimension());
    A_ = A;

    // Treat symmetric matrix as column-major
    char uplo;
    if (Orient == ROW_LOWER_COL_UPPER) {
      uplo = 'U';
    }
    else {
      uplo = 'L';
    }

    std::vector<lapack_int> ipiv(A_.dimension(0));

    //    lapack_int status = LAPACKE_dsytrf(LAPACK_COL_MAJOR, uplo, A_.dimension(),
    //				       A_.data(), A_.offset(), &ipiv[0]);
    lapack_int status = cpplapack_sytrf(uplo, A_.dimension(),
					A_.data(), A_.offset(), &ipiv[0]);
    if (status != 0) {
      std::stringstream s;
      s << "Failed to factorize symmetric matrix: LAPACK ?sytrf returned code " << status;
      throw(matrix_ill_conditioned(s.str() ADEPT_EXCEPTION_LOCATION));
    }

    //    status = LAPACKE_dsytri(LAPACK_COL_MAJOR, uplo, A_.dimension(),
    //			    A_.data(), A_.offset(), &ipiv[0]);
    status = cpplapack_sytri(uplo, A_.dimension(),
			     A_.data(), A_.offset(), &ipiv[0]);
    if (status != 0) {
      std::stringstream s;
      s << "Failed to invert symmetric matrix: LAPACK ?sytri returned code " << status;
      throw(matrix_ill_conditioned(s.str() ADEPT_EXCEPTION_LOCATION));
    }
    return A_;
  }

}

#else // LAPACK not available
    
namespace adept {

  // -------------------------------------------------------------------
  // Invert general square matrix A
  // -------------------------------------------------------------------
  template <typename Type>
  Array<2,Type,false> 
  inv(const Array<2,Type,false>& A) {
    throw feature_not_available("Cannot invert matrix because compiled without LAPACK");
  }

  // -------------------------------------------------------------------
  // Invert symmetric matrix A
  // -------------------------------------------------------------------
  template <typename Type, SymmMatrixOrientation Orient>
  SpecialMatrix<Type,SymmEngine<Orient>,false> 
  inv(const SpecialMatrix<Type,SymmEngine<Orient>,false>& A) {
    throw feature_not_available("Cannot invert matrix because compiled without LAPACK");
  }
  
}

#endif

namespace adept {
  // -------------------------------------------------------------------
  // Explicit instantiations
  // -------------------------------------------------------------------
#define ADEPT_EXPLICIT_INV(TYPE)					\
  template Array<2,TYPE,false>						\
  inv(const Array<2,TYPE,false>& A);					\
  template SpecialMatrix<TYPE,SymmEngine<ROW_LOWER_COL_UPPER>,false>	\
  inv(const SpecialMatrix<TYPE,SymmEngine<ROW_LOWER_COL_UPPER>,false>&); \
  template SpecialMatrix<TYPE,SymmEngine<ROW_UPPER_COL_LOWER>,false>	\
  inv(const SpecialMatrix<TYPE,SymmEngine<ROW_UPPER_COL_LOWER>,false>&)

  ADEPT_EXPLICIT_INV(float);
  ADEPT_EXPLICIT_INV(double);

#undef ADEPT_EXPLICIT_INV
  
}




// =================================================================
// Contents of jacobian.cpp
// =================================================================

/* jacobian.cpp -- Computation of Jacobian matrix

    Copyright (C) 2012-2014 University of Reading
    Copyright (C) 2015-2016 European Centre for Medium-Range Weather Forecasts

    Author: Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.

*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include "adept/Stack.h"
#include "adept/Packet.h"
#include "adept/traits.h"

namespace adept {

  namespace internal {
    static const int MULTIPASS_SIZE = ADEPT_REAL_PACKET_SIZE == 1 ? ADEPT_MULTIPASS_SIZE : ADEPT_REAL_PACKET_SIZE;
  }

  using namespace internal;

  template <typename T>
  T _check_long_double() {
    // The user may have requested Real to be of type "long double" by
    // specifying ADEPT_REAL_TYPE_SIZE=16. If the present system can
    // only support double then sizeof(long double) will be 8, but
    // Adept will not be emitting the best code for this, so it is
    // probably better to fail forcing the user to specify
    // ADEPT_REAL_TYPE_SIZE=8.
    ADEPT_STATIC_ASSERT(ADEPT_REAL_TYPE_SIZE != 16 || ADEPT_REAL_TYPE_SIZE == sizeof(Real),
			COMPILER_DOES_NOT_SUPPORT_16_BYTE_LONG_DOUBLE);
    return 1;
  }

  /*
  void
  Stack::jacobian_forward_kernel(Real* gradient_multipass_b) const
  {
    static const int MULTIPASS_SIZE = Packet<Real>::size;

    // Loop forward through the derivative statements
    for (uIndex ist = 1; ist < n_statements_; ist++) {
      const Statement& statement = statement_[ist];
      // We copy the LHS to "a" in case it appears on the RHS in any
      // of the following statements
      Block<MULTIPASS_SIZE,Real> a; // Initialized to zero automatically
      
      // Loop through operations
      for (uIndex iop = statement_[ist-1].end_plus_one;
	   iop < statement.end_plus_one; iop++) {
	Real* __restrict grad = gradient_multipass_b+index_[iop]*MULTIPASS_SIZE;
	// Loop through columns within this block; we hope the
	// compiler can optimize this loop. Note that it is faster
	// to always use MULTIPASS_SIZE, always known at
	// compile time, than to use block_size, which is not, even
	// though in the last iteration this may involve redundant
	// computations.
	if (multiplier_[iop] == 1.0) {
	  //	    if (__builtin_expect(multiplier_[iop] == 1.0,0)) {
	  for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	    //	      for (uIndex i = 0; i < block_size; i++) {
	    a[i] += grad[i];
	  }
	}
	else {
	  for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	    //	      for (uIndex i = 0; i < block_size; i++) {
	    a[i] += multiplier_[iop]*grad[i];
	  }
	}
      }
      // Copy the results
      for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	gradient_multipass_b[statement.index*MULTIPASS_SIZE+i] = a[i];
      }
    } // End of loop over statements
  }    
  */

#if ADEPT_REAL_PACKET_SIZE > 1
  void
  Stack::jacobian_forward_kernel(Real* __restrict gradient_multipass_b) const
  {

    // Loop forward through the derivative statements
    for (uIndex ist = 1; ist < n_statements_; ist++) {
      const Statement& statement = statement_[ist];
      // We copy the LHS to "a" in case it appears on the RHS in any
      // of the following statements
      Packet<Real> a; // Zeroed automatically
      // Loop through operations
      for (uIndex iop = statement_[ist-1].end_plus_one;
	   iop < statement.end_plus_one; iop++) {
	Packet<Real> g(gradient_multipass_b+index_[iop]*MULTIPASS_SIZE);
	Packet<Real> m(multiplier_[iop]);
	a += m * g;
      }
      // Copy the results
      a.put(gradient_multipass_b+statement.index*MULTIPASS_SIZE);
    } // End of loop over statements
  }    
#else
  void
  Stack::jacobian_forward_kernel(Real* __restrict gradient_multipass_b) const
  {

    // Loop forward through the derivative statements
    for (uIndex ist = 1; ist < n_statements_; ist++) {
      const Statement& statement = statement_[ist];
      // We copy the LHS to "a" in case it appears on the RHS in any
      // of the following statements
      Block<MULTIPASS_SIZE,Real> a; // Zeroed automatically
      // Loop through operations
      for (uIndex iop = statement_[ist-1].end_plus_one;
	   iop < statement.end_plus_one; iop++) {
	for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	  a[i] += multiplier_[iop]*gradient_multipass_b[index_[iop]*MULTIPASS_SIZE+i];
	}
      }
      // Copy the results
      for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	gradient_multipass_b[statement.index*MULTIPASS_SIZE+i] = a[i];
      }
    } // End of loop over statements
  }    
#endif

  void
  Stack::jacobian_forward_kernel_extra(Real* __restrict gradient_multipass_b,
				       uIndex n_extra) const
  {

    // Loop forward through the derivative statements
    for (uIndex ist = 1; ist < n_statements_; ist++) {
      const Statement& statement = statement_[ist];
      // We copy the LHS to "a" in case it appears on the RHS in any
      // of the following statements
      Block<MULTIPASS_SIZE,Real> a; // Zeroed automatically
      // Loop through operations
      for (uIndex iop = statement_[ist-1].end_plus_one;
	   iop < statement.end_plus_one; iop++) {
	for (uIndex i = 0; i < n_extra; i++) {
	  a[i] += multiplier_[iop]*gradient_multipass_b[index_[iop]*MULTIPASS_SIZE+i];
	}
      }
      // Copy the results
      for (uIndex i = 0; i < n_extra; i++) {
	gradient_multipass_b[statement.index*MULTIPASS_SIZE+i] = a[i];
      }
    } // End of loop over statements
  }    



  // Compute the Jacobian matrix, parallelized using OpenMP. Normally
  // the user would call the jacobian or jacobian_forward functions,
  // and the OpenMP version would only be called if OpenMP is
  // available and the Jacobian matrix is large enough for
  // parallelization to be worthwhile.  Note that jacobian_out must be
  // allocated to be of size m*n, where m is the number of dependent
  // variables and n is the number of independents. The independents
  // and dependents must have already been identified with the
  // functions "independent" and "dependent", otherwise this function
  // will fail with FAILURE_XXDEPENDENT_NOT_IDENTIFIED. In the
  // resulting matrix, the "m" dimension of the matrix varies
  // fastest. This is implemented using a forward pass, appropriate
  // for m>=n.
  void
  Stack::jacobian_forward_openmp(Real* jacobian_out) const
  {

    // Number of blocks to cycle through, including a possible last
    // block containing fewer than MULTIPASS_SIZE variables
    int n_block = (n_independent() + MULTIPASS_SIZE - 1)
      / MULTIPASS_SIZE;
    uIndex n_extra = n_independent() % MULTIPASS_SIZE;
    
    int iblock;
    
#pragma omp parallel
    {
      //      std::vector<Block<MULTIPASS_SIZE,Real> > 
      //	gradient_multipass_b(max_gradient_);
      uIndex gradient_multipass_size = max_gradient_*MULTIPASS_SIZE;
      Real* __restrict gradient_multipass_b 
	= alloc_aligned<Real>(gradient_multipass_size);
      
#pragma omp for schedule(static)
      for (iblock = 0; iblock < n_block; iblock++) {
	// Set the index to the dependent variables for this block
	uIndex i_independent =  MULTIPASS_SIZE * iblock;
	
	uIndex block_size = MULTIPASS_SIZE;
	// If this is the last iteration and the number of extra
	// elements is non-zero, then set the block size to the number
	// of extra elements. If the number of extra elements is zero,
	// then the number of independent variables is exactly divisible
	// by MULTIPASS_SIZE, so the last iteration will be the
	// same as all the rest.
	if (iblock == n_block-1 && n_extra > 0) {
	  block_size = n_extra;
	}
	
	// Set the initial gradients all to zero
	for (std::size_t i = 0; i < gradient_multipass_size; i++) {
	  gradient_multipass_b[i] = 0.0;
	}
	// Each seed vector has one non-zero entry of 1.0
	for (uIndex i = 0; i < block_size; i++) {
	  gradient_multipass_b[independent_index_[i_independent+i]*MULTIPASS_SIZE+i] = 1.0;
	}

	jacobian_forward_kernel(gradient_multipass_b);

	// Copy the gradients corresponding to the dependent variables
	// into the Jacobian matrix
	for (uIndex idep = 0; idep < n_dependent(); idep++) {
	  for (uIndex i = 0; i < block_size; i++) {
	    jacobian_out[(i_independent+i)*n_dependent()+idep]
	      = gradient_multipass_b[dependent_index_[idep]*MULTIPASS_SIZE+i];
	  }
	}
      } // End of loop over blocks
      free_aligned(gradient_multipass_b);
    } // End of parallel section
  } // End of jacobian function


  // Compute the Jacobian matrix; note that jacobian_out must be
  // allocated to be of size m*n, where m is the number of dependent
  // variables and n is the number of independents. The independents
  // and dependents must have already been identified with the
  // functions "independent" and "dependent", otherwise this function
  // will fail with FAILURE_XXDEPENDENT_NOT_IDENTIFIED. In the
  // resulting matrix, the "m" dimension of the matrix varies
  // fastest. This is implemented using a forward pass, appropriate
  // for m>=n.
  void
  Stack::jacobian_forward(Real* jacobian_out)
  {
    if (independent_index_.empty() || dependent_index_.empty()) {
      throw(dependents_or_independents_not_identified());
    }
#ifdef _OPENMP
    if (have_openmp_ 
	&& !openmp_manually_disabled_
	&& n_independent() > MULTIPASS_SIZE
	&& omp_get_max_threads() > 1) {
      // Call the parallel version
      jacobian_forward_openmp(jacobian_out);
      return;
    }
#endif

    // For optimization reasons, we process a block of
    // MULTIPASS_SIZE columns of the Jacobian at once; calculate
    // how many blocks are needed and how many extras will remain
    uIndex n_block = n_independent() / MULTIPASS_SIZE;
    uIndex n_extra = n_independent() % MULTIPASS_SIZE;

    ///gradient_multipass_.resize(max_gradient_);
    uIndex gradient_multipass_size = max_gradient_*MULTIPASS_SIZE;
    Real* __restrict gradient_multipass_b 
      = alloc_aligned<Real>(gradient_multipass_size);

    // Loop over blocks of MULTIPASS_SIZE columns
    for (uIndex iblock = 0; iblock < n_block; iblock++) {
      // Set the index to the dependent variables for this block
      uIndex i_independent =  MULTIPASS_SIZE * iblock;

      // Set the initial gradients all to zero
      ///zero_gradient_multipass();
      for (std::size_t i = 0; i < gradient_multipass_size; i++) {
	gradient_multipass_b[i] = 0.0;
      }

      // Each seed vector has one non-zero entry of 1.0
      for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	gradient_multipass_b[independent_index_[i_independent+i]*MULTIPASS_SIZE+i] = 1.0;
      }

      jacobian_forward_kernel(gradient_multipass_b);

      // Copy the gradients corresponding to the dependent variables
      // into the Jacobian matrix
      for (uIndex idep = 0; idep < n_dependent(); idep++) {
	for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	  jacobian_out[(i_independent+i)*n_dependent()+idep] 
	    = gradient_multipass_b[dependent_index_[idep]*MULTIPASS_SIZE+i];
	}
      }
      i_independent += MULTIPASS_SIZE;
    } // End of loop over blocks
    
    // Now do the same but for the remaining few columns in the matrix
    if (n_extra > 0) {
      uIndex i_independent =  MULTIPASS_SIZE * n_block;
      ///zero_gradient_multipass();
      for (std::size_t i = 0; i < gradient_multipass_size; i++) {
	gradient_multipass_b[i] = 0.0;
      }

      for (uIndex i = 0; i < n_extra; i++) {
	gradient_multipass_b[independent_index_[i_independent+i]*MULTIPASS_SIZE+i] = 1.0;
      }

      jacobian_forward_kernel_extra(gradient_multipass_b, n_extra);

      for (uIndex idep = 0; idep < n_dependent(); idep++) {
	for (uIndex i = 0; i < n_extra; i++) {
	  jacobian_out[(i_independent+i)*n_dependent()+idep] 
	    = gradient_multipass_b[dependent_index_[idep]*MULTIPASS_SIZE+i];
	}
      }
    }

    free_aligned(gradient_multipass_b);
  }


  // Compute the Jacobian matrix, parallelized using OpenMP.  Normally
  // the user would call the jacobian or jacobian_reverse functions,
  // and the OpenMP version would only be called if OpenMP is
  // available and the Jacobian matrix is large enough for
  // parallelization to be worthwhile.  Note that jacobian_out must be
  // allocated to be of size m*n, where m is the number of dependent
  // variables and n is the number of independents. The independents
  // and dependents must have already been identified with the
  // functions "independent" and "dependent", otherwise this function
  // will fail with FAILURE_XXDEPENDENT_NOT_IDENTIFIED. In the
  // resulting matrix, the "m" dimension of the matrix varies
  // fastest. This is implemented using a reverse pass, appropriate
  // for m<n.
  void
  Stack::jacobian_reverse_openmp(Real* jacobian_out) const
  {

    // Number of blocks to cycle through, including a possible last
    // block containing fewer than MULTIPASS_SIZE variables
    int n_block = (n_dependent() + MULTIPASS_SIZE - 1)
      / MULTIPASS_SIZE;
    uIndex n_extra = n_dependent() % MULTIPASS_SIZE;
    
    int iblock;

    // Inside the OpenMP loop, the "this" pointer may be NULL if the
    // adept::Stack pointer is declared as thread-local and if the
    // OpenMP memory model uses thread-local storage for private
    // data. If this is the case then local pointers to or copies of
    // the following members of the adept::Stack object may need to be
    // made: dependent_index_ n_statements_ statement_ multiplier_
    // index_ independent_index_ n_dependent() n_independent().
    // Limited testing implies this is OK though.

#pragma omp parallel
    {
      std::vector<Block<MULTIPASS_SIZE,Real> > 
	gradient_multipass_b(max_gradient_);
      
#pragma omp for schedule(static)
      for (iblock = 0; iblock < n_block; iblock++) {
	// Set the index to the dependent variables for this block
	uIndex i_dependent =  MULTIPASS_SIZE * iblock;
	
	uIndex block_size = MULTIPASS_SIZE;
	// If this is the last iteration and the number of extra
	// elements is non-zero, then set the block size to the number
	// of extra elements. If the number of extra elements is zero,
	// then the number of independent variables is exactly divisible
	// by MULTIPASS_SIZE, so the last iteration will be the
	// same as all the rest.
	if (iblock == n_block-1 && n_extra > 0) {
	  block_size = n_extra;
	}

	// Set the initial gradients all to zero
	for (std::size_t i = 0; i < gradient_multipass_b.size(); i++) {
	  gradient_multipass_b[i].zero();
	}
	// Each seed vector has one non-zero entry of 1.0
	for (uIndex i = 0; i < block_size; i++) {
	  gradient_multipass_b[dependent_index_[i_dependent+i]][i] = 1.0;
	}

	// Loop backward through the derivative statements
	for (uIndex ist = n_statements_-1; ist > 0; ist--) {
	  const Statement& statement = statement_[ist];
	  // We copy the RHS to "a" in case it appears on the LHS in any
	  // of the following statements
	  Real a[MULTIPASS_SIZE];
#if MULTIPASS_SIZE > MULTIPASS_SIZE_ZERO_CHECK
	  // For large blocks, we only process the ones where a[i] is
	  // non-zero
	  uIndex i_non_zero[MULTIPASS_SIZE];
#endif
	  uIndex n_non_zero = 0;
	  for (uIndex i = 0; i < block_size; i++) {
	    a[i] = gradient_multipass_b[statement.index][i];
	    gradient_multipass_b[statement.index][i] = 0.0;
	    if (a[i] != 0.0) {
#if MULTIPASS_SIZE > MULTIPASS_SIZE_ZERO_CHECK
	      i_non_zero[n_non_zero++] = i;
#else
	      n_non_zero = 1;
#endif
	    }
	  }

	  // Only do anything for this statement if any of the a values
	  // are non-zero
	  if (n_non_zero) {
	    // Loop through the operations
	    for (uIndex iop = statement_[ist-1].end_plus_one;
		 iop < statement.end_plus_one; iop++) {
	      // Try to minimize pointer dereferencing by making local
	      // copies
	      Real multiplier = multiplier_[iop];
	      Real* __restrict gradient_multipass 
		= &(gradient_multipass_b[index_[iop]][0]);
#if MULTIPASS_SIZE > MULTIPASS_SIZE_ZERO_CHECK
	      // For large blocks, loop over only the indices
	      // corresponding to non-zero a
	      for (uIndex i = 0; i < n_non_zero; i++) {
		gradient_multipass[i_non_zero[i]] += multiplier*a[i_non_zero[i]];
	      }
#else
	      // For small blocks, do all indices
	      for (uIndex i = 0; i < block_size; i++) {
	      //	      for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
		gradient_multipass[i] += multiplier*a[i];
	      }
#endif
	    }
	  }
	} // End of loop over statement
	// Copy the gradients corresponding to the independent
	// variables into the Jacobian matrix
	for (uIndex iindep = 0; iindep < n_independent(); iindep++) {
	  for (uIndex i = 0; i < block_size; i++) {
	    jacobian_out[iindep*n_dependent()+i_dependent+i] 
	      = gradient_multipass_b[independent_index_[iindep]][i];
	  }
	}
      } // End of loop over blocks
    } // end #pragma omp parallel
  } // end jacobian_reverse_openmp


  // Compute the Jacobian matrix; note that jacobian_out must be
  // allocated to be of size m*n, where m is the number of dependent
  // variables and n is the number of independents. The independents
  // and dependents must have already been identified with the
  // functions "independent" and "dependent", otherwise this function
  // will fail with FAILURE_XXDEPENDENT_NOT_IDENTIFIED. In the
  // resulting matrix, the "m" dimension of the matrix varies
  // fastest. This is implemented using a reverse pass, appropriate
  // for m<n.
  void
  Stack::jacobian_reverse(Real* jacobian_out)
  {
    if (independent_index_.empty() || dependent_index_.empty()) {
      throw(dependents_or_independents_not_identified());
    }
#ifdef _OPENMP
    if (have_openmp_ 
	&& !openmp_manually_disabled_
	&& n_dependent() > MULTIPASS_SIZE
	&& omp_get_max_threads() > 1) {
      // Call the parallel version
      jacobian_reverse_openmp(jacobian_out);
      return;
    }
#endif

    //    gradient_multipass_.resize(max_gradient_);
    std::vector<Block<MULTIPASS_SIZE,Real> > 
      gradient_multipass_b(max_gradient_);

    // For optimization reasons, we process a block of
    // MULTIPASS_SIZE rows of the Jacobian at once; calculate
    // how many blocks are needed and how many extras will remain
    uIndex n_block = n_dependent() / MULTIPASS_SIZE;
    uIndex n_extra = n_dependent() % MULTIPASS_SIZE;
    uIndex i_dependent = 0; // uIndex of first row in the block we are
			    // currently computing
    // Loop over the of MULTIPASS_SIZE rows
    for (uIndex iblock = 0; iblock < n_block; iblock++) {
      // Set the initial gradients all to zero
      //      zero_gradient_multipass();
      for (std::size_t i = 0; i < gradient_multipass_b.size(); i++) {
	gradient_multipass_b[i].zero();
      }

      // Each seed vector has one non-zero entry of 1.0
      for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	gradient_multipass_b[dependent_index_[i_dependent+i]][i] = 1.0;
      }
      // Loop backward through the derivative statements
      for (uIndex ist = n_statements_-1; ist > 0; ist--) {
	const Statement& statement = statement_[ist];
	// We copy the RHS to "a" in case it appears on the LHS in any
	// of the following statements
	Real a[MULTIPASS_SIZE];
#if MULTIPASS_SIZE > MULTIPASS_SIZE_ZERO_CHECK
	// For large blocks, we only process the ones where a[i] is
	// non-zero
	uIndex i_non_zero[MULTIPASS_SIZE];
#endif
	uIndex n_non_zero = 0;
	for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	  a[i] = gradient_multipass_b[statement.index][i];
	  gradient_multipass_b[statement.index][i] = 0.0;
	  if (a[i] != 0.0) {
#if MULTIPASS_SIZE > MULTIPASS_SIZE_ZERO_CHECK
	    i_non_zero[n_non_zero++] = i;
#else
	    n_non_zero = 1;
#endif
	  }
	}
	// Only do anything for this statement if any of the a values
	// are non-zero
	if (n_non_zero) {
	  // Loop through the operations
	  for (uIndex iop = statement_[ist-1].end_plus_one;
	       iop < statement.end_plus_one; iop++) {
	    // Try to minimize pointer dereferencing by making local
	    // copies
	    Real multiplier = multiplier_[iop];
	    Real* __restrict gradient_multipass 
	      = &(gradient_multipass_b[index_[iop]][0]);
#if MULTIPASS_SIZE > MULTIPASS_SIZE_ZERO_CHECK
	    // For large blocks, loop over only the indices
	    // corresponding to non-zero a
	    for (uIndex i = 0; i < n_non_zero; i++) {
	      gradient_multipass[i_non_zero[i]] += multiplier*a[i_non_zero[i]];
	    }
#else
	    // For small blocks, do all indices
	    for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	      gradient_multipass[i] += multiplier*a[i];
	    }
#endif
	  }
	}
      } // End of loop over statement
      // Copy the gradients corresponding to the independent variables
      // into the Jacobian matrix
      for (uIndex iindep = 0; iindep < n_independent(); iindep++) {
	for (uIndex i = 0; i < MULTIPASS_SIZE; i++) {
	  jacobian_out[iindep*n_dependent()+i_dependent+i] 
	    = gradient_multipass_b[independent_index_[iindep]][i];
	}
      }
      i_dependent += MULTIPASS_SIZE;
    } // End of loop over blocks
    
    // Now do the same but for the remaining few rows in the matrix
    if (n_extra > 0) {
      for (std::size_t i = 0; i < gradient_multipass_b.size(); i++) {
	gradient_multipass_b[i].zero();
      }
      //      zero_gradient_multipass();
      for (uIndex i = 0; i < n_extra; i++) {
	gradient_multipass_b[dependent_index_[i_dependent+i]][i] = 1.0;
      }
      for (uIndex ist = n_statements_-1; ist > 0; ist--) {
	const Statement& statement = statement_[ist];
	Real a[MULTIPASS_SIZE];
#if MULTIPASS_SIZE > MULTIPASS_SIZE_ZERO_CHECK
	uIndex i_non_zero[MULTIPASS_SIZE];
#endif
	uIndex n_non_zero = 0;
	for (uIndex i = 0; i < n_extra; i++) {
	  a[i] = gradient_multipass_b[statement.index][i];
	  gradient_multipass_b[statement.index][i] = 0.0;
	  if (a[i] != 0.0) {
#if MULTIPASS_SIZE > MULTIPASS_SIZE_ZERO_CHECK
	    i_non_zero[n_non_zero++] = i;
#else
	    n_non_zero = 1;
#endif
	  }
	}
	if (n_non_zero) {
	  for (uIndex iop = statement_[ist-1].end_plus_one;
	       iop < statement.end_plus_one; iop++) {
	    Real multiplier = multiplier_[iop];
	    Real* __restrict gradient_multipass 
	      = &(gradient_multipass_b[index_[iop]][0]);
	    //	    if (index_[iop] > max_gradient_-1
	    //		|| index_[iop] < 0) {
	    //	    std::cerr << "AAAAAA: iop=" << iop << " index_[iop]=" << index_[iop] << " max_gradient_=" << max_gradient_ << " ist=" << ist << "\n";
	      //	    }
#if MULTIPASS_SIZE > MULTIPASS_SIZE_ZERO_CHECK
	    for (uIndex i = 0; i < n_non_zero; i++) {
	      gradient_multipass[i_non_zero[i]] += multiplier*a[i_non_zero[i]];
	    }
#else
	    for (uIndex i = 0; i < n_extra; i++) {
	      //	      std::cerr << "BBBBB: i=" << i << " gradient_multipass[i]=" << gradient_multipass[i] << " multiplier=" << multiplier << " a[i]=" << a[i] << "\n";
	      gradient_multipass[i] += multiplier*a[i];
	    }
#endif
	  }
	}
      }
      for (uIndex iindep = 0; iindep < n_independent(); iindep++) {
	for (uIndex i = 0; i < n_extra; i++) {
	  jacobian_out[iindep*n_dependent()+i_dependent+i] 
	    = gradient_multipass_b[independent_index_[iindep]][i];
	}
      }
    }
  }

  // Compute the Jacobian matrix; note that jacobian_out must be
  // allocated to be of size m*n, where m is the number of dependent
  // variables and n is the number of independents. In the resulting
  // matrix, the "m" dimension of the matrix varies fastest. This is
  // implemented by calling one of jacobian_forward and
  // jacobian_reverse, whichever would be faster.
  void
  Stack::jacobian(Real* jacobian_out)
  {
    //    std::cout << ">>> Computing " << n_dependent() << "x" << n_independent()
    //	      << " Jacobian from " << n_statements_ << " statements, "
    //	      << n_operations() << " operations and " << max_gradient_ << " gradients\n";

    if (n_independent() <= n_dependent()) {
      jacobian_forward(jacobian_out);
    }
    else {
      jacobian_reverse(jacobian_out);
    }
  }
  
} // End namespace adept


// =================================================================
// Contents of settings.cpp
// =================================================================

/* settings.cpp -- View/change the overall Adept settings

    Copyright (C) 2016 European Centre for Medium-Range Weather Forecasts

    Author: Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.

*/

#include <sstream>
#include <cstring>

#include <adept/base.h>
#include <adept/settings.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_OPENBLAS_CBLAS_HEADER
#include <cblas.h>
#endif

namespace adept {

  // -------------------------------------------------------------------
  // Get compile-time settings
  // -------------------------------------------------------------------

  // Return the version of Adept at compile time
  std::string
  version()
  {
    return ADEPT_VERSION_STR;
  }

  // Return the compiler used to compile the Adept library (e.g. "g++
  // [4.3.2]" or "Microsoft Visual C++ [1800]")
  std::string
  compiler_version()
  {
#ifdef CXX
    std::string cv = CXX; // Defined in config.h
#elif defined(_MSC_VER)
    std::string cv = "Microsoft Visual C++";
#else
    std::string cv = "unknown";
#endif

#ifdef __GNUC__

#define STRINGIFY3(A,B,C) STRINGIFY(A) "." STRINGIFY(B) "." STRINGIFY(C)
#define STRINGIFY(A) #A
    cv += " [" STRINGIFY3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__) "]";
#undef STRINGIFY
#undef STRINGIFY3

#elif defined(_MSC_VER)

#define STRINGIFY1(A) STRINGIFY(A)
#define STRINGIFY(A) #A
    cv += " [" STRINGIFY1(_MSC_VER) "]";
#undef STRINGIFY
#undef STRINGIFY1

#endif
    return cv;
  }

  // Return the compiler flags used when compiling the Adept library
  // (e.g. "-Wall -g -O3")
  std::string
  compiler_flags()
  {
#ifdef CXXFLAGS
    return CXXFLAGS; // Defined in config.h
#else
    return "unknown";
#endif
  }

  // Return a multi-line string listing numerous aspects of the way
  // Adept has been configured.
  std::string
  configuration()
  {
    std::stringstream s;
    s << "Adept version " << adept::version() << ":\n";
    s << "  Compiled with " << adept::compiler_version() << "\n";
    s << "  Compiler flags \"" << adept::compiler_flags() << "\"\n";
#ifdef BLAS_LIBS
    if (std::strlen(BLAS_LIBS) > 2) {
      const char* blas_libs = BLAS_LIBS + 2;
      s << "  BLAS support from " << blas_libs << " library\n";
    }
    else {
      s << "  BLAS support from built-in library\n";
    }
#endif
#ifdef HAVE_OPENBLAS_CBLAS_HEADER
    s << "  Number of BLAS threads may be specified up to maximum of "
      << max_blas_threads() << "\n";
#endif
    s << "  Jacobians processed in blocks of size " 
      << ADEPT_MULTIPASS_SIZE << "\n";
    return s.str();
  }


  // -------------------------------------------------------------------
  // Get/set number of threads for array operations
  // -------------------------------------------------------------------

  // Get the maximum number of threads available for BLAS operations
  int
  max_blas_threads()
  {
#ifdef HAVE_OPENBLAS_CBLAS_HEADER
    return openblas_get_num_threads();
#else
    return 1;
#endif
  }

  // Set the maximum number of threads available for BLAS operations
  // (zero means use the maximum sensible number on the current
  // system), and return the number actually set. Note that OpenBLAS
  // uses pthreads and the Jacobian calculation uses OpenMP - this can
  // lead to inefficient behaviour so if you are computing Jacobians
  // then you may get better performance by setting the number of
  // array threads to one.
  int
  set_max_blas_threads(int n)
  {
#ifdef HAVE_OPENBLAS_CBLAS_HEADER
    openblas_set_num_threads(n);
    return openblas_get_num_threads();
#else
    return 1;
#endif
  }

  // Was the library compiled with matrix multiplication support (from
  // BLAS)?
  bool
  have_matrix_multiplication() {
#ifdef HAVE_BLAS
    return true;
#else
    return false;
#endif
  }

  // Was the library compiled with linear algebra support (e.g. inv
  // and solve from LAPACK)
  bool
  have_linear_algebra() {
#ifdef HAVE_LAPACK
    return true;
#else
    return false;
#endif
  }

} // End namespace adept


// =================================================================
// Contents of solve.cpp
// =================================================================

/* solve.cpp -- Solve systems of linear equations using LAPACK

    Copyright (C) 2015-2016 European Centre for Medium-Range Weather Forecasts

    Author: Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.
*/
                             

#include <vector>


#include <adept/solve.h>
#include <adept/Array.h>
#include <adept/SpecialMatrix.h>

// If ADEPT_SOURCE_H is defined then we are in a header file generated
// from all the source files, so cpplapack.h will already have been
// included
#ifndef AdeptSource_H
#include "cpplapack.h"
#endif

#ifdef HAVE_LAPACK

namespace adept {

  // -------------------------------------------------------------------
  // Solve Ax = b for general square matrix A
  // -------------------------------------------------------------------
  template <typename T>
  Array<1,T,false> 
  solve(const Array<2,T,false>& A, const Array<1,T,false>& b) {
    Array<2,T,false> A_;
    Array<1,T,false> b_;

    // LAPACKE is more efficient with column-major input
    // if (A.is_row_contiguous()) {
      A_.resize_column_major(A.dimensions());
      A_ = A;
    // }
    // else {
    //   A_.link(A);
    // }

    // if (b_.offset(0) != 0) {
      b_ = b;
    // }
    // else {
    //   b_.link(b);
    // }

    std::vector<lapack_int> ipiv(A_.dimension(0));

    //    lapack_int status = LAPACKE_dgesv(LAPACK_COL_MAJOR, A_.dimension(0), 1,
    //				      A_.data(), A_.offset(1), &ipiv[0],
    //				      b_.data(), b_.dimension(0));
    lapack_int status = cpplapack_gesv(A_.dimension(0), 1,
				       A_.data(), A_.offset(1), &ipiv[0],
				       b_.data(), b_.dimension(0));

    if (status != 0) {
      std::stringstream s;
      s << "Failed to solve general system of equations: LAPACK ?gesv returned code " << status;
      throw(matrix_ill_conditioned(s.str() ADEPT_EXCEPTION_LOCATION));
    }
    return b_;    
  }

  // -------------------------------------------------------------------
  // Solve AX = B for general square matrix A and rectangular matrix B
  // -------------------------------------------------------------------
  template <typename T>
  Array<2,T,false> 
  solve(const Array<2,T,false>& A, const Array<2,T,false>& B) {
    Array<2,T,false> A_;
    Array<2,T,false> B_;
    
    // LAPACKE is more efficient with column-major input
    // if (A.is_row_contiguous()) {
      A_.resize_column_major(A.dimensions());
      A_ = A;
    // }
    // else {
    //   A_.link(A);
    // }

    // if (B.is_row_contiguous()) {
      B_.resize_column_major(B.dimensions());
      B_ = B;
    // }
    // else {
    //   B_.link(B);
    // }

    std::vector<lapack_int> ipiv(A_.dimension(0));

    //    lapack_int status = LAPACKE_dgesv(LAPACK_COL_MAJOR, A_.dimension(0), B.dimension(1),
    //				      A_.data(), A_.offset(1), &ipiv[0],
    //				      B_.data(), B_.offset(1));
    lapack_int status = cpplapack_gesv(A_.dimension(0), B.dimension(1),
				       A_.data(), A_.offset(1), &ipiv[0],
				       B_.data(), B_.offset(1));
    if (status != 0) {
      std::stringstream s;
      s << "Failed to solve general system of equations for matrix RHS: LAPACK ?gesv returned code " << status;
      throw(matrix_ill_conditioned(s.str() ADEPT_EXCEPTION_LOCATION));
    }
    return B_;    
  }


  // -------------------------------------------------------------------
  // Solve Ax = b for symmetric square matrix A
  // -------------------------------------------------------------------
  template <typename T, SymmMatrixOrientation Orient>
  Array<1,T,false>
  solve(const SpecialMatrix<T,SymmEngine<Orient>,false>& A,
	const Array<1,T,false>& b) {
    SpecialMatrix<T,SymmEngine<Orient>,false> A_;
    Array<1,T,false> b_;

    // Not sure why the original code copies A...
    A_.resize(A.dimension());
    A_ = A;
    // A_.link(A);

    // if (b.offset(0) != 1) {
      b_ = b;
    // }
    // else {
    //   b_.link(b);
    // }

    // Treat symmetric matrix as column-major
    char uplo;
    if (Orient == ROW_LOWER_COL_UPPER) {
      uplo = 'U';
    }
    else {
      uplo = 'L';
    }

    std::vector<lapack_int> ipiv(A_.dimension());

    //    lapack_int status = LAPACKE_dsysv(LAPACK_COL_MAJOR, uplo, A_.dimension(0), 1,
    //				      A_.data(), A_.offset(), &ipiv[0],
    //				      b_.data(), b_.dimension(0));
    lapack_int status = cpplapack_sysv(uplo, A_.dimension(0), 1,
				       A_.data(), A_.offset(), &ipiv[0],
				       b_.data(), b_.dimension(0));

    if (status != 0) {
      //      std::stringstream s;
      //      s << "Failed to solve symmetric system of equations: LAPACK ?sysv returned code " << status;
      //      throw(matrix_ill_conditioned(s.str() ADEPT_EXCEPTION_LOCATION));
      std::cerr << "Warning: LAPACK solve symmetric system failed (?sysv): trying general (?gesv)\n";
      return solve(Array<2,T,false>(A_),b_);
    }
    return b_;    
  }


  // -------------------------------------------------------------------
  // Solve AX = B for symmetric square matrix A
  // -------------------------------------------------------------------
  template <typename T, SymmMatrixOrientation Orient>
  Array<2,T,false>
  solve(const SpecialMatrix<T,SymmEngine<Orient>,false>& A,
	const Array<2,T,false>& B) {
    SpecialMatrix<T,SymmEngine<Orient>,false> A_;
    Array<2,T,false> B_;

    A_.resize(A.dimension());
    A_ = A;
    // A_.link(A);

    // if (B.is_row_contiguous()) {
      B_.resize_column_major(B.dimensions());
      B_ = B;
    // }
    // else {
    //   B_.link(B);
    // }

    // Treat symmetric matrix as column-major
    char uplo;
    if (Orient == ROW_LOWER_COL_UPPER) {
      uplo = 'U';
    }
    else {
      uplo = 'L';
    }

    std::vector<lapack_int> ipiv(A_.dimension());

    //    lapack_int status = LAPACKE_dsysv(LAPACK_COL_MAJOR, uplo, A_.dimension(0), B.dimension(1),
    //				      A_.data(), A_.offset(), &ipiv[0],
    //				      B_.data(), B_.offset(1));
    lapack_int status = cpplapack_sysv(uplo, A_.dimension(0), B.dimension(1),
				       A_.data(), A_.offset(), &ipiv[0],
				       B_.data(), B_.offset(1));

    if (status != 0) {
      std::stringstream s;
      s << "Failed to solve symmetric system of equations with matrix RHS: LAPACK ?sysv returned code " << status;
      throw(matrix_ill_conditioned(s.str() ADEPT_EXCEPTION_LOCATION));
    }
    return B_;
  }

}

#else

namespace adept {
  
  // -------------------------------------------------------------------
  // Solve Ax = b for general square matrix A
  // -------------------------------------------------------------------
  template <typename T>
  Array<1,T,false> 
  solve(const Array<2,T,false>& A, const Array<1,T,false>& b) {
    throw feature_not_available("Cannot solve linear equations because compiled without LAPACK");
  }

  // -------------------------------------------------------------------
  // Solve AX = B for general square matrix A and rectangular matrix B
  // -------------------------------------------------------------------
  template <typename T>
  Array<2,T,false> 
  solve(const Array<2,T,false>& A, const Array<2,T,false>& B) {
    throw feature_not_available("Cannot solve linear equations because compiled without LAPACK");
  }

  // -------------------------------------------------------------------
  // Solve Ax = b for symmetric square matrix A
  // -------------------------------------------------------------------
  template <typename T, SymmMatrixOrientation Orient>
  Array<1,T,false>
  solve(const SpecialMatrix<T,SymmEngine<Orient>,false>& A,
	const Array<1,T,false>& b) {
    throw feature_not_available("Cannot solve linear equations because compiled without LAPACK");
  }

  // -------------------------------------------------------------------
  // Solve AX = B for symmetric square matrix A
  // -------------------------------------------------------------------
  template <typename T, SymmMatrixOrientation Orient>
  Array<2,T,false>
  solve(const SpecialMatrix<T,SymmEngine<Orient>,false>& A,
	const Array<2,T,false>& B) {
    throw feature_not_available("Cannot solve linear equations because compiled without LAPACK");
  }

}

#endif


namespace adept {

  // -------------------------------------------------------------------
  // Explicit instantiations
  // -------------------------------------------------------------------
#define ADEPT_EXPLICIT_SOLVE(TYPE,RRANK)				\
  template Array<RRANK,TYPE,false>					\
  solve(const Array<2,TYPE,false>& A, const Array<RRANK,TYPE,false>& b); \
  template Array<RRANK,TYPE,false>					\
  solve(const SpecialMatrix<TYPE,SymmEngine<ROW_LOWER_COL_UPPER>,false>& A, \
	const Array<RRANK,TYPE,false>& b);					\
  template Array<RRANK,TYPE,false>					\
  solve(const SpecialMatrix<TYPE,SymmEngine<ROW_UPPER_COL_LOWER>,false>& A, \
	const Array<RRANK,TYPE,false>& b);

  ADEPT_EXPLICIT_SOLVE(float,1)
  ADEPT_EXPLICIT_SOLVE(float,2)
  ADEPT_EXPLICIT_SOLVE(double,1)
  ADEPT_EXPLICIT_SOLVE(double,2)
#undef ADEPT_EXPLICIT_SOLVE

}



// =================================================================
// Contents of vector_utilities.cpp
// =================================================================

/* vector_utilities.cpp -- Vector utility functions

    Copyright (C) 2016 European Centre for Medium-Range Weather Forecasts

    Author: Robin Hogan <r.j.hogan@ecmwf.int>

    This file is part of the Adept library.

*/

#include <adept/vector_utilities.h>

namespace adept {

  Array<1,Real,false>
  linspace(Real x1, Real x2, Index n) {
    Array<1,Real,false> ans(n);
    if (n > 1) {
      for (Index i = 0; i < n; ++i) {
	ans(i) = x1 + (x2-x1)*i / static_cast<Real>(n-1);
      }
    }
    else if (n == 1 && x1 == x2) {
      ans(0) = x1;
      return ans;
    }
    else if (n == 1) {
      throw(invalid_operation("linspace(x1,x2,n) with n=1 only valid if x1=x2"));
    }
    return ans;
  }

}



#endif

