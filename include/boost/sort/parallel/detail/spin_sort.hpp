//----------------------------------------------------------------------------
/// @file spin_sort.hpp
/// @brief Spin Sort algorithm
///
/// @author Copyright (c) 2016 Francisco José Tapia (fjtapia@gmail.com )\n
///         Distributed under the Boost Software License, Version 1.0.\n
///         ( See accompanying file LICENSE_1_0.txt or copy at
///           http://www.boost.org/LICENSE_1_0.txt  )
/// @version 0.1
///
/// @remarks
//-----------------------------------------------------------------------------
#ifndef __BOOST_SORT_PARALLEL_ALGORITHM_SPIN_SORT_HPP
#define __BOOST_SORT_PARALLEL_ALGORITHM_SPIN_SORT_HPP

#include <boost/sort/parallel/detail/indirect.hpp>
#include <boost/sort/parallel/detail/insertion_sort.hpp>
#include <boost/sort/parallel/detail/util/compare_traits.hpp>
#include <boost/sort/parallel/detail/util/nbits.hpp>
#include <boost/sort/parallel/detail/util/range.hpp>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <memory>
#include <type_traits>
#include <vector>

namespace boost
{
namespace sort
{
namespace parallel
{
namespace detail
{
//
//----------------------------------------------------------------------------
//                USING SENTENCES
//----------------------------------------------------------------------------
using util::range;
using util::nbits64;
using std::iterator_traits;

//-----------------------------------------------------------------------------
//  function : range_sort
/// @brief this function divide r_input in two parts, sort it,and merge moving
///        the elements to range_buf
/// @param range_input : range with the elements to sort
/// @param range_buffer : range with the elements sorted
/// @param comp : object for to compare two elements
/// @param level : when is 1, sort with the insertion_sort algorithm
///                if not make a recursive call splitting the ranges
//-----------------------------------------------------------------------------
template < class Iter1_t, class Iter2_t, class Compare >
void range_sort (const range< Iter1_t > &range_input,
                 const range< Iter2_t > &range_buffer, Compare comp,
                 uint32_t level)
{
    typedef range< Iter1_t > range_it1;
    typedef range< Iter2_t > range_it2;
    assert (range_input.size ( ) == range_buffer.size ( ) and level != 0);

    size_t nelem1 = (range_input.size ( ) + 1) >> 1;
    range_it1 range_input1 (range_input.first, range_input.first + nelem1),
        range_input2 (range_input.first + nelem1, range_input.last);

    if (level < 2) {
        insertion_sort (range_input1.first, range_input1.last, comp);
        insertion_sort (range_input2.first, range_input2.last, comp);
    }
    else
    {
        range_sort (range_it2 (range_buffer.first, range_buffer.first + nelem1),
                    range_input1, comp, level - 1);

        range_sort (range_it2 (range_buffer.first + nelem1, range_buffer.last),
                    range_input2, comp, level - 1);
    };

    full_merge (range_buffer, range_input1, range_input2, comp);
};

//---------------------------------------------------------------------------
/// @struct spin_sort
/// @brief  This class implement s stable sort algorithm with 1 thread, with
///         an auxiliary memory of N/2 elements
//----------------------------------------------------------------------------
template < class Iter_t, typename Compare = util::compare_iter< Iter_t > >
class spin_sort
{
    //------------------------------------------------------------------------
    //               DEFINITIONS AND CONSTANTS
    //------------------------------------------------------------------------
    typedef typename iterator_traits< Iter_t >::value_type value_t;
    typedef range< Iter_t > range_it;
    typedef range< value_t * > range_buf;
    // When the number of elements to sort is smaller than Sort_min, are sorted
    // by the insertion sort algorithm
    static const uint32_t Sort_min = 36;

    //------------------------------------------------------------------------
    //                      VARIABLES
    //------------------------------------------------------------------------
    // Pointer to the auxiliary memory
    value_t *ptr;

    // Number of elements in the auxiliary memory
    size_t nptr;

    // construct indicate if the auxiliary memory in initialized, and owner
    // indicate if the auxiliary memory had been created inside the object or
    // had
    // been received as a parameter
    bool construct = false, owner = false;

    //------------------------------------------------------------------------
    //                   PRIVATE FUNCTIONS
    //-------------------------------------------------------------------------
    spin_sort (Iter_t first, Iter_t last, Compare comp, value_t *paux,
               size_t naux);

public:
    //------------------------------------------------------------------------
    //                   PUBLIC FUNCTIONS
    //-------------------------------------------------------------------------
    spin_sort (Iter_t first, Iter_t last, Compare comp = Compare ( ))
        : spin_sort (first, last, comp, nullptr, 0){};

    spin_sort (Iter_t first, Iter_t last, Compare comp, range_buf range_aux)
        : spin_sort (first, last, comp, range_aux.first, range_aux.size ( )){};
    //
    //-----------------------------------------------------------------------
    //  function :~spin_sort
    /// @brief destructor of the struct. Destroy the elements if construct is
    /// true,
    ///        and return the memory if owner is true
    //-----------------------------------------------------------------------
    ~spin_sort (void)
    {
        if (construct) {
            destroy (range< value_t * > (ptr, ptr + nptr));
            construct = false;
        };
        if (owner and ptr != nullptr) std::return_temporary_buffer (ptr);
    };
}; //        End of class spin_sort
//
//############################################################################
//                                                                          ##
//              N O N    I N L I N E      F U N C T I O N S                 ##
//                                                                          ##
//                                                                          ##
//############################################################################
//
//-----------------------------------------------------------------------------
//  function : spin_sort
/// @brief constructor of the struct
//
/// @param first : iterator to the first element of the range to sort
/// @param last : iterator after the last element to the range to sort
/// @param comp : object for to compare two elements pointed by Iter_t
///               iterators
/// @param paux : pointer to the auxiliary memory provided. If nullptr, the
///               memory is created inside the class
/// @param naux : number of elements pointed by paux
//-----------------------------------------------------------------------------
template < class Iter_t, class Compare >
spin_sort< Iter_t, Compare >
  ::spin_sort (Iter_t first, Iter_t last, Compare comp, value_t *paux,
               size_t naux)
    : ptr (paux), nptr (naux), construct (false), owner (false)
{
    range< Iter_t > range_input (first, last);
    assert (range_input.valid ( ));

    size_t nelem = range_input.size ( );
    owner = construct = false;

    nptr = (nelem + 1) >> 1;
    size_t nelem_1 = nptr;
    size_t nelem_2 = nelem - nelem_1;

    if (nelem <= (Sort_min << 1)) {
        insertion_sort (range_input.first, range_input.last, comp);
        return;
    };
    //------------------- check if sort --------------------------------------
    bool sw = true;
    for (Iter_t it1 = range_input.first, it2 = range_input.first + 1;
         it2 != range_input.last and (sw = not comp (*it2, *it1)); it1 = it2++)
        ;
    if (sw) return;

    if (ptr == nullptr) {
        ptr = std::get_temporary_buffer< value_t > (nptr).first;
        if (ptr == nullptr) throw std::bad_alloc ( );
        owner = true;
    };
    range_buf range_aux (ptr, (ptr + nptr));

    //------------------------------------------------------------------------
    //                  Process
    //------------------------------------------------------------------------
    uint32_t nlevel = nbits64 (((nelem + Sort_min - 1) / Sort_min) - 1) - 1;
    assert (nlevel != 0);

    if ((nlevel & 1) == 1) {
        //---------------------------------------------------------------------
        // if the number of levels is odd, the data are in the first parameter
        // of range_sort, and the results appear in the second parameter
        //---------------------------------------------------------------------
        range_it range_1 (first, first + nelem_2),
            range_2 (first + nelem_2, last);
        range_aux = uninit_move (range_aux, range_2);
        construct = true;

        range_sort (range_aux, range_2, comp, nlevel);
        range_buf rng_bx (range_aux.first, range_aux.first + nelem_2);

        range_sort (range_1, rng_bx, comp, nlevel);
        half_merge (range_input, rng_bx, range_2, comp);
    }
    else
    {
        //--------------------------------------------------------------------
        // If the number of levels is even, the data are in the second
        // parameter of range_sort, and the results are in the same parameter
        //---------------------------------------------------------------------
        range_it range_1 (first, first + nelem_1),
            range_2 (first + nelem_1, last);
        range_aux = uninit_move (range_aux, range_1);
        construct = true;

        range_sort (range_1, range_aux, comp, nlevel);

        range_1.last = range_1.first + range_2.size ( );
        range_sort (range_1, range_2, comp, nlevel);
        half_merge (range_input, range_aux, range_2, comp);
    };
};

//****************************************************************************
}; //    End namespace algorithm
}; //    End namespace parallel
}; //    End namespace sort
}; //    End namepspace boost
//****************************************************************************
//
#endif
