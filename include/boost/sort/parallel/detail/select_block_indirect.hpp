//----------------------------------------------------------------------------
/// @file select_block_indirect.hpp
/// @brief block indirect sort algorithm
///
/// @author Copyright (c) 2016 Francisco Jose Tapia (fjtapia@gmail.com )\n
///         Distributed under the Boost Software License, Version 1.0.\n
///         ( See accompanying file LICENSE_1_0.txt or copy at
///           http://www.boost.org/LICENSE_1_0.txt  )
/// @version 0.1
///
/// @remarks
//-----------------------------------------------------------------------------
#ifndef __BOOST_SORT_PARALLEL_DETAIL_SELECT_BLOCK_INDIRECT_HPP
#define __BOOST_SORT_PARALLEL_DETAIL_SELECT_BLOCK_INDIRECT_HPP

#include <boost/sort/parallel/detail/block_indirect_sort.hpp>
#include <boost/sort/parallel/detail/util/compare_traits.hpp>
#include <boost/sort/parallel/detail/util/nbits.hpp>

namespace boost
{
namespace sort
{
namespace parallel
{
namespace detail
{
using boost::sort::parallel::detail::util::tmsb;
//
///---------------------------------------------------------------------------
//  function select_block_indirect
/// @brief This class is select the block size in the block_indirect_sort
///        algorithm depending of the type and size of the data to sort
///
//----------------------------------------------------------------------------
template < class Iter_t, class Compare,
           util::enable_if_string< util::value_iter< Iter_t > > * = nullptr >
void select_block_indirect (Iter_t first, Iter_t last, Compare cmp,
                           uint32_t nthr = std::thread::hardware_concurrency ())
{
    block_indirect_sort< 128, 128, Iter_t, Compare > (first, last, cmp, nthr);
};

template < size_t Size >
struct block_size
{
    static constexpr const uint32_t BitsSize =
        (Size == 0) ? 0 : (Size > 256) ? 9 : tmsb[Size - 1];
    static constexpr const uint32_t sz[10] = {4096, 4096, 4096, 4096, 2048,
                                              1024, 768,  512,  256,  128};
    static constexpr const uint32_t data = sz[BitsSize];
};
//
///---------------------------------------------------------------------------
/// @struct select_block_indirect_sort
/// @brief This class is select the block size in the block_indirect_sort
///        algorithm depending of the type and size of the data to sort
///
//----------------------------------------------------------------------------
template < class Iter_t, class Compare,
          util::enable_if_not_string< util::value_iter< Iter_t > >* = nullptr >

void select_block_indirect (Iter_t first, Iter_t last, Compare cmp,
			  uint32_t nthr = std::thread::hardware_concurrency ())
{
    block_indirect_sort<block_size< sizeof (util::value_iter< Iter_t >) >::data,
	                28, Iter_t, Compare > (first, last, cmp, nthr);
};

//
//****************************************************************************
}; //    End namespace detail
}; //    End namespace parallel
}; //    End namespace sort
}; //    End namespace boost
//****************************************************************************
//
#endif
