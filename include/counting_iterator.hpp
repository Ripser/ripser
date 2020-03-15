#pragma once

#include <iterator>

namespace mrzv
{

class counting_iterator
{
    public:
        typedef size_t value_type;
        typedef const size_t& reference;
        typedef const size_t* pointer;
        typedef ptrdiff_t difference_type;
        typedef std::random_access_iterator_tag iterator_category;

        explicit counting_iterator(size_t x): m_inc(x)  {}
        size_t const& base() const { return m_inc; }
        reference operator*() const { return m_inc; }
        counting_iterator& operator++() { ++m_inc; return *this; }
        counting_iterator& operator--() { --m_inc; return *this; }

        counting_iterator& operator+=(difference_type i) { m_inc += i; return *this; }
        counting_iterator& operator-=(difference_type i) { m_inc -= i; return *this; }

        bool operator==(const counting_iterator& other) const   { return m_inc == other.m_inc; }
        bool operator!=(const counting_iterator& other) const   { return m_inc != other.m_inc; }

        counting_iterator operator+(difference_type x) const { counting_iterator result = *this; result += x; return result; }

        bool operator<(const counting_iterator& other) const { return m_inc < other.m_inc; }

        friend
        difference_type operator-(const counting_iterator& x, const counting_iterator& y)    { return x.m_inc - y.m_inc; }

        // and loads of others
    private:
        size_t m_inc;
};

}
