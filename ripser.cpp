#include <vector>
#include <iostream>
#include <iomanip>  
#include <fstream>
#include <iterator>
#include <string>
#include <cassert>
#include <queue>
#include <unordered_map>
#include "prettyprint.hpp"

typedef float value_t;
typedef long index_t;
typedef long coefficient_t;

//#define PRECOMPUTE_DIAMETERS
//#define PRECOMPUTE_DIAMETERS_IN_TOP_DIMENSION

#define USE_BINARY_SEARCH
#define USE_EXPONENTIAL_SEARCH

//#define ASSEMBLE_REDUCTION_MATRIX
//#define USE_COEFFICIENTS

//#define INDICATE_PROGRESS
#define PRINT_PERSISTENCE_PAIRS

#define FILE_FORMAT_DIPHA
//#define FILE_FORMAT_UPPER_TRIANGULAR_CSV
//#define FILE_FORMAT_LOWER_TRIANGULAR_CSV

class binomial_coeff_table {
    std::vector<std::vector<index_t> > B;
    index_t n_max, k_max;
    
public:
    binomial_coeff_table(index_t n, index_t k) {
        n_max = n;
        k_max = k;
    
        B.resize(n + 1);
        for( index_t i = 0; i <= n; i++ ) {
            B[i].resize(k + 1);
            for ( index_t j = 0; j <= std::min(i, k); j++ )
            {
                if (j == 0 || j == i)
                    B[i][j] = 1;
                else
                    B[i][j] = B[i-1][j-1] + B[i-1][j];
            }
        }
    }
    
    inline index_t operator()(index_t n, index_t k) const {
        assert(n <= n_max);
        assert(k <= k_max);
        return B[n][k];
    }
};

//
// https://comeoncodeon.wordpress.com/2011/10/09/modular-multiplicative-inverse/
//
std::vector<coefficient_t> multiplicative_inverse_vector (const coefficient_t m) {
    std::vector<coefficient_t> mod_inverse(m);
    mod_inverse[1] = 1;
    for(coefficient_t i = 2; i < m; ++i) {
        mod_inverse[i] = (-(m/i) * mod_inverse[m % i]) % m + m;
    }
    return mod_inverse;
}

template<typename OutputIterator>
inline OutputIterator get_simplex_vertices( index_t idx, const index_t dim, index_t n, const binomial_coeff_table& binomial_coeff, OutputIterator out  )
{
    --n;
    
    for( index_t k = dim + 1; k > 0; --k ) {
 
        #ifdef USE_BINARY_SEARCH
        if ( binomial_coeff( n , k ) > idx ) {
            index_t count;
            
            #ifdef USE_EXPONENTIAL_SEARCH
            for (count = 1; (binomial_coeff( n - count , k ) > idx); count = std::min(count << 1, n));
            #else
            count = n;
            #endif
            
            while (count > 0) {
                index_t i = n;
                index_t step = count >> 1;
                i -= step;
                if (binomial_coeff( i , k ) > idx) {
                    n = --i;
                    count -= step + 1;
                } else count = step;
            }
        }
        #else
        while( binomial_coeff( n , k ) > idx )
            --n;
        #endif
        
        assert( binomial_coeff( n , k ) <= idx );
        assert( binomial_coeff( n + 1, k ) > idx );
        
        *out++ = n;
        
        idx -= binomial_coeff( n , k );

    }

    return out;
}

std::vector<index_t> vertices_of_simplex(const index_t simplex_index, const index_t dim, const index_t n, const binomial_coeff_table& binomial_coeff) {
    std::vector<index_t> vertices;
    get_simplex_vertices( simplex_index, dim, n, binomial_coeff, std::back_inserter(vertices) );
    return vertices;
}

template <typename DistanceMatrix>
class rips_filtration_comparator {
public:
    const DistanceMatrix& dist;
    
    const index_t dim;

private:
    mutable std::vector<index_t> vertices;
    
    const binomial_coeff_table& binomial_coeff;
    
public:
    rips_filtration_comparator(
        const DistanceMatrix& _dist,
        const index_t _dim,
        const binomial_coeff_table& _binomial_coeff
    ): dist(_dist), dim(_dim), vertices(_dim + 1), binomial_coeff(_binomial_coeff) {};
    
    inline value_t diameter(const index_t index) const {
        value_t diam = 0;
        get_simplex_vertices(index, dim, dist.size(), binomial_coeff, vertices.begin() );
        
        for (index_t i = 0; i <= dim; ++i)
            for (index_t j = 0; j < i; ++j) {
                diam = std::max(diam, dist(vertices[i], vertices[j]));
            }
        return diam;
    }

    inline bool operator()(const index_t a, const index_t b) const
    {
        assert(a < binomial_coeff(dist.size(), dim + 1));
        assert(b < binomial_coeff(dist.size(), dim + 1));
        
        value_t a_diam = 0, b_diam = 0;
    
        b_diam = diameter(b);
        
        get_simplex_vertices(a, dim, dist.size(), binomial_coeff, vertices.begin() );
        for (index_t i = 0; i <= dim; ++i)
            for (index_t j = i + 1; j <= dim; ++j) {
                a_diam = std::max(a_diam, dist(vertices[i], vertices[j]));
                if (a_diam > b_diam)
                    return true;
            }
       
        return (a_diam == b_diam) && (a > b);
    }
    
    template <typename Entry>
    inline bool operator()(const Entry& a, const Entry& b) const
    {
        return operator()(get_index(a), get_index(b));
    }

};

#ifdef USE_COEFFICIENTS
struct entry_t {
    index_t index;
    coefficient_t value;
    
    entry_t(index_t _index, coefficient_t _value) :
        index(_index), value(_value) {}

    entry_t(index_t _index) :
        index(_index), value(1) {}
    
    entry_t() :
    index(0), value(1) {}
};

inline entry_t make_entry(index_t _index, coefficient_t _value) {
    return entry_t(_index, _value);
}

inline index_t get_index(entry_t e) {
    return e.index;
}

inline index_t get_coefficient(entry_t e) {
    return e.value;
}

std::ostream& operator<< (std::ostream& stream, const entry_t& e) {
    stream << get_index(e) << ":" << get_coefficient(e);
    return stream;
}

#else

typedef index_t entry_t;

inline index_t get_index(entry_t i) {
    return i;
}

inline index_t get_coefficient(entry_t i) {
    return 1;
}

inline entry_t make_entry(index_t _index, coefficient_t _value) {
    return entry_t(_index);
}


#endif

class simplex_coboundary_enumerator {
private:
    index_t idx, modified_idx, dim, n, k;
    
    const binomial_coeff_table& binomial_coeff;
    
public:
    inline simplex_coboundary_enumerator (
        index_t _idx, index_t _dim, index_t _n,
        const binomial_coeff_table& _binomial_coeff):
        idx(_idx), modified_idx(_idx), dim(_dim), k(dim + 1), n(_n - 1), binomial_coeff(_binomial_coeff) {}
    
    inline bool has_next()
    {
        while ( (k != -1 && n != -1) && (binomial_coeff( n , k ) <= idx) ) {
            idx -= binomial_coeff( n , k );
            
            modified_idx -= binomial_coeff( n , k );
            modified_idx += binomial_coeff( n , k + 1 );
            
            --n;
            --k;
        }
        return k != -1 && n != -1;
    }
    
    inline entry_t next()
    {
        while ( binomial_coeff( n , k ) <= idx ) {
            idx -= binomial_coeff( n , k );
            
            modified_idx -= binomial_coeff( n , k );
            modified_idx += binomial_coeff( n , k + 1 );
            
            --n;
        }
        return make_entry(
            modified_idx + binomial_coeff( n-- , k + 1 ),
            k & 1 ? 1 : -1
        );
    }
};

template <typename DistanceMatrix>
std::vector<value_t> get_diameters (
    const DistanceMatrix& dist, const index_t dim,
    const std::vector<value_t>& previous_diameters,
    const binomial_coeff_table& binomial_coeff
    )
{
    index_t n = dist.size();
    
    std::vector<value_t> diameters(binomial_coeff(n, dim + 1), 0);
    
    std::vector<index_t> coboundary;

    for (index_t simplex = 0; simplex < previous_diameters.size(); ++simplex) {
        coboundary.clear();
        
        #ifdef INDICATE_PROGRESS
        std::cout << "\033[Kpropagating diameter of simplex " << simplex + 1 << "/" << previous_diameters.size() << std::flush << "\r";
        #endif
        
        simplex_coboundary_enumerator coboundaries(simplex, dim - 1, n, binomial_coeff);
        while (coboundaries.has_next()) {
            index_t coface = get_index(coboundaries.next());
            diameters[coface] = std::max( diameters[coface], previous_diameters[simplex]);
        }
    }
    
    #ifdef INDICATE_PROGRESS
    std::cout << "\033[K";
    #endif
    
    return diameters;
}


class rips_filtration_diameter_comparator {
private:
    const std::vector<value_t>& diameters;
    
    const index_t dim;

public:
    std::vector<index_t> vertices;
    
    typedef value_t dist_t;
    
    const binomial_coeff_table& binomial_coeff;
    
public:
    rips_filtration_diameter_comparator(
        const std::vector<value_t>& _diameters,
        const index_t _dim,
        const binomial_coeff_table& _binomial_coeff
        ):
        diameters(_diameters), dim(_dim),
        binomial_coeff(_binomial_coeff) {}
    
    inline value_t diameter(const index_t a) const {
        assert(a < diameters.size());
        return diameters[a];
    }

    inline bool operator()(const index_t a, const index_t b) const
    {
        assert(a < diameters.size());
        assert(b < diameters.size());
        
        dist_t a_diam = diameters[a], b_diam = diameters[b];

        return ((a_diam > b_diam) || ((a_diam == b_diam) && (a > b)));
    }
    
    template <typename Entry>
    inline bool operator()(const Entry& a, const Entry& b) const
    {
        return operator()(get_index(a), get_index(b));
    }
};


class distance_matrix  {
public:
    typedef value_t value_type;
    std::vector<std::vector<value_t> > distances;
    inline value_t operator()(const index_t a, const index_t b) const {
        return distances[a][b];
    }
    
    inline size_t size() const {
        return distances.size();
    }
};


class compressed_upper_distance_matrix_adapter  {
public:
    typedef value_t value_type;
    std::vector<value_t>& distances;
    std::vector<value_t*> rows;
    
    index_t n;
    
    void init_distances () {
        distances.resize(n * (n-1) / 2);
    }
        
    void init_rows () {
        rows.resize(n);
        value_t* pointer = &distances[0] - 1;
        for (index_t i = 0; i < n - 1; ++i) {
            rows[i] = pointer;
            pointer += n - i - 2;
        }
    }
    
    compressed_upper_distance_matrix_adapter(std::vector<value_t>& _distances) :
        distances(_distances)
    {
        n = (1 + std::sqrt(1+ 8 * _distances.size())) / 2;
        
        assert( distances.size() == n * (n-1) / 2 );
        
        init_rows();
    }
    
    inline value_t operator()(const index_t a, const index_t b) const {
        if (a < b)
            return rows[a][b];
        else if (a > b)
            return rows[b][a];
        else
            return 0;
    }
    
    inline size_t size() const {
        return n;
    }
};


class compressed_lower_distance_matrix_adapter  {
public:
    typedef value_t value_type;
    std::vector<value_t>& distances;
    std::vector<value_t*> rows;
    
    index_t n;
    
    void init_distances () {
        distances.resize(n * (n-1) / 2);
    }
        
    void init_rows () {
        rows.resize(n);
        value_t* pointer = &distances[0];
        for (index_t i = 1; i < n; ++i) {
            rows[i] = pointer;
            pointer += i;
        }
    }
    
    compressed_lower_distance_matrix_adapter(
    std::vector<value_t>& _distances, const index_t _n) :
        distances(_distances), n(_n)
    {
        init_distances();
        init_rows();
    }

    compressed_lower_distance_matrix_adapter(std::vector<value_t>& _distances) :
        distances(_distances)
    {
        n = (1 + std::sqrt(1+ 8 * distances.size())) / 2;
        assert( distances.size() == n * (n-1) / 2 );
        
        init_rows();
    }
    
    template <typename DistanceMatrix>
    compressed_lower_distance_matrix_adapter(
    std::vector<value_t>& _distances, const DistanceMatrix& mat) :
        distances(_distances), n(mat.size()) {
        init_distances();
        init_rows();
        
        for (index_t i = 1; i < n; ++i) {
            for (index_t j = 0; j < i; ++j) {
                rows[i][j] = mat(i, j);
            }
        }
    }
    
    inline value_t operator()(const index_t i, const index_t j) const {
        if (i > j)
            return rows[i][j];
        else if (i < j)
            return rows[j][i];
        else
            return 0;
    }
    
    inline size_t size() const {
        return n;
    }
};

#ifdef USE_COEFFICIENTS
template <typename Heap>
inline entry_t pop_pivot(Heap& column, coefficient_t modulus)
{
    
    if( column.empty() )
        return entry_t(-1);
    else {
        index_t pivot_index = get_index(column.top());

        coefficient_t coefficient = 0;
        while( !column.empty() && get_index(column.top()) == pivot_index ) {
            coefficient_t new_coefficient = (coefficient + get_coefficient(column.top())) % modulus;
            assert(new_coefficient >= 0);
            coefficient = new_coefficient;
            column.pop();
            
            if( coefficient == 0 )
                pivot_index = column.empty() ? -1 : get_index(column.top());
        }
        return make_entry(pivot_index, coefficient);
    }
}

template <typename Heap>
inline entry_t get_pivot(Heap& column, coefficient_t modulus)
{
    entry_t max_element = pop_pivot(column, modulus);
    if( get_index(max_element) != -1 )
        column.push( max_element );
    return max_element;
}

#else

template <typename Heap>
inline entry_t pop_pivot(Heap& column, coefficient_t modulus)
{
    
    if( column.empty() )
        return -1;
    else {
        index_t pivot_index = get_index(column.top());
        column.pop();
        while( !column.empty() && column.top() == pivot_index ) {
            column.pop();
            if( column.empty() )
                return -1;
            else {
                pivot_index = column.top();
                column.pop();
            }
        }
        return pivot_index;
    }
}

template <typename Heap>
inline entry_t get_pivot(Heap& column, coefficient_t modulus)
{
    entry_t max_element = pop_pivot(column, modulus);
    if( get_index(max_element) != -1 )
        column.push( max_element );
    return max_element;
}

#endif

template <typename Comparator>
void assemble_columns_to_reduce (
    std::vector<index_t>& columns_to_reduce,
    std::unordered_map<index_t, index_t>& pivot_column_index,
    const Comparator& comp,
    index_t dim, index_t n,
    value_t threshold,
    const binomial_coeff_table& binomial_coeff
) {
    index_t num_simplices = binomial_coeff(n, dim + 2);
    
    columns_to_reduce.clear();
    
    for (index_t index = 0; index < num_simplices; ++index) {

        if (comp.diameter(index) <= threshold && pivot_column_index.find(index) == pivot_column_index.end()) {
            columns_to_reduce.push_back(index);
        }
    }

    std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), comp);
}


template <typename ValueType>
class compressed_sparse_matrix  {
public:
    std::vector<size_t> bounds;
    std::vector<ValueType> entries;
    
    
    inline size_t size() const {
        return bounds.size();
    }
    
    inline typename std::vector<ValueType>::const_iterator cbegin(size_t index) const {
        assert(index < size());
        return index == 0 ? entries.cbegin() : entries.cbegin() + bounds[index - 1];
    }

    inline typename std::vector<ValueType>::const_iterator cend(size_t index) const {
        assert(index < size());
        return entries.cbegin() + bounds[index];
    }
    
    template <typename Iterator>
    inline void append(Iterator begin, Iterator end) {
        for (Iterator it = begin; it != end; ++it) {
            entries.push_back(*it);
        }
        bounds.push_back(entries.size());
    }

    inline void append() {
        bounds.push_back(entries.size());
    }

    inline void push_back(ValueType e) {
        assert(0 < size());
        entries.push_back(e);
        ++bounds.back();
    }
    
    inline void pop_back() {
        assert(0 < size());
        entries.pop_back();
        --bounds.back();
    }
    
    template <typename Collection>
    inline void append(const Collection collection) {
        append(collection.cbegin(), collection.cend());
    }

};


template <typename Heap>
inline std::vector<entry_t> get_column_vector(Heap column, coefficient_t modulus)
{
    std::vector<entry_t> temp_col;
    entry_t pivot = pop_pivot( column, modulus );
    while( get_index(pivot) != -1 ) {
        temp_col.push_back( pivot );
        pivot = pop_pivot( column, modulus );
    }
    return temp_col;
}


template <typename Heap>
inline std::vector<entry_t> get_heap_vector(Heap heap)
{
    std::vector<entry_t> temp_col;
    while( !heap.empty() ) {
        temp_col.push_back( heap.top() );
        heap.pop();
    }
    return temp_col;
}





template <typename ComparatorCofaces, typename Comparator>
void compute_pairs(
    std::vector<index_t>& columns_to_reduce,
    std::unordered_map<index_t, index_t>& pivot_column_index,
    const ComparatorCofaces& comp, const Comparator& comp_prev,
    index_t dim, index_t n,
    value_t threshold, coefficient_t modulus,
    const binomial_coeff_table& binomial_coeff,
    const std::vector<coefficient_t>& multiplicative_inverse
) {
    
    std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
    
    
    #ifdef ASSEMBLE_REDUCTION_MATRIX
    compressed_sparse_matrix <entry_t> reduction_matrix;
    #else
    #ifdef USE_COEFFICIENTS
    std::vector <entry_t> reduction_coefficients;
    #endif
    #endif
    
    for (index_t i = 0; i < columns_to_reduce.size(); ++i) {
        index_t column_to_reduce = columns_to_reduce[i];
        
        #ifdef ASSEMBLE_REDUCTION_MATRIX
        std::priority_queue<entry_t, std::vector<entry_t>, decltype(comp_prev) > reduction_column(comp_prev);
        #endif
        
        std::priority_queue<entry_t, std::vector<entry_t>, decltype(comp) > working_coboundary(comp);
        
        #ifdef INDICATE_PROGRESS
        std::cout << "\033[K" << "reducing column " << i + 1 << "/" << columns_to_reduce.size()
        << " (diameter " << comp_prev.diameter(column_to_reduce) << ")"
        << std::flush << "\r";
        #endif
        
        index_t j = i;
        
        index_t column_to_add = column_to_reduce;
        
        
        // start with a pivot entry with coefficient -1 in order to initialize working_coboundary
        // with the coboundary of the simplex with index column_to_reduce
        
        entry_t pivot = make_entry(column_to_reduce, -1);
        
        
        std::vector<index_t> coboundary;
        
//        std::cout << "reducing " << column_to_reduce << ": pivot ";
        
        #ifdef ASSEMBLE_REDUCTION_MATRIX
        reduction_matrix.append();
        #endif
        

        // initialize reduction_matrix as identity matrix
        #ifdef ASSEMBLE_REDUCTION_MATRIX
        reduction_matrix.push_back(make_entry(column_to_reduce, 1));
        #else
        #ifdef USE_COEFFICIENTS
        reduction_coefficients.push_back(make_entry(column_to_reduce, 1));
        #endif
        #endif
        
        do {
            
            const coefficient_t factor = modulus - get_coefficient(pivot);


//            std::priority_queue<index_t, std::vector<entry_t>, decltype(comp) > eliminating_coboundary(comp);

//            std::cout << "w:" << get_column_vector(working_coboundary, modulus) << std::endl;

            #ifdef ASSEMBLE_REDUCTION_MATRIX
            
            for (auto it = reduction_matrix.cbegin(j); it != reduction_matrix.cend(j); ++it) {
                const entry_t& simplex = *it;
    
                reduction_column.push( simplex );
            
                simplex_coboundary_enumerator cofaces(get_index(simplex), dim, n, binomial_coeff);
                while (cofaces.has_next()) {
                    entry_t coface = cofaces.next();
                    
                    index_t coface_index = get_index(coface);
                    if (comp.diameter(coface_index) <= threshold) {
                        coefficient_t coface_coefficient = get_coefficient(coface) + modulus;
                        coface_coefficient %= modulus;
                        
                        coface_coefficient *= get_coefficient(simplex);
                        coface_coefficient %= modulus;
                        
                        coface_coefficient *= factor;
                        coface_coefficient %= modulus;
                        
                        assert(coface_coefficient >= 0);
                        
                        entry_t e = make_entry(coface_index, coface_coefficient);
                        working_coboundary.push(e);
//                        eliminating_coboundary.push(e);
                    }
                }
            }

            #else
            
            #ifdef USE_COEFFICIENTS
            const entry_t& simplex = reduction_coefficients[j];
            
            simplex_coboundary_enumerator cofaces(get_index(simplex), dim, n, binomial_coeff);
            while (cofaces.has_next()) {
                entry_t coface = cofaces.next();
                
                index_t coface_index = get_index(coface);
                if (comp.diameter(coface_index) <= threshold) {
                    coefficient_t coface_coefficient = get_coefficient(coface) + modulus;
                    coface_coefficient %= modulus;
                    
                    coface_coefficient *= get_coefficient(simplex);
                    coface_coefficient %= modulus;
                    
                    coface_coefficient *= factor;
                    coface_coefficient %= modulus;
                    
                    entry_t e = make_entry(coface_index, coface_coefficient);
                    working_coboundary.push(e);
//                    eliminating_coboundary.push(e);
                }
            }
            
            #else
            
            simplex_coboundary_enumerator cofaces(column_to_add, dim, n, binomial_coeff);
            while (cofaces.has_next()) {
                index_t coface_index = cofaces.next();
                if (comp.diameter(coface_index) <= threshold) {
                    working_coboundary.push(coface_index);
//                    eliminating_coboundary.push(e);
                }
            }
            #endif
            
            #endif


            
//            std::cout << get_heap_vector(working_coboundary) << std::endl;

//            std::cout << "e:" << get_column_vector(eliminating_coboundary, modulus) << std::endl;
//            std::cout << "w:" << get_column_vector(working_coboundary, modulus) << std::endl << std::endl;
            
            pivot = get_pivot(working_coboundary, modulus);
            
//            std::cout << get_index(pivot) << " ";
        
            if (get_index(pivot) != -1) {
                auto pair = pivot_column_index.find(get_index(pivot));
                
                if (pair == pivot_column_index.end()) {
//                    std::cout << std::endl;
                    
                    pivot_column_index.insert(std::make_pair(get_index(pivot), i));
                    
                    #ifdef USE_COEFFICIENTS
                    const coefficient_t inverse = multiplicative_inverse[ get_coefficient( pivot ) ];
                    #else
                    const coefficient_t inverse = 1;
                    #endif
                    
                    // replace current column of reduction_matrix (with a single diagonal 1 entry)
                    // by reduction_column (possibly with a different entry on the diagonal)
                    #ifdef ASSEMBLE_REDUCTION_MATRIX
                    reduction_matrix.pop_back();
                    while (true) {
                        entry_t e = pop_pivot(reduction_column, modulus);
                        index_t index = get_index(e);
                        if (index == -1)
                            break;
                        const coefficient_t coefficient = inverse * get_coefficient(e) % modulus;
                        assert(coefficient > 0);
                        reduction_matrix.push_back(make_entry(index, coefficient));
                    }
                    #else
                    #ifdef USE_COEFFICIENTS
                    reduction_coefficients.pop_back();
                    reduction_coefficients.push_back(make_entry(column_to_reduce, inverse));
                    #endif
                    #endif
                    
                    #ifdef PRINT_PERSISTENCE_PAIRS
                    value_t birth = comp_prev.diameter(column_to_reduce), death = comp.diameter(get_index(pivot));
                    if (birth != death) {
                        #ifdef INDICATE_PROGRESS
                        std::cout << "\033[K";
                        #endif
                        std::cout << " [" << birth << "," << death << ")" << std::endl << std::flush;
                    }
                    #endif

                    break;
                }
                
                j = pair->second;
                column_to_add = columns_to_reduce[j];
            }

        } while ( get_index(pivot) != -1 );
        
        #ifdef PRINT_PERSISTENCE_PAIRS
        if ( get_index(pivot) == -1 ) {
//            std::cout << std::endl;
        
            value_t birth = comp_prev.diameter(column_to_reduce);
            #ifdef INDICATE_PROGRESS
            std::cout << "\033[K";
            #endif
            std::cout << " [" << birth << ", )" << std::endl << std::flush;
        }
        #endif

    
//        #ifdef ASSEMBLE_REDUCTION_MATRIX
//        std::cout << "reduction matrix fill-in: " << i + 1 << " + " << reduction_matrix.entries.size() - (i + 1) << std::endl;
//        #endif
    
    }

    #ifdef INDICATE_PROGRESS
    std::cout << "\033[K";
    #endif

}

bool is_prime(const long n) {
    bool is_prime = true;
    for (int i = 2; i <= n/2; ++i)
        if (n%i == 0) {
            is_prime = false;
            break;
        }
    return is_prime;
}

void print_usage_and_exit(int exit_code)
{
    std::cerr << "Usage: " << "ripser " << "[options] filename" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  --help           print this screen" << std::endl;
    std::cerr << "  --top_dim <k>    compute persistent homology up to dimension <k>" << std::endl;
    std::cerr << "  --threshold <t>  compute Rips complexes up to diameter <t>" << std::endl;
    #ifdef USE_COEFFICIENTS
    std::cerr << "  --modulus <p>    compute homology with coefficients in the prime field Z/<p>Z" << std::endl;
    #endif
    
    exit(exit_code);
}

int main( int argc, char** argv ) {
    
    if( argc < 2 ) print_usage_and_exit(-1);

    const char *filename = nullptr;
    
    index_t dim_max = 1;
    value_t threshold = std::numeric_limits<value_t>::max();

    #ifdef USE_COEFFICIENTS
    coefficient_t modulus = 2;
    #else
    const coefficient_t modulus = 2;
    #endif

    for( index_t i = 1; i < argc; ++i ) {
        const std::string arg(argv[ i ]);
        if( arg == "--help" ) {
            print_usage_and_exit(0);
        } else if( arg == "--top_dim" ) {
            std::string parameter = std::string( argv[ ++i ] );
            size_t next_pos;
            dim_max = std::stol( parameter, &next_pos );
            if( next_pos != parameter.size() )
                print_usage_and_exit( -1 );
        } else if( arg == "--threshold" ) {
            std::string parameter = std::string( argv[ ++i ] );
            size_t next_pos;
            threshold = std::stof( parameter, &next_pos );
            if( next_pos != parameter.size() )
                print_usage_and_exit( -1 );
        #ifdef USE_COEFFICIENTS
        } else if( arg == "--modulus" ) {
            std::string parameter = std::string( argv[ ++i ] );
            size_t next_pos;
            modulus = std::stol( parameter, &next_pos );
            if( next_pos != parameter.size() || !is_prime(modulus) )
                print_usage_and_exit( -1 );
        #endif
        } else {
            if (filename) {
                print_usage_and_exit( -1 );
            }
            filename = argv[i];
        }
    }


    std::ifstream input_stream( filename, std::ios_base::binary | std::ios_base::in );
    if( input_stream.fail( ) ) {
        std::cerr << "couldn't open file" << filename << std::endl;
        exit(-1);
    }
    
    std::vector<std::vector<value_t>> diameters(dim_max + 2);
    
    #ifdef FILE_FORMAT_DIPHA
    
    int64_t magic_number;
    input_stream.read( reinterpret_cast<char*>(&magic_number), sizeof( int64_t ) );
    if( magic_number != 8067171840 ) {
        std::cerr << filename << " is not a Dipha file (magic number: 8067171840)" << std::endl;
        exit(-1);
    }

    int64_t file_type;
    input_stream.read( reinterpret_cast<char*>(&file_type), sizeof( int64_t ) );
    if( file_type != 7 ) {
        std::cerr << filename << " is not a Dipha distance matrix (file type: 7)" << std::endl;
        exit(-1);
    }
   
    int64_t n;
    input_stream.read( reinterpret_cast<char*>(&n), sizeof( int64_t ) );
    
    distance_matrix dist_full;
    dist_full.distances = std::vector<std::vector<value_t>>(n, std::vector<value_t>(n));
    
    for( int i = 0; i < n; ++i ) {
        for( int j = 0; j < n; ++j ) {
            double val;
            input_stream.read( reinterpret_cast<char*>(&val), sizeof(double) );
            dist_full.distances[i][j] = val;
        }
    }
    
    std::cout << "distance matrix with " << n << " points" << std::endl;

    compressed_lower_distance_matrix_adapter dist(diameters[1], dist_full);

    std::cout << "distance matrix transformed to lower triangular form" << std::endl;

    #endif


    #ifdef FILE_FORMAT_UPPER_TRIANGULAR_CSV
    
    std::vector<value_t> distances;
    std::string value_string;
    while(std::getline(input_stream, value_string, ','))
        distances.push_back(std::stod(value_string));
    
    compressed_upper_distance_matrix_adapter dist_upper(distances);
    
    index_t n = dist_upper.size();
    
    std::cout << "distance matrix with " << n << " points" << std::endl;
 
    compressed_lower_distance_matrix_adapter dist(diameters[1], dist_upper);

    std::cout << "distance matrix transformed to lower triangular form" << std::endl;
    
    #endif
    

    #ifdef FILE_FORMAT_LOWER_TRIANGULAR_CSV
    
    std::vector<value_t>& distances = diameters[1];
    std::string value_string;
    while(std::getline(input_stream, value_string, ','))
        distances.push_back(std::stod(value_string));
    
    compressed_lower_distance_matrix_adapter dist(distances);
    
    index_t n = dist.size();
    
    std::cout << "distance matrix with " << n << " points" << std::endl;
     
    #endif



    auto result = std::minmax_element(dist.distances.begin(), dist.distances.end());
    std::cout << "value range: [" << *result.first << "," << *result.second << "]" << std::endl;
    

    assert(dim_max < n - 2);
    
    binomial_coeff_table binomial_coeff(n, dim_max + 2);
    
    std::vector<coefficient_t> multiplicative_inverse(multiplicative_inverse_vector(modulus));
    
    
    std::vector<index_t> columns_to_reduce;
    
    
    for (index_t index = n; index-- > 0; ) {
        columns_to_reduce.push_back(index);
    }

    

    index_t dim = 0;

    {
        rips_filtration_diameter_comparator comp(diameters[1], dim + 1, binomial_coeff);
        rips_filtration_comparator<decltype(dist)> comp_prev(dist, dim, binomial_coeff);

        std::unordered_map<index_t, index_t> pivot_column_index;

        compute_pairs(
            columns_to_reduce,
            pivot_column_index,
            comp, comp_prev,
            dim, n,
            threshold, modulus,
            binomial_coeff, multiplicative_inverse
        );
        
        assemble_columns_to_reduce(
            columns_to_reduce,
            pivot_column_index,
            comp,
            dim, n,
            threshold,
            binomial_coeff
        );
    }
    
    #ifdef PRECOMPUTE_DIAMETERS
    #ifdef PRECOMPUTE_DIAMETERS_IN_TOP_DIMENSION
    for (dim = 1; dim <= dim_max; ++dim) {
    #else
    for (dim = 1; dim < dim_max; ++dim) {
    #endif
    #else
    for (dim = 1; dim <= dim_max; ++dim) {
    #endif
    
        #ifdef PRECOMPUTE_DIAMETERS
        
        #ifdef INDICATE_PROGRESS
        std::cout << "precomputing " << dim + 1 << "-simplex diameters" << std::endl;
        #endif
        diameters[dim + 1] = get_diameters( dist, dim + 1, diameters[dim], binomial_coeff );
    
        rips_filtration_diameter_comparator comp(diameters[dim + 1], dim + 1, binomial_coeff);
        rips_filtration_diameter_comparator comp_prev(diameters[dim], dim, binomial_coeff);
        
        #else
        
        rips_filtration_comparator<decltype(dist)> comp(dist, dim + 1, binomial_coeff);
        rips_filtration_comparator<decltype(dist)> comp_prev(dist, dim, binomial_coeff);
        
        #endif

        std::unordered_map<index_t, index_t> pivot_column_index;
    
        compute_pairs(
            columns_to_reduce,
            pivot_column_index,
            comp, comp_prev,
            dim, n,
            threshold, modulus,
            binomial_coeff, multiplicative_inverse
        );
        
        assemble_columns_to_reduce(
            columns_to_reduce,
            pivot_column_index,
            comp,
            dim, n,
            threshold,
            binomial_coeff
        );
        
//        if ( dim > 1 )
//            diameters[dim] = std::vector<value_t>();
    }
        
    #ifdef PRECOMPUTE_DIAMETERS
    #ifndef PRECOMPUTE_DIAMETERS_IN_TOP_DIMENSION
    {
        dim = dim_max;
        
        rips_filtration_diameter_comparator comp_prev(diameters[dim], dim, binomial_coeff);
        rips_filtration_comparator<decltype(dist)> comp(dist, dim + 1, binomial_coeff);
        
        std::unordered_map<index_t, index_t> pivot_column_index;

        compute_pairs(
            columns_to_reduce,
            pivot_column_index,
            comp, comp_prev,
            dim, n,
            threshold, modulus,
            binomial_coeff, multiplicative_inverse
        );
    }
    #endif
    #endif
    
}
