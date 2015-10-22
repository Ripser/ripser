#include <vector>
#include <iostream>
#include <iomanip>  
#include <fstream>
#include <string>
#include <cassert>
#include <queue>
#include <unordered_map>
#include "prettyprint.hpp"
//#include <boost/numeric/ublas/symmetric.hpp>

typedef float value_t;
typedef long index_t;

#define PRECOMPUTE_DIAMETERS

#define PRECOMPUTE_DIAMETERS_IN_TOP_DIMENSION

#define USE_BINARY_SEARCH

#define INDICATE_PROGRESS

#define PRINT_PERSISTENCE_PAIRS


//#define ASSEMBLE_REDUCTION_COLUMN

//#define FILE_FORMAT_DIPHA
#define FILE_FORMAT_CSV


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
//        std::cout << "B(" << n << "," << k << ")\n";
        assert(n <= n_max);
        assert(k <= k_max);
        return B[n][k];
    }
};



template<typename OutputIterator>
OutputIterator get_simplex_vertices( index_t idx, const index_t dim, index_t n, const binomial_coeff_table& binomial_coeff, OutputIterator out  )
{
    --n;
    
    for( index_t k = dim + 1; k > 0; --k ) {
 
        #ifdef USE_BINARY_SEARCH
        if ( binomial_coeff( n , k ) > idx ) {
            index_t count;
            
            for (count = 1; (binomial_coeff( n - count , k ) > idx); count = std::min(count << 1, n));
            
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
    
    typedef decltype(dist(0,0)) dist_t;
    
    bool reverse;
    
    const binomial_coeff_table& binomial_coeff;
    
public:
    rips_filtration_comparator(
    const DistanceMatrix& _dist,
    const index_t _dim,
    const binomial_coeff_table& _binomial_coeff
//    ,
//    bool _reverse = true
    ):
    dist(_dist), dim(_dim), vertices(_dim + 1),
    binomial_coeff(_binomial_coeff)
//    ,
//    reverse(_reverse)
    {};
    
    inline dist_t diameter(const index_t index) const {
        dist_t diam = 0;
        get_simplex_vertices(index, dim, dist.size(), binomial_coeff, vertices.begin() );
        
        for (index_t i = 0; i <= dim; ++i)
            for (index_t j = i + 1; j <= dim; ++j) {
                diam = std::max(diam, dist(vertices[i], vertices[j]));
            }
        return diam;
    }

    inline bool operator()(const index_t a, const index_t b) const
    {
        assert(a < binomial_coeff(dist.size(), dim + 1));
        assert(b < binomial_coeff(dist.size(), dim + 1));
        
        dist_t a_diam = 0, b_diam = 0;

    
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
    
};

template <typename DistanceMatrix, index_t dim>
class rips_filtration_comparator_dim;

template <typename DistanceMatrix>
class rips_filtration_comparator_dim <DistanceMatrix, 2> {
public:
    const DistanceMatrix& dist;
    
private:
    mutable std::vector<index_t> vertices;
    
    typedef decltype(dist(0,0)) dist_t;
    
    bool reverse;
    
    const binomial_coeff_table& binomial_coeff;
    
public:
    rips_filtration_comparator_dim(
    const DistanceMatrix& _dist,
    const binomial_coeff_table& _binomial_coeff
    ):
    dist(_dist), vertices(3),
    binomial_coeff(_binomial_coeff)
    {};
    
    inline dist_t diameter(const index_t index) const {
        dist_t diam = 0;
        get_simplex_vertices(index, 2, dist.size(), binomial_coeff, vertices.begin() );
        
        diam = std::max(diam, dist(vertices[0], vertices[1]));
        diam = std::max(diam, dist(vertices[0], vertices[2]));
        diam = std::max(diam, dist(vertices[0], vertices[3]));
        diam = std::max(diam, dist(vertices[1], vertices[2]));
        diam = std::max(diam, dist(vertices[1], vertices[3]));
        diam = std::max(diam, dist(vertices[2], vertices[3]));

        return diam;
    }

    inline bool operator()(const index_t a, const index_t b) const
    {
        assert(a < binomial_coeff(dist.size(), 3));
        assert(b < binomial_coeff(dist.size(), 3));
        
        dist_t a_diam = 0, b_diam = 0;

    
        b_diam = diameter(b);
        
        get_simplex_vertices(a, 2, dist.size(), binomial_coeff, vertices.begin() );
        if ((a_diam = dist(vertices[0], vertices[1])) >= b_diam) return (a_diam > b_diam) || (a > b);
        if ((a_diam = dist(vertices[0], vertices[2])) >= b_diam) return (a_diam > b_diam) || (a > b);
        if ((a_diam = dist(vertices[0], vertices[3])) >= b_diam) return (a_diam > b_diam) || (a > b);
        if ((a_diam = dist(vertices[1], vertices[2])) >= b_diam) return (a_diam > b_diam) || (a > b);
        if ((a_diam = dist(vertices[1], vertices[3])) >= b_diam) return (a_diam > b_diam) || (a > b);
        if ((a_diam = dist(vertices[2], vertices[3])) >= b_diam) return (a_diam > b_diam) || (a > b);
        
        return false;
    }
    
};


template<typename OutputIterator>
void inline get_simplex_coboundary( index_t idx, index_t dim, index_t n, const binomial_coeff_table& binomial_coeff, OutputIterator coboundary )
{

    --n;
    
    index_t modified_idx = idx;
    
    for( index_t k = dim + 1; k != -1 && n != -1; --k ) {
        while( binomial_coeff( n , k ) > idx ) {
            *coboundary++ = modified_idx + binomial_coeff( n , k + 1 );
            if (n==0) break;
            --n;
        }
        
        idx -= binomial_coeff( n , k );
            
        modified_idx -= binomial_coeff( n , k );
        modified_idx += binomial_coeff( n , k + 1 );
        
        --n;
    }
 
    return;
}

class simplex_coboundary_enumerator {
private:
    index_t idx, modified_idx, dim, n, k;
    
    const binomial_coeff_table& binomial_coeff;
    
public:
    inline simplex_coboundary_enumerator (
        index_t _idx,
        index_t _dim,
        index_t _n,
        const binomial_coeff_table& _binomial_coeff):
        idx(_idx), modified_idx(_idx), dim(_dim), k(dim + 1), n(_n - 1), binomial_coeff(_binomial_coeff)
    {
        
    }
    
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
    
    inline index_t next()
    {
        while ( binomial_coeff( n , k ) <= idx ) {
            idx -= binomial_coeff( n , k );
            
            modified_idx -= binomial_coeff( n , k );
            modified_idx += binomial_coeff( n , k + 1 );
            
            --n;
        }
        return modified_idx + binomial_coeff( n-- , k + 1 );
    }
};

template <typename DistanceMatrix>
std::vector<value_t> get_diameters (
    const DistanceMatrix& dist,
    const index_t dim,
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
            index_t coface = coboundaries.next();
            diameters[coface] = std::max( diameters[coface], previous_diameters[simplex]);
        }
    }
    
    #ifdef INDICATE_PROGRESS
    std::cout << "\033[K";
    #endif
    
//    rips_filtration_comparator<decltype(dist)> comp(dist, dim, binomial_coeff);
//    for (index_t simplex = 0; simplex < diameters.size(); ++simplex) {
//        assert(diameters[simplex] == comp.diameter(simplex));
//    }

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
        binomial_coeff(_binomial_coeff)
    {};
    
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


class compressed_upper_distance_matrix  {
public:
    typedef value_t value_type;
    std::vector<value_t> distances;
    std::vector<value_t*> rows;
    
    index_t n;
    
    compressed_upper_distance_matrix(std::vector<value_t>&& row_vector) noexcept {
        n = (1 + std::sqrt(1+ 8 * row_vector.size())) / 2;
        
        distances = std::move(row_vector);
        
        rows.resize(n);
        
        value_t* pointer = &distances[0] - 1;
        for (index_t i = 0; i < n - 1; ++i) {
            rows[i] = pointer;
            pointer += n - i - 2;
        }
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


class compressed_lower_distance_matrix  {
public:
    typedef value_t value_type;
    std::vector<value_t> distances;
    std::vector<value_t*> rows;
    
    index_t n;
    
    void init () {
        rows.resize(n);
        distances.resize(n * (n-1) / 2);
        
        value_t* pointer = &distances[0];
        for (index_t i = 1; i < n; ++i) {
            rows[i] = pointer;
            pointer += i;
        }
    }
    
    compressed_lower_distance_matrix(const index_t _n) : n(_n) {
        init();
    }

    compressed_lower_distance_matrix(std::vector<value_t>&& row_vector) {
        n = (1 + std::sqrt(1+ 8 * row_vector.size())) / 2;
        
        distances = std::move(row_vector);
        
        init();
    }
    
//    template <typename DistanceMatrix>
//    compressed_lower_distance_matrix(const DistanceMatrix& mat) : compressed_lower_distance_matrix(mat.size()) {
    compressed_lower_distance_matrix(const compressed_upper_distance_matrix& mat) {
        n = mat.size();
        
        init();
    
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

    
    template <typename Collection>
    inline void append(const Collection collection) {
        append(collection.cbegin(), collection.cend());
    }

};

template <typename Heap>
inline index_t pop_pivot(Heap& column)
{
    if( column.empty() )
        return -1;
    else {
        index_t max_element = column.top();
        column.pop();
        while( !column.empty() && column.top() == max_element ) {
            column.pop();
            if( column.empty() )
                return -1;
            else {
                max_element = column.top();
                column.pop();
            }
        }
        return max_element;
    }
}

template <typename Heap>
inline index_t get_pivot(Heap& column)
{
    index_t max_element = pop_pivot(column);
    if( max_element != -1 )
        column.push( max_element );
    return max_element;
}

template <typename Heap>
inline std::vector<index_t> move_to_column_vector(Heap& column)
{
    std::vector<index_t> temp_col;
    index_t pivot = pop_pivot( column );
    while( pivot != -1 ) {
        temp_col.push_back( pivot );
        pivot = pop_pivot( column );
    }
    return temp_col;
}



void print_help_and_exit()
{
    std::cerr << "Usage: " << "ripser " << "[options] input_filename output_filename" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << std::endl;
    std::cerr << "--help        --  prints this screen" << std::endl;
    std::cerr << "--top_dim N   --  maximal dimension to compute" << std::endl;
    std::cerr << "--threshold D --  maximal diameter to compute" << std::endl;
    
    exit(-1);
}

template <typename ComparatorCofaces, typename Comparator>
void compute_pairs(
    std::vector<index_t>& columns_to_reduce,
    std::unordered_map<index_t, index_t>& pivot_column_index,
    const ComparatorCofaces& comp, const Comparator& comp_prev,
    index_t dim, index_t dim_max, index_t n,
    value_t threshold,
    const binomial_coeff_table& binomial_coeff
) {
    
        std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
        
        for (index_t i = 0; i < columns_to_reduce.size(); ++i) {
            index_t index = columns_to_reduce[i];
            
            #ifdef ASSEMBLE_REDUCTION_COLUMN
            std::priority_queue<index_t, std::vector<index_t>, decltype(comp) > reduction_column(comp);
            #endif
            
            std::priority_queue<index_t, std::vector<index_t>, decltype(comp) > working_coboundary(comp);
            
            #ifdef INDICATE_PROGRESS
            std::cout << "\033[K" << "reducing column " << i + 1 << "/" << columns_to_reduce.size()
            << " (diameter " << comp_prev.diameter(index) << ")"
            << std::flush << "\r";
            #endif
            
            index_t pivot, column = index;
            
            std::vector<index_t> coboundary;

            do {
            
                #ifdef ASSEMBLE_REDUCTION_COLUMN
                reduction_column.push( column );
                #endif

                simplex_coboundary_enumerator coboundaries(column, dim, n, binomial_coeff);
                while (coboundaries.has_next()) {
                    index_t coface = coboundaries.next();
                    if (comp.diameter(coface) <= threshold)
                        working_coboundary.push(coface);
                }
                               
                pivot = get_pivot(working_coboundary);
                

                //add boundary column at index birth to working_coboundary
                
                //since the boundary column is not stored explicitly,
                //add the boundary of the column at index birth in the reduction matrix instead
                
                //space-efficient variant: add just the boundary column at index birth instead
                //this avoids having to store the reduction matrix
                
                if (pivot != -1) {
                    auto pair = pivot_column_index.find(pivot);
                    
                    if (pair == pivot_column_index.end()) {
                        pivot_column_index.insert(std::make_pair(pivot, index));
                        
                        #ifdef PRINT_PERSISTENCE_PAIRS
                        value_t birth = comp_prev.diameter(index), death = comp.diameter(pivot);
                        if (birth != death)
                            std::cout << "\033[K" << " [" << birth << "," << death << ")" << std::endl << std::flush;
                        #endif

                        break;
                    }
                    
                    column = pair->second;
                }

            } while ( pivot != -1 );
            
            #ifdef PRINT_PERSISTENCE_PAIRS
            if ( pivot == -1 ) {
                value_t birth = comp_prev.diameter(index);
                std::cout << "\033[K" << " [" << birth << ", )" << std::endl << std::flush;
            }
            #endif

        }
    
        std::cout << "\033[K";

}




int main( int argc, char** argv ) {
    
    if( argc < 3 ) print_help_and_exit();

    std::string input_filename = argv[ argc - 2 ];
    std::string output_filename = argv[ argc - 1 ];
    
    index_t dim_max = 1;
    value_t threshold = std::numeric_limits<value_t>::max();

    for( index_t idx = 1; idx < argc - 2; idx++ ) {
        const std::string option = argv[ idx ];
        if( option == "--help" ) {
            print_help_and_exit();
        } else if( option == "--top_dim" ) {
            idx++;
            if( idx >= argc - 2 )
                print_help_and_exit();
            std::string parameter = std::string( argv[ idx ] );
            size_t pos_last_digit;
            dim_max = std::stoll( parameter, &pos_last_digit );
            if( pos_last_digit != parameter.size() )
                print_help_and_exit();
        } else if( option == "--threshold" ) {
            idx++;
            if( idx >= argc - 2 )
                print_help_and_exit();
            std::string parameter = std::string( argv[ idx ] );
            size_t pos_last_digit;
            threshold = std::stof( parameter, &pos_last_digit );
            if( pos_last_digit != parameter.size() )
                print_help_and_exit();
	} else print_help_and_exit();
    }


    std::ifstream input_stream( input_filename.c_str( ), std::ios_base::binary | std::ios_base::in );
    if( input_stream.fail( ) ) {
        std::cerr << "couldn't open file" << input_filename << std::endl;
        exit(-1);
    }
    
    #ifdef FILE_FORMAT_DIPHA
    
    int64_t magic_number;
    input_stream.read( (char*)&magic_number, sizeof( int64_t ) );
    if( magic_number != 8067171840 ) {
        std::cerr << input_filename << " is not a Dipha file (magic number: 8067171840)" << std::endl;
        exit(-1);
    }

   
    int64_t file_type;
    input_stream.read( (char*)&file_type, sizeof( int64_t ) );
    if( file_type != 7 ) {
        std::cerr << input_filename << " is not a Dipha distance matrix (file type: 7)" << std::endl;
        exit(-1);
    }
   
    int64_t n;
    input_stream.read( (char*)&n, sizeof( int64_t ) );
    

    //std::vector<value_t> distances =

    distance_matrix dist;
    dist.distances = std::vector<std::vector<value_t>>(n, std::vector<value_t>(n));
    
    for( int i = 0; i < n; ++i ) {
        input_stream.read( (char*)&dist.distances[i][0], n * sizeof(int64_t) );
    }
    
    //std::cout << dist.distances << std::endl;
    
    #endif


    #ifdef FILE_FORMAT_CSV
    
    std::vector<value_t> distances;
    std::string value_string;
    while(std::getline(input_stream, value_string, ','))
        distances.push_back(std::stod(value_string));
    
//    std::cout << "distances: " << distances << std::endl;
    
    compressed_upper_distance_matrix dist_upper(std::move(distances));
    
    index_t n = dist_upper.size();
    
    #endif
    
    std::cout << "distance matrix with " << n << " points" << std::endl;
    
//    std::vector<std::vector<value_t>> distance_matrix_full(n, std::vector<value_t>(n));
//    for (index_t i = 0; i < n; ++i)
//        for (index_t j = 0; j < n; ++j)
//            distance_matrix_full[i][j] = dist_upper(i, j);
//    std::cout << distance_matrix_full << std::endl;



    compressed_lower_distance_matrix dist(dist_upper);

    std::cout << "distance matrix transformed to lower triangular form" << std::endl;

//    for (index_t i = 0; i < n; ++i)
//        for (index_t j = 0; j < n; ++j)
//            distance_matrix_full[i][j] = dist(i, j);
//    std::cout << distance_matrix_full << std::endl;
    
    
//    std::cout << "distances (lower): " << dist.distances << std::endl;

//    return 0;
    
    assert(dim_max < n - 2);
    
    binomial_coeff_table binomial_coeff(n, dim_max + 2);
    
//    std::vector<index_t> coboundary;
//    get_simplex_coboundary( 1, 2, n, binomial_coeff, std::back_inserter(coboundary) );
//    std::cout << coboundary <<std::endl;
//    
//    coboundary.clear();
//    simplex_coboundary_enumerator coboundaries(1, 2, n, binomial_coeff);
//    while (coboundaries.has_next())
//        coboundary.push_back(coboundaries.next());
//    std::cout << coboundary <<std::endl;
   
    
    std::vector<index_t> columns_to_reduce;
    
    
    for (index_t index = n; index-- > 0; ) {
        columns_to_reduce.push_back(index);
    }

    
    #ifdef PRECOMPUTE_DIAMETERS
    std::vector<value_t> previous_diameters( n , 0 );


    std::vector<value_t>& diameters = dist.distances;
    
//    std::cout << "precomputing 1-simplex diameters" << std::endl;
//    
//    std::vector<value_t> diameters( binomial_coeff( n, 2 ) );
//    
//    std::vector<index_t> edge_vertices(2);
//    for (index_t edge = 0; edge < diameters.size(); ++edge) {
//        #ifdef INDICATE_PROGRESS
//        std::cout << "\033[Kstoring diameter of simplex " << edge + 1 << "/" << diameters.size() << std::flush << "\r";
//        #endif
//        
//        edge_vertices.clear();
//        get_simplex_vertices( edge, 1, n, binomial_coeff, std::back_inserter(edge_vertices) );
//        diameters[edge] = dist(edge_vertices[0], edge_vertices[1]);
//    }
//    #ifdef INDICATE_PROGRESS
//    std::cout << "\033[K";
//    #endif

    #endif

    
//    std::cout << "diameters: " << diameters << std::endl;
    
    
    for (index_t dim = 0; dim < dim_max; ++dim) {
    
        std::unordered_map<index_t, index_t> pivot_column_index;
    
        #ifdef PRECOMPUTE_DIAMETERS
        rips_filtration_diameter_comparator comp(diameters, dim + 1, binomial_coeff);
        rips_filtration_diameter_comparator comp_prev(previous_diameters, dim, binomial_coeff);
        #else
        rips_filtration_comparator<decltype(dist)> comp(dist, dim + 1, binomial_coeff);
        rips_filtration_comparator<decltype(dist)> comp_prev(dist, dim, binomial_coeff);
        #endif

        compute_pairs(
            columns_to_reduce,
            pivot_column_index,
            comp, comp_prev,
            dim, dim_max, n,
            threshold,
            binomial_coeff
        );
        
        
        index_t num_simplices = binomial_coeff(n, dim + 2);
        
        columns_to_reduce.clear();
        
        for (index_t index = 0; index < num_simplices; ++index) {

            if (comp.diameter(index) <= threshold && pivot_column_index.find(index) == pivot_column_index.end()) {
                columns_to_reduce.push_back(index);
            }
        }
        
        std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), comp);

        
        
        #ifdef PRECOMPUTE_DIAMETERS
        previous_diameters = std::move(diameters);
        
        #ifndef PRECOMPUTE_DIAMETERS_IN_TOP_DIMENSION
        if (dim == dim_max - 1)
            break;
        #endif
        
        std::cout << "precomputing " << dim + 2 << "-simplex diameters" << std::endl;
        diameters = get_diameters( dist, dim + 2, previous_diameters, binomial_coeff );
        #endif
        
    }
    
    index_t dim = dim_max;
    
    #ifdef PRECOMPUTE_DIAMETERS
    rips_filtration_diameter_comparator comp_prev(previous_diameters, dim, binomial_coeff);
    #else
    rips_filtration_comparator<decltype(dist)> comp_prev(dist, dim, binomial_coeff);
    #endif
    
    #ifdef PRECOMPUTE_DIAMETERS_IN_TOP_DIMENSION
    rips_filtration_diameter_comparator comp(diameters, dim + 1, binomial_coeff);
    #else
    rips_filtration_comparator<decltype(dist)> comp(dist, dim + 1, binomial_coeff);
    #endif
    
    std::unordered_map<index_t, index_t> pivot_column_index;

    compute_pairs(
        columns_to_reduce,
        pivot_column_index,
        comp, comp_prev,
        dim, dim_max, n,
        threshold,
        binomial_coeff
    );
    
    
}
