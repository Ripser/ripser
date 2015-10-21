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


#define PRECOMPUTE_DIAMETERS

#define PRECOMPUTE_DIAMETERS_IN_TOP_DIMENSION


//#define ASSEMBLE_REDUCTION_COLUMN

//#define FILE_FORMAT_DIPHA
#define FILE_FORMAT_CSV


class binomial_coeff_table {
    std::vector<std::vector<long> > B;
    long n_max, k_max;
    
public:
    binomial_coeff_table(long n, long k) {
        n_max = n;
        k_max = k;
    
        B.resize(n + 1);
        for( long i = 0; i <= n; i++ ) {
            B[i].resize(k + 1);
            for ( long j = 0; j <= std::min(i, k); j++ )
            {
                if (j == 0 || j == i)
                    B[i][j] = 1;

                else
                    B[i][j] = B[i-1][j-1] + B[i-1][j];
            }
        }
    }
    
    long operator()(long n, long k) const {
//        std::cout << "B(" << n << "," << k << ")\n";
        assert(n <= n_max);
        assert(k <= k_max);
        return B[n][k];
    }
};



template<typename OutputIterator>
OutputIterator get_simplex_vertices( long idx, long dim, long n, const binomial_coeff_table& binomial_coeff, OutputIterator out  )
{
    --n;
    
    for( long k = dim + 1; k > 0; --k ) {
        while( binomial_coeff( n , k ) > idx ) {
            --n;
        }
        *out++ = n;
        
        idx -= binomial_coeff( n , k );

    }

    return out;
}

std::vector<long> vertices_of_simplex(long simplex_index, long dim, long n, const binomial_coeff_table& binomial_coeff) {
    std::vector<long> vertices;
    get_simplex_vertices( simplex_index, dim, n, binomial_coeff, std::back_inserter(vertices) );
    return vertices;
}

template <typename DistanceMatrix>
class rips_filtration_comparator {
public:
    const DistanceMatrix& dist;
    
    const long dim;

private:
    mutable std::vector<long> vertices;
    
    typedef decltype(dist(0,0)) dist_t;
    
    bool reverse;
    
    const binomial_coeff_table& binomial_coeff;
    
public:
    rips_filtration_comparator(
    const DistanceMatrix& _dist,
    const long _dim,
    const binomial_coeff_table& _binomial_coeff
//    ,
//    bool _reverse = true
    ):
    dist(_dist), dim(_dim), vertices(_dim + 1),
    binomial_coeff(_binomial_coeff)
//    ,
//    reverse(_reverse)
    {};
    
    dist_t diameter(const long index) const {
        dist_t diam = 0;
        get_simplex_vertices(index, dim, dist.size(), binomial_coeff, vertices.begin() );
        
        for (long i = 0; i <= dim; ++i)
            for (long j = i + 1; j <= dim; ++j) {
                diam = std::max(diam, dist(vertices[i], vertices[j]));
            }
        return diam;
    }

    bool operator()(const long a, const long b) const
    {
        assert(a < binomial_coeff(dist.size(), dim + 1));
        assert(b < binomial_coeff(dist.size(), dim + 1));
        
        dist_t a_diam = 0, b_diam = 0;

    
        b_diam = diameter(b);
        
        get_simplex_vertices(a, dim, dist.size(), binomial_coeff, vertices.begin() );
        for (long i = 0; i <= dim; ++i)
            for (long j = i + 1; j <= dim; ++j) {
                a_diam = std::max(a_diam, dist(vertices[i], vertices[j]));
                if (a_diam > b_diam)
//                   (((a_diam < b_diam) && !reverse) ||
//                    ((a_diam > b_diam) && reverse))
                    return true;
            }
       
        return (a_diam == b_diam) && (a > b);
//        return (a_diam == b_diam) && (((a < b) && !reverse) || ((a > b) && reverse));
    }
    
};


template<typename OutputIterator>
void get_simplex_coboundary( long idx, long dim, long n, const binomial_coeff_table& binomial_coeff, OutputIterator coboundary )
{

    --n;
    
    long modified_idx = idx;
    
    for( long k = dim + 1; k >= 0 && n >= 0; --k ) {
        while( binomial_coeff( n , k ) > idx ) {
//            std::cout << "binomial_coeff(" << n << ", " << k << ") = " << binomial_coeff( n , k ) << " > " << idx << std::endl;
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

template <typename DistanceMatrix>
std::vector<double> get_diameters (
    const DistanceMatrix& dist,
    const long dim,
    const std::vector<double>& previous_diameters,
    const binomial_coeff_table& binomial_coeff
    )
{
    long n = dist.size();
    
    std::vector<double> diameters(binomial_coeff(n, dim + 1), 0);
    
    std::vector<long> coboundary;

    for (long simplex = 0; simplex < previous_diameters.size(); ++simplex) {
        coboundary.clear();
        
        std::cout << "\033[Kpropagating diameter of simplex " << simplex + 1 << "/" << previous_diameters.size() << std::flush << "\r";

        
        get_simplex_coboundary( simplex, dim - 1, n, binomial_coeff, std::back_inserter(coboundary) );
        for (long coface: coboundary) {
            diameters[coface] = std::max( diameters[coface], previous_diameters[simplex]);
        }
    }
    
    std::cout << "\033[K";
    
//    rips_filtration_comparator<decltype(dist)> comp(dist, dim, binomial_coeff);
//    for (long simplex = 0; simplex < diameters.size(); ++simplex) {
//        assert(diameters[simplex] == comp.diameter(simplex));
//    }

    return diameters;
}


class rips_filtration_diameter_comparator {
private:
    const std::vector<double>& diameters;
    
    const long dim;

public:
    std::vector<long> vertices;
    
    typedef double dist_t;
    
    
    const binomial_coeff_table& binomial_coeff;
    
public:
    rips_filtration_diameter_comparator(
        const std::vector<double>& _diameters,
        const long _dim,
        const binomial_coeff_table& _binomial_coeff
        ):
        diameters(_diameters), dim(_dim),
        binomial_coeff(_binomial_coeff)
    {};
    
    double diameter(const long a) const {
        assert(a < diameters.size());
        return diameters[a];
    }

    bool operator()(const long a, const long b) const
    {
        assert(a < diameters.size());
        assert(b < diameters.size());
        
        dist_t a_diam = diameters[a], b_diam = diameters[b];

        return ((a_diam > b_diam) || ((a_diam == b_diam) && (a > b)));
    }
    
};


class distance_matrix  {
public:
    typedef double value_type;
    std::vector<std::vector<double> > distances;
    double operator()(const long a, const long b) const {
        return distances[a][b];
    }
    
    size_t size() const {
        return distances.size();
    }

};


class compressed_distance_matrix  {
public:
    typedef double value_type;
    std::vector<double> distances;
    std::vector<double*> rows;
    
    long n;
    
    compressed_distance_matrix(std::vector<double>&& row_vector) noexcept {
        n = (1 + std::sqrt(1+ 8 * row_vector.size())) / 2;
        
        distances = std::move(row_vector);
        
        rows.resize(n);
        
        double* pointer = &distances[0] - 1;
        for (long i = 0; i < n - 1; ++i) {
            rows[i] = pointer;
            pointer += n - i - 2;
        }
    }
    
    double operator()(const long a, const long b) const {
        if (a < b)
            return rows[a][b];
        else if (a > b)
            return rows[b][a];
        else
            return 0;
    }
    
    size_t size() const {
        return n;
    }

};

template <typename ValueType>
class compressed_sparse_matrix  {
public:
    std::vector<size_t> bounds;
    std::vector<ValueType> entries;
    
    
    size_t size() const {
        return bounds.size();
    }
    
    typename std::vector<ValueType>::const_iterator cbegin(size_t index) {
        assert(index < size());
        return index == 0 ? entries.cbegin() : entries.cbegin() + bounds[index - 1];
    }

    typename std::vector<ValueType>::const_iterator cend(size_t index) {
        assert(index < size());
        return entries.cbegin() + bounds[index];
    }
    
    template <typename Iterator>
    void append(Iterator begin, Iterator end) {
        for (Iterator it = begin; it != end; ++it) {
            entries.push_back(*it);
        }
        bounds.push_back(entries.size());
    }

    
    template <typename Collection>
    void append(const Collection collection) {
        append(collection.cbegin(), collection.cend());
    }

};

template <typename Heap>
long pop_pivot(Heap& column)
{
    if( column.empty() )
        return -1;
    else {
        long max_element = column.top();
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
long get_pivot(Heap& column)
{
    long max_element = pop_pivot(column);
    if( max_element != -1 )
        column.push( max_element );
    return max_element;
}

template <typename Heap>
std::vector<long> move_to_column_vector(Heap& column)
{
    std::vector<long> temp_col;
    long pivot = pop_pivot( column );
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
    std::vector<long>& columns_to_reduce,
    std::unordered_map<long, long>& pivot_column_index,
    const ComparatorCofaces& comp, const Comparator& comp_prev,
    long dim, long dim_max, long n,
    double threshold,
    const binomial_coeff_table& binomial_coeff
) {
    
        std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
        
        for (long i = 0; i < columns_to_reduce.size(); ++i) {
            long index = columns_to_reduce[i];
            
            #ifdef ASSEMBLE_REDUCTION_COLUMN
            std::priority_queue<long, std::vector<long>, decltype(comp) > reduction_column(comp);
            #endif
            
            std::priority_queue<long, std::vector<long>, decltype(comp) > working_coboundary(comp);
            
            std::cout << "\033[K" << "reducing column " << i + 1 << "/" << columns_to_reduce.size()
            << " (diameter " << comp_prev.diameter(index) << ")"
            << std::flush << "\r";
            
            long pivot, column = index;
            
            std::vector<long> coboundary;

            do {
            
                #ifdef ASSEMBLE_REDUCTION_COLUMN
                reduction_column.push( column );
                #endif

                coboundary.clear();
                get_simplex_coboundary( column, dim, n, binomial_coeff, std::back_inserter(coboundary) );
                
                for (long e: coboundary) if (comp.diameter(e) <= threshold) working_coboundary.push(e); // std::cout << "push " << e << std::endl;}
                               
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
                        
                        double birth = comp_prev.diameter(index), death = comp.diameter(pivot);
                        if (birth != death)
                            std::cout << "\033[K" << " [" << birth << "," << death << ")" << std::endl << std::flush;

                        break;
                    }
                    
                    column = pair->second;
                }

            } while ( pivot != -1 );
            
            
            if ( pivot == -1 ) {
                double birth = comp_prev.diameter(index);
                std::cout << "\033[K" << " [" << birth << ", )" << std::endl << std::flush;
            }

        }
    
        std::cout << "\033[K";

}




int main( int argc, char** argv ) {
    
    if( argc < 3 ) print_help_and_exit();

    std::string input_filename = argv[ argc - 2 ];
    std::string output_filename = argv[ argc - 1 ];
    
    long dim_max = 1;
    double threshold = std::numeric_limits<double>::max();

    for( long idx = 1; idx < argc - 2; idx++ ) {
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
            threshold = std::stod( parameter, &pos_last_digit );
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
    

    //std::vector<double> distances =

    distance_matrix dist;
    dist.distances = std::vector<std::vector<double>>(n, std::vector<double>(n));
    
    for( int i = 0; i < n; ++i ) {
        input_stream.read( (char*)&dist.distances[i][0], n * sizeof(int64_t) );
    }
    
    //std::cout << dist.distances << std::endl;
    
    #endif


    #ifdef FILE_FORMAT_CSV
    
    std::vector<double> distances;
    std::string value_string;
    while(std::getline(input_stream, value_string, ','))
        distances.push_back(std::stod(value_string));
    
//    std::cout << "distances: " << distances << std::endl;
    
    compressed_distance_matrix dist(std::move(distances));
    
    long n = dist.size();
    
    #endif
    
    std::cout << "distance matrix with " << n << " points" << std::endl;
    


//    std::vector<std::vector<double>> distance_matrix_full(n, std::vector<double>(n));
//    for (long i = 0; i < n; ++i)
//        for (long j = 0; j < n; ++j)
//            distance_matrix_full[i][j] = dist(i, j);
//    
//    std::cout << distance_matrix_full << std::endl;
//    
//    std::cout << "distances: " << distances << std::endl;
//
//    return 0;

    
    assert(dim_max < n - 2);
    
    binomial_coeff_table binomial_coeff(n, dim_max + 2);
    

    std::vector<long> columns_to_reduce;
    
    
    for (long index = n - 1; index >= 0; --index) {
        columns_to_reduce.push_back(index);
    }

    #ifdef PRECOMPUTE_DIAMETERS
    std::vector<double> previous_diameters( n , 0 );

    std::cout << "precomputing 1-simplex diameters" << std::endl;
    
    std::vector<double> diameters( binomial_coeff( n, 2 ) );
    
//    rips_filtration_comparator<decltype(dist)> comp1(dist, 1, binomial_coeff);

    std::vector<long> edge_vertices(2);
    for (long edge = 0; edge < diameters.size(); ++edge) {
//        diameters[edge] = comp1.diameter(edge);
        
        std::cout << "\033[Kstoring diameter of simplex " << edge + 1 << "/" << diameters.size() << std::flush << "\r";

        edge_vertices.clear();
        get_simplex_vertices( edge, 1, n, binomial_coeff, std::back_inserter(edge_vertices) );
        diameters[edge] = dist(edge_vertices[0], edge_vertices[1]);
    }
    std::cout << "\033[K";
    #endif

    
    
    
    for (long dim = 0; dim < dim_max; ++dim) {
    
        std::unordered_map<long, long> pivot_column_index;
    
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
        
        
        long num_simplices = binomial_coeff(n, dim + 2);
        
        columns_to_reduce.clear();
        
        for (long index = 0; index < num_simplices; ++index) {

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
    
    long dim = dim_max;
    
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
    
    std::unordered_map<long, long> pivot_column_index;

    compute_pairs(
        columns_to_reduce,
        pivot_column_index,
        comp, comp_prev,
        dim, dim_max, n,
        threshold,
        binomial_coeff
    );
    
    
}