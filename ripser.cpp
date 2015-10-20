#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <queue>
#include <unordered_map>
#include "prettyprint.hpp"
//#include <boost/numeric/ublas/symmetric.hpp>



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
    std::vector<long> vertices;
    
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
    
    dist_t diameter(const long index) {
        dist_t diam = 0;
        get_simplex_vertices(index, dim, dist.size(), binomial_coeff, vertices.begin() );
        
        for (long i = 0; i <= dim; ++i)
            for (long j = i + 1; j <= dim; ++j) {
                auto d = dist(vertices[i], vertices[j]);
                diam = std::max(diam, dist(vertices[i], vertices[j]));
            }
        return diam;
    }

    bool operator()(const long a, const long b)
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
    
    std::vector<double> diameters(binomial_coeff(n, dim + 1));
    
    std::vector<long> coboundary;

    for (long simplex = 0; simplex < n; ++simplex) {
        coboundary.clear();
        
        get_simplex_coboundary( simplex, dim - 1, n, binomial_coeff, std::back_inserter(coboundary) );
        for (long coface: coboundary) {
            diameters[coface] = std::max( diameters[coface], previous_diameters[simplex]);
        }
    }
    
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
    
    double diameter(const long a) {
        assert(a < diameters.size());
        return diameters[a];
    }

    bool operator()(const long a, const long b)
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

int main( int argc, char** argv ) {
    
    std::cout << "sizeof(long) = " << sizeof(long) << std::endl;
    
    if( argc < 3 ) print_help_and_exit();

    std::string input_filename = argv[ argc - 2 ];
    std::string output_filename = argv[ argc - 1 ];
    
    long dim_max = 2;
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
    
    std::vector<double> distances;
    std::string value_string;
    while(std::getline(input_stream, value_string, ','))
        distances.push_back(std::stod(value_string));
    
//    std::cout << "distances: " << distances << std::endl;
    
    compressed_distance_matrix dist(std::move(distances));
    
    long n = dist.size();
    
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

    
//    dist.distances = {
//        {0,1,3,4,3,1},
//        {1,0,1,3,4,3},
//        {3,1,0,1,3,4},
//        {4,3,1,0,1,3},
//        {3,4,3,1,0,1},
//        {1,3,4,3,1,0}
//    };
//    n = dist.size();

    
    assert(dim_max < n - 1);
    
    binomial_coeff_table binomial_coeff(n, dim_max + 2);
    
    
//    std::cout << dist (0,1) << std::endl;
// 
//    rips_filtration_comparator<distance_matrix> comp1(dist, 1, binomial_coeff);
//    rips_filtration_comparator<distance_matrix> comp2(dist, 2, binomial_coeff);
//    
//    std::cout << (comp1(0,1) ? "0<1" : "0>=1") << std::endl;
//    std::cout << (comp1(1,0) ? "1<0" : "1>=0") << std::endl;
//    
//    std::cout << (comp1(0,2) ? "0<2" : "0>=2") << std::endl;
//    std::cout << (comp1(1,2) ? "1<2" : "1>=2") << std::endl;
//    
//    std::vector<long> edges = {0,1,2,3,4,5};
//    
//    std::sort(edges.begin(), edges.end(), comp1);
//    
//    std::cout << "sorted edges: " << edges << std::endl;
//   
//   
//    std::vector<long> triangles = {0,1,2,3};
//    
//    std::sort(triangles.begin(), triangles.end(), comp2);
//    
//    std::cout << "sorted triangles: " << triangles << std::endl;
//    
//    
//    long dim = 1;
//    long simplex_index = 2;
//    
//    double threshold = 7;
//    
//    std::vector<long> vertices;
//    
//    get_simplex_vertices( simplex_index, dim, n, binomial_coeff, std::back_inserter(vertices) );
//
//    
//    std::cout << "coboundary of simplex " << vertices << ":" << std::endl;
//    
//    std::vector<long> coboundary;
//    get_simplex_coboundary( simplex_index, dim, n, binomial_coeff, std::back_inserter(coboundary) );
//    
//    
//    for (long coboundary_simplex_index: coboundary) {
//        std::vector<long> vertices;
//    
//        get_simplex_vertices( coboundary_simplex_index, dim + 1, n, binomial_coeff, std::back_inserter(vertices) );
//        std::cout << " " << coboundary_simplex_index << " " << vertices << " (" << comp1.diameter(coboundary_simplex_index) << ")" << std::endl;
//
//    }
//    
//    
//    compressed_sparse_matrix<long> csm;
//    
//    csm.append(std::vector<long>({1,2,3}));
//    
//    csm.append(std::vector<long>({5,6,7,8}));
//
//    csm.append(std::vector<long>({10,11}));
//
//    csm.append(std::vector<long>());
//
//    csm.append(std::vector<long>({2}));
//    
//    std::cout << "compressed sparse matrix: " << std::endl;
//
//    for (long i = 0; i < csm.size(); ++i) {
//        std::cout << " " << std::vector<long>(csm.cbegin(i), csm.cend(i)) << std::endl;
//    }
//    
//    std::cout << "bounds: " << csm.bounds << std::endl;
//
//    std::cout << "entries: " << csm.entries << std::endl;
//    
//    
//    std::priority_queue<long, std::vector<long>, rips_filtration_comparator<distance_matrix> > queue(comp1);
//    
//    for (long e: coboundary) queue.push(e);
//    
//    std::cout << "pivot of coboundary: " << queue.top() << std::endl;
//    
//    std::cout << (comp1(3,6) ? "3<6" : "3>=6") << std::endl;
//    std::cout << (comp1(0,6) ? "0<6" : "0>=6") << std::endl;
//    
//    
//    
//    std::vector<long> columns_to_reduce = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
//
//    std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), comp1);
//    
//    std::cout << "sorted 1-simplex columns to reduce: " << columns_to_reduce << std::endl;
//    
//    for (long index: columns_to_reduce) {
//        std::vector<long> vertices;
//
//        get_simplex_vertices( index, 1, n, binomial_coeff, std::back_inserter(vertices) );
//        std::cout << " " << index << " " << vertices << " (" << comp1.diameter(index) << ")" << std::endl;
//    }
//    
//    
//    std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), comp2);
//    
//    std::cout << "sorted 2-simplex columns to reduce: " << columns_to_reduce << std::endl;
//    
//    for (long index: columns_to_reduce) {
//        std::vector<long> vertices;
//
//        get_simplex_vertices( index, 2, n, binomial_coeff, std::back_inserter(vertices) );
//        std::cout << " " << index << " " << vertices << " (" << comp2.diameter(index) << ")" << std::endl;
//    }
//
//
//
//
//

    std::vector<long> columns_to_reduce;
    
    std::vector<long> coboundary;
    
    for (long index = n - 1; index >= 0; --index) {
        columns_to_reduce.push_back(index);
    }

    std::vector<double> previous_diameters( n );

    std::cout << "computing edge diameters" << std::endl;
    
    std::vector<double> diameters( binomial_coeff( n, 2 ) );
    
//    rips_filtration_comparator<decltype(dist)> comp1(dist, 1, binomial_coeff);

    std::vector<long> edge_vertices(2);
    for (long edge = 0; edge < diameters.size(); ++edge) {
//        diameters[edge] = comp1.diameter(edge);
        
        edge_vertices.clear();
        get_simplex_vertices( edge, 1, n, binomial_coeff, std::back_inserter(edge_vertices) );
        diameters[edge] = dist(edge_vertices[0], edge_vertices[1]);
    }
    
    for (long dim = 0; dim < dim_max; ++dim) {
        

        //compressed_sparse_matrix<long> reduction_matrix;
        
        //rips_filtration_comparator<decltype(dist)> comp(dist, dim + 1, binomial_coeff);
        
        rips_filtration_diameter_comparator comp(diameters, dim + 1, binomial_coeff);
        
        std::unordered_map<long, long> pivot_column_index;
        
        std::cout << "reduce columns in dim " << dim << std::endl;
        
//        std::cout << "reduce columns in dim " << dim << ": " << columns_to_reduce << std::endl;
//        std::cout << " rows in dim " << dim + 1 << ": " << columns_to_reduce << std::endl;
        
//        std::vector<long> rows(binomial_coeff(n, dim + 2));
//        for (long simplex_index = 0; simplex_index < binomial_coeff(n, dim + 2); ++simplex_index) {
//            rows[simplex_index] = simplex_index;
//        }
//        std::sort(rows.begin(), rows.end(), comp);
//        
//        for (long simplex_index: rows) {
//            std::vector<long> vertices;
//            
//            get_simplex_vertices( simplex_index, dim + 1, n, binomial_coeff, std::back_inserter(vertices) );
//            std::cout << " " << simplex_index << " " << vertices << " (" << comp.diameter(simplex_index) << ")" << std::endl;
//
//        }
        
        for (long i = 0; i < columns_to_reduce.size(); ++i) {
            long index = columns_to_reduce[i];
            
            std::priority_queue<long, std::vector<long>, decltype(comp) > reduction_column(comp);
            
            std::priority_queue<long, std::vector<long>, decltype(comp) > working_coboundary(comp);
            
//            std::cout << "reduce column " << index << " (" << i + 1 << "/" << columns_to_reduce.size() << ")" << std::endl;

        
            long pivot, column = index;
            
            do {
            
                reduction_column.push( column );

                coboundary.clear();
                get_simplex_coboundary( column, dim, n, binomial_coeff, std::back_inserter(coboundary) );
                
                std::vector<long> sorted_coboundary = coboundary; std::sort(sorted_coboundary.begin(), sorted_coboundary.end(), comp);
//                std::cout << "add " << sorted_coboundary << " to working col" << std::endl;
//                for (long coboundary_simplex_index: coboundary) {
//                    std::vector<long> vertices;
//                
//                    get_simplex_vertices( coboundary_simplex_index, dim + 1, n, binomial_coeff, std::back_inserter(vertices) );
//                    std::cout << " " << coboundary_simplex_index << " " << vertices << " (" << comp.diameter(coboundary_simplex_index) << ")" << std::endl;
//                }
                
                for (long e: coboundary) if (comp.diameter(e) <= threshold) working_coboundary.push(e); // std::cout << "push " << e << std::endl;}
                
//            std::cout << "=" << std::endl;
//                auto working_coboundary_copy = working_coboundary;
//                while (!working_coboundary_copy.empty()) {
//                    std::cout << " " << working_coboundary_copy.top() << std::endl;
//                    working_coboundary_copy.pop();
//                }
//            std::cout << "=" << std::endl;

                
                pivot = get_pivot(working_coboundary);
                

                //add boundary column at index birth to working_coboundary
                
                //since the boundary column is not stored explicitly,
                //add the boundary of the column at index birth in the reduction matrix instead
                
                //space-efficient variant: add just the boundary column at index birth instead
                //this avoids having to store the reduction matrix
                
                if (pivot != -1) {
//                    std::cout << "pivot: " << pivot << std::endl;
                    auto pair = pivot_column_index.find(pivot);
                    
                    if (pair == pivot_column_index.end()) {
                        pivot_column_index.insert(std::make_pair(pivot, index));
                        break;
                    }
                    
                    column = pair->second;
                }
                
                /*
                coboundary.clear();
                get_simplex_coboundary( birth, dim, n, binomial_coeff, std::back_inserter(coboundary) );
                for (long e: coboundary) if (comp2.diameter(e) <= threshold) working_coboundary.push(e);

                //space-efficient variant: drop this part (and the reduction_matrix)
                
                for (long col = reduction_matrix.cbegin()) {
                    coboundary.clear();
                    get_simplex_coboundary( col, dim, n, binomial_coeff, std::back_inserter(coboundary) );
                    for (long e: coboundary) if (comp2.diameter(e) <= threshold) working_coboundary.push(e);
                }
                 */
            } while ( pivot != -1 );
            
//            std::cout << "size of working column heap: " << working_coboundary.size() << std::endl;
//            
//            std::vector<long> cochain = move_to_column_vector( reduction_column );
//            std::cout << "reduction cochain: " << cochain << std::endl;
//
//            if ( pivot != -1 ) {
//                std::vector<long> coboundary = move_to_column_vector( working_coboundary );
//                std::cout << "reduced coboundary: " << coboundary << std::endl;
//            }
//
//            std::cout << "fill-in: " << cochain.size()-1 << std::endl;



        }
        
        std::cout << "dimension " << dim << " pairs:" << std::endl;
//        std::cout << pivot_column_index << std::endl;
        
        
//        rips_filtration_comparator<decltype(dist)> comp_prev(dist, dim, binomial_coeff);
        rips_filtration_diameter_comparator comp_prev(previous_diameters, dim, binomial_coeff);

        
        for (std::pair<long,long> pair: pivot_column_index) {
            double birth = comp_prev.diameter(pair.second), death = comp.diameter(pair.first);
//            std::cout << vertices_of_simplex(pair.second, dim, n, binomial_coeff) << "," <<
//                        vertices_of_simplex(pair.first, dim + 1, n, binomial_coeff) << std::endl;
            if (birth != death)
                std::cout << " [" << birth << "," << death << ")" << std::endl;
        }

        if (dim == dim_max - 1)
            break;
        
        
        long num_simplices = binomial_coeff(n, dim + 2);
        
        columns_to_reduce.clear();
        
//        std::cout << "columns to reduce in dim " << dim + 1 << " (" << num_simplices << " total)" << std::endl;
        
        for (long index = 0; index < num_simplices; ++index) {
        
//            if (comp.diameter(index) > threshold) {
//                std::cout << " " << vertices_of_simplex(index, dim + 1, n, binomial_coeff) << ": " << comp.diameter(index) << " above threshold" << std::endl;
//            } else if (pivot_column_index.find(index) != pivot_column_index.end()) {
//                std::cout << " " << vertices_of_simplex(index, dim + 1, n, binomial_coeff) << " appears in pair" << std::endl;
//            }
            
            if (comp.diameter(index) <= threshold && pivot_column_index.find(index) == pivot_column_index.end()) {
                columns_to_reduce.push_back(index);
            }
        }
        std::sort(columns_to_reduce.begin(), columns_to_reduce.end(), comp);

//        std::cout << "sorted " << dim + 1 << "-columns to reduce: " << columns_to_reduce << std::endl;
//        
//        for (long index: columns_to_reduce) {
//            std::vector<long> vertices;
//
//            get_simplex_vertices( index, dim, n, binomial_coeff, std::back_inserter(vertices) );
//            std::cout << " " << index << " " << vertices << " (" << comp.diameter(index) << ")" << std::endl;
//        }
//        
//        std::cout << std::endl;
    
        
        previous_diameters = diameters;
        diameters = get_diameters( dist, dim + 2, previous_diameters, binomial_coeff);
        


    }
    
}