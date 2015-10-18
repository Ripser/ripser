#include <vector>
#include <iostream>
#include <cassert>
#include "prettyprint.hpp"
//#include <boost/numeric/ublas/symmetric.hpp>



class binomial_coeff_table {
    std::vector<std::vector<int> > B;
    int n_max, k_max;
    
public:
    binomial_coeff_table(int n, int k) {
        n_max = n;
        k_max = k;
    
        B.resize(n + 1);
        for( int i = 0; i <= n; i++ ) {
            B[i].resize(k + 1);
            for ( int j = 0; j <= std::min(i, k); j++ )
            {
                if (j == 0 || j == i)
                    B[i][j] = 1;

                else
                    B[i][j] = B[i-1][j-1] + B[i-1][j];
            }
        }
    }
    
    int operator()(int n, int k) const {
//        std::cout << "B(" << n << "," << k << ")\n";

        return B[n][k];
    }
};



template<typename OutputIterator>
OutputIterator get_simplex_vertices( int idx, int dim, int n, const binomial_coeff_table& binomial_coeff, OutputIterator out  )
{
    --n;
    
    for( int k = dim + 1; k > 0; --k ) {
        while( binomial_coeff( n , k ) > idx ) {
            --n;
        }
        *out++ = n;
        
        idx -= binomial_coeff( n , k );

    }

    return out;
}

template <typename DistanceMatrix>
class rips_filtration_comparator {
private:
    std::vector<int> vertices;
    
    const binomial_coeff_table& binomial_coeff;
    
public:
    const DistanceMatrix& dist;
    
    const int dim;

    rips_filtration_comparator(const DistanceMatrix& _dist, const int _dim, const binomial_coeff_table& _binomial_coeff): dist(_dist), dim(_dim), vertices(_dim + 1),
    binomial_coeff(_binomial_coeff) {};

    bool operator()(const int a, const int b)
    {
        assert(a < binomial_coeff(dist.size(), dim + 1));
        assert(b < binomial_coeff(dist.size(), dim + 1));
        
        decltype(dist(0,0)) a_diam = 0, b_diam = 0;

    
        //vertices.clear(); //std::back_inserter(vertices)
        get_simplex_vertices(a, dim, dist.size(), binomial_coeff, vertices.begin() );
        
        for (int i = 0; i <= dim; ++i)
            for (int j = i + 1; j <= dim; ++j) {
                auto d = dist(vertices[i], vertices[j]);
                a_diam = std::max(a_diam, dist(vertices[i], vertices[j]));
//                std::cout << " i = " << i << " j = " << j << " d = " << d << std::endl;
//                std::cout << " a_vertices[i] = " << vertices[i] << " a_vertices[j] = " << vertices[j] << std::endl;
//                std::cout << " a_diam = " << a_diam << std::endl;
            }
        
//        std::cout << "a_diam = " << a_diam << std::endl;
        
        //vertices.clear(); //std::back_inserter(vertices)
        get_simplex_vertices(b, dim, dist.size(), binomial_coeff, vertices.begin() );

        for (int i = 0; i <= dim; ++i)
            for (int j = i + 1; j <= dim; ++j) {
                b_diam = std::max(b_diam, dist(vertices[i], vertices[j]));
                if (a_diam < b_diam)
                    return true;
            }
        
//        std::cout << "b_diam = " << b_diam << std::endl;
       
        return (a_diam == b_diam && a <= b);
    }
    
};


template<typename OutputIterator>
void get_simplex_coboundary( int idx, int dim, int n, const binomial_coeff_table& binomial_coeff, OutputIterator coboundary )
{

    --n;
    
    int modified_idx = idx;
    
    for( int k = dim + 1; k >= 0 && n >= 0; --k ) {
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



class distance_matrix  {
public:
    typedef double value_type;
    std::vector<std::vector<double> > distances;
    double operator()(const int a, const int b) const {
        return distances[a][b];
    }
    
    size_t size() const {
        return distances.size();
    }

};


int main() {

    distance_matrix dist;
    dist.distances = {
        {0,1,2,3},
        {1,0,1,2},
        {2,1,0,1},
        {3,2,1,0}
    };
    dist.distances = {
        {0,1,1,1,1},
        {1,0,1,1,1},
        {1,1,0,1,1},
        {1,1,1,0,1},
        {1,1,1,1,0}
    };
    
    binomial_coeff_table binomial_coeff(dist.size(), 4);
    
    
    std::cout << dist (0,1) << std::endl;
 
    rips_filtration_comparator<distance_matrix> comp1(dist, 1, binomial_coeff);
    rips_filtration_comparator<distance_matrix> comp2(dist, 2, binomial_coeff);
    
    std::cout << (comp1(0,1) ? "0<1" : "0>=1") << std::endl;
    std::cout << (comp1(1,0) ? "1<0" : "1>=0") << std::endl;
    
    std::cout << (comp1(0,2) ? "0<2" : "0>=2") << std::endl;
    std::cout << (comp1(1,2) ? "1<2" : "1>=2") << std::endl;
    
    std::vector<int> edges = {0,1,2,3,4,5};
    
    std::sort(edges.begin(), edges.end(), comp1);
    
    std::cout << "sorted edges: " << edges << std::endl;
   
   
    std::vector<int> triangles = {0,1,2,3};
    
    std::sort(triangles.begin(), triangles.end(), comp2);
    
    std::cout << "sorted triangles: " << triangles << std::endl;
    
    
    int dim = 1;
    int simplex_index = 2;
    
    
    std::vector<int> vertices;
    
    get_simplex_vertices( simplex_index, dim, dist.size(), binomial_coeff, std::back_inserter(vertices) );

    
    std::cout << "coboundary of simplex " << vertices << ":" << std::endl;
    
    std::vector<int> coboundary;
    get_simplex_coboundary( simplex_index, dim, dist.size(), binomial_coeff, std::back_inserter(coboundary) );
    
    
    for (int coboundary_simplex_index: coboundary) {
        std::vector<int> vertices;
    
        get_simplex_vertices( coboundary_simplex_index, dim + 1, dist.size(), binomial_coeff, std::back_inserter(vertices) );
        std::cout << " " << vertices << std::endl;

    }

}