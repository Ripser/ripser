/*

 Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes

 MIT License

 Copyright (c) 2015–2019 Ulrich Bauer

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 You are under no obligation whatsoever to provide any bug fixes, patches, or
 upgrades to the features, functionality or performance of the source code
 ("Enhancements") to anyone; however, if you choose to make your Enhancements
 available either publicly, or directly to the author of this software, without
 imposing a separate written license agreement for such Enhancements, then you
 hereby grant the following license: a non-exclusive, royalty-free perpetual
 license to install, use, modify, prepare derivative works, incorporate into
 other computer software, distribute, and sublicense such enhancements or
 derivative works thereof, in binary and source code form.

*/

//#define USE_COEFFICIENTS

#define INDICATE_PROGRESS
//#define PRINT_PERSISTENCE_PAIRS

//#define USE_SERIAL
//#define USE_SHUFFLED_SERIAL

#if defined(USE_SERIAL) || defined(USE_SHUFFLED_SERIAL)
#define USING_SERIAL
#define USE_SERIAL_ATOMIC_REF
#endif

#define USE_CONCURRENT_PIVOTS

//#define USE_GOOGLE_HASHMAP

//#define USE_TBB_HASHMAP

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>

#include <thread>
#include <future>
#include <atomic_ref.hpp>
#include <reclamation.hpp>

#include <counting_iterator.hpp>

#ifndef USING_SERIAL
#include <boost/sort/parallel/sort.hpp>
#endif

#ifdef USE_SHUFFLED_SERIAL
#include <random>
#endif

#ifdef USE_PARALLEL_STL
#include <execution>
#endif

#ifdef USE_TBB
#include <tbb/tbb.h>
#include <tbb/parallel_sort.h>
#endif

#ifdef USE_GOOGLE_HASHMAP
#include <sparsehash/dense_hash_map>
template <class Key, class T, class H, class E>
class hash_map : public google::dense_hash_map<Key, T, H, E> {
public:
	explicit hash_map() : google::dense_hash_map<Key, T, H, E>() { this->set_empty_key(-1); }
	inline void reserve(size_t hint) { this->resize(hint); }
};
#elif defined(USE_CONCURRENT_PIVOTS)
#include <trivial_concurrent_hash_map.hpp>
template <class Key, class T, class H, class E>
class hash_map : public mrzv::trivial_concurrent_hash_map<Key, T, H, E> {};
#elif defined(USE_TBB_HASHMAP)
#include <tbb/concurrent_unordered_map.h>
template <class Key, class T, class H, class E>
class hash_map : public tbb::concurrent_unordered_map<Key, T, H, E>
{
    public:
        using Parent = tbb::concurrent_unordered_map<Key, T, H, E>;
        using iterator = typename Parent::iterator;

        Key key(iterator it) const { return it->first; }
        T   value(iterator it) const { return it->second; }

        bool update(iterator it, T& expected, T desired) { it->second = desired; return true; }

        template<class F>
        void foreach(const F& f) const  { for(auto& x : (*this)) f(x); }

        void reserve(size_t hint) {}
};
#else
template <class Key, class T, class H, class E>
class hash_map : public std::unordered_map<Key, T, H, E>
{
    public:
        using Parent = std::unordered_map<Key,T,H,E>;
        using iterator = typename Parent::iterator;

        Key key(iterator it) const { return it->first; }
        T   value(iterator it) const { return it->second; }

        bool update(iterator it, T& expected, T desired) { it->second = desired; return true; }

        template<class F>
        void foreach(const F& f) const  { for(auto& x : (*this)) f(x); }
};
#endif

typedef float value_t;
typedef int64_t index_t;
typedef uint16_t coefficient_t;

#ifdef INDICATE_PROGRESS
static const std::chrono::milliseconds time_step(40);
#endif

static const std::string clear_line("\r\033[K");

static const size_t num_coefficient_bits = 8;
static const index_t coefficient_mask = (static_cast<index_t>(1) << num_coefficient_bits) - 1;

static const index_t max_simplex_index =
    (1l << (8 * sizeof(index_t) - 1 - num_coefficient_bits)) - 1;

void check_overflow(index_t i) {
	if
#ifdef USE_COEFFICIENTS
	    (i > max_simplex_index)
#else
	    (i < 0)
#endif
		throw std::overflow_error("simplex index " + std::to_string((uint64_t)i) +
		                          " in filtration is larger than maximum index " +
		                          std::to_string(max_simplex_index));
}

class binomial_coeff_table {
	std::vector<std::vector<index_t>> B;

public:
	binomial_coeff_table(index_t n, index_t k) : B(n + 1) {
		for (index_t i = 0; i <= n; ++i) {
			B[i].resize(k + 1, 0);
			B[i][0] = 1;
			for (index_t j = 1; j < std::min(i, k + 1); ++j)
				B[i][j] = B[i - 1][j - 1] + B[i - 1][j];
			if (i <= k) B[i][i] = 1;
			check_overflow(B[i][std::min(i >> 1, k)]);
		}
	}

	index_t operator()(index_t n, index_t k) const {
		assert(n < B.size() && k < B[n].size() && n >= k - 1);
		return B[n][k];
	}
};

bool is_prime(const coefficient_t n) {
	if (!(n & 1) || n < 2) return n == 2;
	for (coefficient_t p = 3; p <= n / p; p += 2)
		if (!(n % p)) return false;
	return true;
}

std::vector<coefficient_t> multiplicative_inverse_vector(const coefficient_t m) {
	std::vector<coefficient_t> inverse(m);
	inverse[1] = 1;
	// m = a * (m / a) + m % a
	// Multipying with inverse(a) * inverse(m % a):
	// 0 = inverse(m % a) * (m / a) + inverse(a)  (mod m)
	for (coefficient_t a = 2; a < m; ++a) inverse[a] = m - (inverse[m % a] * (m / a)) % m;
	return inverse;
}

#ifdef USE_COEFFICIENTS

typedef index_t entry_t;

entry_t make_entry(index_t i, coefficient_t c) { return (i << num_coefficient_bits) | c; }
index_t get_index(const entry_t& e) { return (e >> num_coefficient_bits); }
index_t get_coefficient(const entry_t& e) { return (e & coefficient_mask); }
void set_coefficient(entry_t& e, const coefficient_t c) { e = (e & ~coefficient_mask) | c; }

//std::ostream& operator<<(std::ostream& stream, const entry_t& e) {
//    stream << get_index(e) << ":" << get_coefficient(e);
//    return stream;
//}

#else

typedef index_t entry_t;
const index_t get_index(const entry_t& i) { return i; }
index_t get_coefficient(const entry_t& i) { return 1; }
entry_t make_entry(index_t _index, coefficient_t _value) { return entry_t(_index); }
void set_coefficient(entry_t& e, const coefficient_t c) {}

#endif

const entry_t& get_entry(const entry_t& e) { return e; }

typedef std::pair<value_t, index_t> diameter_index_t;
value_t get_diameter(const diameter_index_t& i) { return i.first; }
index_t get_index(const diameter_index_t& i) { return i.second; }

typedef std::pair<index_t, value_t> index_diameter_t;
index_t get_index(const index_diameter_t& i) { return i.first; }
value_t get_diameter(const index_diameter_t& i) { return i.second; }

struct diameter_entry_t : std::pair<value_t, entry_t> {
	using std::pair<value_t, entry_t>::pair;
	diameter_entry_t(value_t _diameter, index_t _index, coefficient_t _coefficient)
	    : diameter_entry_t(_diameter, make_entry(_index, _coefficient)) {}
	diameter_entry_t(const diameter_index_t& _diameter_index, coefficient_t _coefficient)
	    : diameter_entry_t(get_diameter(_diameter_index),
	                       make_entry(get_index(_diameter_index), _coefficient)) {}
	diameter_entry_t(const index_t& _index) : diameter_entry_t(0, _index, 0) {}
};

const entry_t& get_entry(const diameter_entry_t& p) { return p.second; }
entry_t& get_entry(diameter_entry_t& p) { return p.second; }
const index_t get_index(const diameter_entry_t& p) { return get_index(get_entry(p)); }
const coefficient_t get_coefficient(const diameter_entry_t& p) {
	return get_coefficient(get_entry(p));
}
const value_t& get_diameter(const diameter_entry_t& p) { return p.first; }
void set_coefficient(diameter_entry_t& p, const coefficient_t c) {
	set_coefficient(get_entry(p), c);
}

template <typename Entry> struct greater_diameter_or_smaller_index {
	bool operator()(const Entry& a, const Entry& b) const {
		return (get_diameter(a) > get_diameter(b)) ||
		       ((get_diameter(a) == get_diameter(b)) && (get_index(a) < get_index(b)));
	}
};

enum compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

template <compressed_matrix_layout Layout> struct compressed_distance_matrix {
	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	compressed_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(std::move(_distances)), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
		assert(distances.size() == size() * (size() - 1) / 2);
		init_rows();
	}

	template <typename DistanceMatrix>
	compressed_distance_matrix(const DistanceMatrix& mat)
	    : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size()) {
		init_rows();

		for (size_t i = 1; i < size(); ++i)
			for (size_t j = 0; j < i; ++j) rows[i][j] = mat(i, j);
	}

	value_t operator()(const index_t i, const index_t j) const;
	size_t size() const { return rows.size(); }
	void init_rows();
};

typedef compressed_distance_matrix<LOWER_TRIANGULAR> compressed_lower_distance_matrix;
typedef compressed_distance_matrix<UPPER_TRIANGULAR> compressed_upper_distance_matrix;

template <> void compressed_lower_distance_matrix::init_rows() {
	value_t* pointer = &distances[0];
	for (size_t i = 1; i < size(); ++i) {
		rows[i] = pointer;
		pointer += i;
	}
}

template <> void compressed_upper_distance_matrix::init_rows() {
	value_t* pointer = &distances[0] - 1;
	for (size_t i = 0; i < size() - 1; ++i) {
		rows[i] = pointer;
		pointer += size() - i - 2;
	}
}

template <>
value_t compressed_lower_distance_matrix::operator()(const index_t i, const index_t j) const {
	return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
}

template <>
value_t compressed_upper_distance_matrix::operator()(const index_t i, const index_t j) const {
	return i == j ? 0 : i > j ? rows[j][i] : rows[i][j];
}

struct sparse_distance_matrix {
	std::vector<std::vector<index_diameter_t>> neighbors;

	index_t num_edges;

	sparse_distance_matrix(std::vector<std::vector<index_diameter_t>>&& _neighbors,
	                       index_t _num_edges)
	    : neighbors(std::move(_neighbors)), num_edges(_num_edges) {}

	template <typename DistanceMatrix>
	sparse_distance_matrix(const DistanceMatrix& mat, const value_t threshold)
	    : neighbors(mat.size()), num_edges(0) {

		for (size_t i = 0; i < size(); ++i)
			for (size_t j = 0; j < size(); ++j)
				if (i != j && mat(i, j) <= threshold) {
					++num_edges;
					neighbors[i].push_back({j, mat(i, j)});
				}
	}

	size_t size() const { return neighbors.size(); }
};

struct euclidean_distance_matrix {
	std::vector<std::vector<value_t>> points;

	euclidean_distance_matrix(std::vector<std::vector<value_t>>&& _points)
	    : points(std::move(_points)) {
		for (auto p : points) { assert(p.size() == points.front().size()); }
	}

	value_t operator()(const index_t i, const index_t j) const {
		assert(i < points.size());
		assert(j < points.size());
		return std::sqrt(std::inner_product(
		    points[i].begin(), points[i].end(), points[j].begin(), value_t(), std::plus<value_t>(),
		    [](value_t u, value_t v) { return (u - v) * (u - v); }));
	}

	size_t size() const { return points.size(); }
};

class union_find {
	std::vector<index_t> parent;
	std::vector<uint8_t> rank;

public:
	union_find(const index_t n) : parent(n), rank(n, 0) {
		for (index_t i = 0; i < n; ++i) parent[i] = i;
	}

	index_t find(index_t x) {
		index_t y = x, z;
		while ((z = parent[y]) != y) y = z;
		while ((z = parent[x]) != y) {
			parent[x] = y;
			x = z;
		}
		return z;
	}

	void link(index_t x, index_t y) {
		if ((x = find(x)) == (y = find(y))) return;
		if (rank[x] > rank[y])
			parent[y] = x;
		else {
			parent[x] = y;
			if (rank[x] == rank[y]) ++rank[y];
		}
	}
};

template <typename T> T begin(std::pair<T, T>& p) { return p.first; }
template <typename T> T end(std::pair<T, T>& p) { return p.second; }

template <typename ValueType> class compressed_sparse_matrix {
    std::vector<std::vector<ValueType>*> columns;

public:
    using Column = std::vector<ValueType>;

    compressed_sparse_matrix(size_t n):
        columns(n, nullptr)                                 {}
    ~compressed_sparse_matrix()                             { for(Column* x : columns) delete x; }

    size_t size() const { return columns.size(); }

    Column* column(const index_t index)
        { return mrzv::atomic_ref<Column*>(columns[index]).load(); }

    Column* exchange(const index_t index, Column* desired)
        { return mrzv::atomic_ref<Column*>(columns[index]).exchange(desired); }

    void store(const index_t index, Column* desired)
        { return mrzv::atomic_ref<Column*>(columns[index]).store(desired); }

    bool update(const index_t index, Column*& expected, Column* desired)
        { return mrzv::atomic_ref<Column*>(columns[index]).compare_exchange_weak(expected, desired); }
};


template <class Predicate>
index_t get_max(index_t top, const index_t bottom, const Predicate pred) {
	if (!pred(top)) {
		index_t count = top - bottom;
		while (count > 0) {
			index_t step = count >> 1, mid = top - step;
			if (!pred(mid)) {
				top = mid - 1;
				count -= step + 1;
			} else
				count = step;
		}
	}
	return top;
}

template <typename DistanceMatrix> class ripser {
	const DistanceMatrix dist;
	const index_t n, dim_max;
	const value_t threshold;
	const float ratio;
	const coefficient_t modulus;
	const unsigned num_threads;
	const binomial_coeff_table binomial_coeff;
	const std::vector<coefficient_t> multiplicative_inverse;

	struct entry_hash {
		std::size_t operator()(const entry_t& e) const {
			return std::hash<index_t>()(::get_index(e));
		}
	};

	struct equal_index {
		bool operator()(const entry_t& e, const entry_t& f) const {
			return ::get_index(e) == ::get_index(f);
		}
	};

	typedef hash_map<entry_t, entry_t, entry_hash, equal_index> entry_hash_map;

	typedef compressed_sparse_matrix<diameter_entry_t> Matrix;
	typedef typename Matrix::Column MatrixColumn;

	typedef std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>,
			                    greater_diameter_or_smaller_index<diameter_entry_t>> WorkingColumn;

public:
	ripser(DistanceMatrix&& _dist, index_t _dim_max, value_t _threshold, float _ratio,
	       coefficient_t _modulus, unsigned _num_threads)
	    : dist(std::move(_dist)), n(dist.size()),
	      dim_max(std::min(_dim_max, index_t(dist.size() - 2))), threshold(_threshold),
	      ratio(_ratio), modulus(_modulus),
#ifndef USING_SERIAL
		  num_threads((_num_threads == 0 ? std::thread::hardware_concurrency() : _num_threads)),
#else
		  num_threads(1),
#endif
		  binomial_coeff(n, dim_max + 2),
	      multiplicative_inverse(multiplicative_inverse_vector(_modulus)) {}

	index_t get_max_vertex(const index_t idx, const index_t k, const index_t n) const {
		return get_max(n, k - 1, [&](index_t w) -> bool { return (binomial_coeff(w, k) <= idx); });
	}

	index_t get_edge_index(const index_t i, const index_t j) const {
		return binomial_coeff(i, 2) + j;
	}

	template <typename OutputIterator>
	OutputIterator get_simplex_vertices(index_t idx, const index_t dim, index_t n,
	                                    OutputIterator out) const {
		--n;
		for (index_t k = dim + 1; k > 0; --k) {
			n = get_max_vertex(idx, k, n);
			*out++ = n;
			idx -= binomial_coeff(n, k);
		}
		return out;
	}

	class simplex_coboundary_enumerator;

	void assemble_columns_to_reduce(std::vector<diameter_index_t>& simplices,
	                                std::vector<diameter_index_t>& columns_to_reduce,
	                                entry_hash_map& pivot_column_index, index_t dim) {

#ifdef INDICATE_PROGRESS
		std::cerr << clear_line << "assembling columns" << std::flush;
		std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + time_step;
#endif

		--dim;
		columns_to_reduce.clear();
		std::vector<diameter_index_t> next_simplices;

		const size_t chunk_size = 1024;
		std::atomic<size_t> achunk {0};
#ifdef INDICATE_PROGRESS
		std::atomic<int> progress {0};
#endif

		std::vector<std::vector<diameter_index_t>> next_simplices_vec(num_threads), columns_to_reduce_vec(num_threads);
#if defined(USE_TBB)
		tbb::parallel_for<unsigned>(0, num_threads,
				[&](unsigned i) {
#else
		std::vector<std::future<void>> handles;
		for (unsigned i = 0; i < num_threads; ++i)
			handles.emplace_back(std::async(std::launch::async, [&,i]() {
#endif
				std::vector<diameter_index_t>& next_simplices    = next_simplices_vec[i];
				std::vector<diameter_index_t>& columns_to_reduce = columns_to_reduce_vec[i];
#ifdef INDICATE_PROGRESS
				int indicate_progress = progress++;
#endif
				size_t cur_chunk = achunk++;
				while(cur_chunk * chunk_size < simplices.size()) {
					size_t from = cur_chunk * chunk_size;
					size_t to   = std::min((cur_chunk + 1) * chunk_size, simplices.size());

					for (size_t idx = from; idx < to; ++idx) {
						auto& simplex = simplices[idx];
						simplex_coboundary_enumerator cofacets(diameter_entry_t(simplex, 1), dim, *this);
						while (cofacets.has_next(false)) {
#ifdef INDICATE_PROGRESS
							if (indicate_progress == 0) {
								if (std::chrono::steady_clock::now() > next) {
									std::cerr << clear_line << "assembling columns (processing "
											  << idx << "/" << simplices.size() << " simplices)" << std::flush;
									next = std::chrono::steady_clock::now() + time_step;
								}
							}
#endif
							auto cofacet = cofacets.next();
							if (get_diameter(cofacet) <= threshold) {

								next_simplices.push_back({get_diameter(cofacet), get_index(cofacet)});

								if (pivot_column_index.find(get_entry(cofacet)) == pivot_column_index.end())
									columns_to_reduce.push_back({get_diameter(cofacet), get_index(cofacet)});
							}
						}
					}
					cur_chunk = achunk++;
				}
#if defined(USE_TBB)
            });
#else
			}));
        handles.clear();
#endif

		// figure out offsets to put everything together and resize
		std::vector<size_t> simplices_prefix { 0 }, columns_to_reduce_prefix { 0 };
		for (unsigned i = 0; i < num_threads; ++i) {
			simplices_prefix.push_back(simplices_prefix.back() + next_simplices_vec[i].size());
			columns_to_reduce_prefix.push_back(columns_to_reduce_prefix.back() + columns_to_reduce_vec[i].size());
		}
		next_simplices.resize(simplices_prefix.back());
		columns_to_reduce.resize(columns_to_reduce_prefix.back());

		// copy into the arrays
#if defined(USE_TBB)
		tbb::parallel_for<unsigned>(0, num_threads,
				[&](unsigned i) {
#else
		for (unsigned i = 0; i < num_threads; ++i)
			handles.emplace_back(std::async(std::launch::async, [&,i]() {
#endif
				size_t k = 0;
				for (size_t j = simplices_prefix[i]; j < simplices_prefix[i+1]; ++j)
					next_simplices[j] = next_simplices_vec[i][k++];
				
				k = 0;
				for (size_t j = columns_to_reduce_prefix[i]; j < columns_to_reduce_prefix[i+1]; ++j)
					columns_to_reduce[j] = columns_to_reduce_vec[i][k++];
#if defined(USE_TBB)
            });
#else
			}));
		handles.clear();		// force execution
#endif

		simplices.swap(next_simplices);

#ifdef INDICATE_PROGRESS
		std::cerr << clear_line << "sorting " << columns_to_reduce.size() << " columns"
		          << std::flush;
#endif

#ifdef USING_SERIAL
		std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
		          greater_diameter_or_smaller_index<diameter_index_t>());
#elif defined(USE_TBB)
		tbb::parallel_sort(columns_to_reduce.begin(), columns_to_reduce.end(),
		          greater_diameter_or_smaller_index<diameter_index_t>());
#else
		boost::sort::parallel::parallel_sort(columns_to_reduce.begin(), columns_to_reduce.end(),
				  greater_diameter_or_smaller_index<diameter_index_t>(), num_threads);
#endif
#ifdef INDICATE_PROGRESS
		std::cerr << clear_line << std::flush;
#endif
	}

	void compute_dim_0_pairs(std::vector<diameter_index_t>& edges,
	                         std::vector<diameter_index_t>& columns_to_reduce) {
#ifdef PRINT_PERSISTENCE_PAIRS
		std::cout << "persistence intervals in dim 0:" << std::endl;
#endif

		union_find dset(n);

		edges = get_edges();
		std::sort(edges.rbegin(), edges.rend(),
		          greater_diameter_or_smaller_index<diameter_index_t>());
		std::vector<index_t> vertices_of_edge(2);
		for (auto e : edges) {
			get_simplex_vertices(get_index(e), 1, n, vertices_of_edge.rbegin());
			index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

			if (u != v) {
#ifdef PRINT_PERSISTENCE_PAIRS
				if (get_diameter(e) != 0)
					std::cout << " [0," << get_diameter(e) << ")" << std::endl;
#endif
				dset.link(u, v);
			} else
				columns_to_reduce.push_back(e);
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

#ifdef PRINT_PERSISTENCE_PAIRS
		for (index_t i = 0; i < n; ++i)
			if (dset.find(i) == i) std::cout << " [0, )" << std::endl;
#endif
	}

	diameter_entry_t pop_pivot(WorkingColumn& column) {
		diameter_entry_t pivot(-1);
#ifdef USE_COEFFICIENTS
		while (!column.empty()) {
			if (get_coefficient(pivot) == 0)
				pivot = column.top();
			else if (get_index(column.top()) != get_index(pivot))
				return pivot;
			else
				set_coefficient(pivot,
				                (get_coefficient(pivot) + get_coefficient(column.top())) % modulus);
			column.pop();
		}
		return (get_coefficient(pivot) == 0) ? -1 : pivot;
#else
		while (!column.empty()) {
			pivot = column.top();
			column.pop();
			if (column.empty() || get_index(column.top()) != get_index(pivot)) return pivot;
			column.pop();
		}
		return -1;
#endif
	}

	diameter_entry_t get_pivot(WorkingColumn& column) {
		diameter_entry_t result = pop_pivot(column);
		if (get_index(result) != -1) column.push(result);
		return result;
	}

	std::pair<diameter_entry_t,bool> init_coboundary_and_get_pivot(const diameter_entry_t simplex,
	                                               WorkingColumn& working_coboundary, const index_t& dim,
	                                               entry_hash_map& pivot_column_index, Matrix& reduction_matrix,
												   const size_t index_column_to_reduce) {
		thread_local static std::vector<diameter_entry_t> cofacet_entries;
		bool check_for_emergent_pair = true;
		cofacet_entries.clear();
		simplex_coboundary_enumerator cofacets(simplex, dim, *this);
		while (cofacets.has_next()) {
			diameter_entry_t cofacet = cofacets.next();
			if (get_diameter(cofacet) <= threshold) {
				cofacet_entries.push_back(cofacet);
				if (check_for_emergent_pair && (get_diameter(simplex) == get_diameter(cofacet))) {
					if (pivot_column_index.find(get_entry(cofacet)) == pivot_column_index.end()) {
						if (pivot_column_index.insert({get_entry(cofacet), make_entry(index_column_to_reduce, get_coefficient(cofacet))}).second)
							return { cofacet, true };
					}
					check_for_emergent_pair = false;
				}
			}
		}
		for (auto cofacet : cofacet_entries) working_coboundary.push(cofacet);
		return { get_pivot(working_coboundary), false };
	}

	void add_simplex_coboundary(const diameter_entry_t simplex, const index_t& dim,
	                            WorkingColumn& working_reduction_column, WorkingColumn& working_coboundary, bool add_diagonal = true) {
		if (add_diagonal) working_reduction_column.push(simplex);
		simplex_coboundary_enumerator cofacets(simplex, dim, *this);
		while (cofacets.has_next()) {
			diameter_entry_t cofacet = cofacets.next();
			if (get_diameter(cofacet) <= threshold) working_coboundary.push(cofacet);
		}
	}

	void add_coboundary(MatrixColumn* reduction_column_to_add,
	                    const std::vector<diameter_index_t>& columns_to_reduce,
	                    const size_t index_column_to_add, const coefficient_t factor, const size_t& dim,
	                    WorkingColumn& working_reduction_column, WorkingColumn& working_coboundary, bool add_diagonal = true) {
		diameter_entry_t column_to_add(columns_to_reduce[index_column_to_add], factor);
		add_simplex_coboundary(column_to_add, dim, working_reduction_column, working_coboundary, add_diagonal);

		if (!reduction_column_to_add) return;
		for (diameter_entry_t simplex : *reduction_column_to_add) {
			set_coefficient(simplex, get_coefficient(simplex) * factor % modulus);
			add_simplex_coboundary(simplex, dim, working_reduction_column, working_coboundary);
		}
	}

	MatrixColumn* generate_column(WorkingColumn&& working_reduction_column) {
		if (working_reduction_column.empty())
			return nullptr;

		MatrixColumn column;
		while (true) {
			diameter_entry_t e = pop_pivot(working_reduction_column);
			if (get_index(e) == -1) break;
			assert(get_coefficient(e) > 0);
			column.push_back(e);
		}

		if (column.empty()) return nullptr;
		return new MatrixColumn(std::move(column));
	}

	template<class F>
	void foreach(const std::vector<diameter_index_t>& columns_to_reduce, const F& f) {
#if defined(INDICATE_PROGRESS) && !defined(USE_SERIAL)
		std::atomic<int> progress(0);
		std::cerr << clear_line << "Starting reduction of " << columns_to_reduce.size() << " columns" << std::endl;
#endif
#ifdef USE_SERIAL
		int epoch_counter = 0;
		mrzv::MemoryManager<MatrixColumn> memory_manager(epoch_counter, 1);

		for (size_t index_column_to_reduce = 0; index_column_to_reduce < columns_to_reduce.size();
		     ++index_column_to_reduce) {
			size_t next = f(index_column_to_reduce, true, memory_manager);
			assert(next == index_column_to_reduce);
		}
#elif defined(USE_SHUFFLED_SERIAL) // meant for debugging
		std::vector<size_t> indices(mrzv::counting_iterator(0), mrzv::counting_iterator(columns_to_reduce.size()));
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(indices.begin(), indices.end(), g);

		int epoch_counter = 0;
		mrzv::MemoryManager<MatrixColumn> memory_manager(epoch_counter, 1);

		size_t count = 0;
		for (size_t index_column_to_reduce : indices) {
			bool first = true;
			size_t next;
			do {
				next = index_column_to_reduce;
				index_column_to_reduce = f(next, first, memory_manager);
				first = false;
			} while (next != index_column_to_reduce);

			if (++count == 1024) {
				memory_manager.quiescent();
				count = 0;
			}
		}
#elif defined(USE_PARALLEL_STL)
		int epoch_counter = 0;
		std::for_each(std::execution::par, mrzv::counting_iterator(0), mrzv::counting_iterator(columns_to_reduce.size()),
				[&](size_t index_column_to_reduce) {
			thread_local int count = 0;
			thread_local mrzv::MemoryManager<MatrixColumn> memory_manager(epoch_counter, std::thread::hardware_concurrency());

			bool first = true;
			size_t next;
			do {
				next = index_column_to_reduce;
				index_column_to_reduce = f(next, first, memory_manager);
				first = false;
			} while (next != index_column_to_reduce);

			if (++count == 1024) {
				memory_manager.quiescent();
				count = 0;
			}
		});
#elif defined(USE_TBB)
		int epoch_counter = 0;
		tbb::parallel_for<size_t>(0, columns_to_reduce.size(),
				[&](size_t index_column_to_reduce) {
			thread_local int count = 0;
			thread_local mrzv::MemoryManager<MatrixColumn> memory_manager(epoch_counter, num_threads);

			bool first = true;
			size_t next;
			do {
				next = index_column_to_reduce;
				index_column_to_reduce = f(next, first, memory_manager);
				first = false;
			} while (next != index_column_to_reduce);

			if (++count == 1024) {
				memory_manager.quiescent();
				count = 0;
			}
		});
#else	// default: hand-rolled chunking
		const size_t chunk_size = 1024;
		size_t chunk = 0;
		unsigned n_threads = num_threads;

		int epoch_counter = 0;
		std::vector<std::thread> threads;
		for (unsigned t = 0; t < n_threads; ++t)
			threads.emplace_back([&]() {
				mrzv::atomic_ref<size_t> achunk(chunk);

				mrzv::MemoryManager<MatrixColumn> memory_manager(epoch_counter, n_threads);

#ifdef INDICATE_PROGRESS
				int indicate_progress = progress++;
				std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + time_step;
#endif
				
				size_t cur_chunk = achunk++;
				while(cur_chunk * chunk_size < columns_to_reduce.size()) {
					size_t from = cur_chunk * chunk_size;
					size_t to   = std::min((cur_chunk + 1) * chunk_size, columns_to_reduce.size());
#ifdef INDICATE_PROGRESS
					if (indicate_progress == 0) {
						if (std::chrono::steady_clock::now() > next) {
							std::cerr << clear_line << "reducing columns " << from << " - " << to
									  << "/" << columns_to_reduce.size()
									  << std::flush;
							next = std::chrono::steady_clock::now() + time_step;
						}
					}
#endif
					for (size_t idx = from; idx < to; ++idx) {
						size_t index_column_to_reduce = idx;
						bool first = true;
						size_t next;
						do {
							next = index_column_to_reduce;
							index_column_to_reduce = f(next, first, memory_manager);
							first = false;
						} while (next != index_column_to_reduce);
					}
					cur_chunk = achunk++;
					memory_manager.quiescent();
				}
			});

		for (auto& thread : threads)
			thread.join();
#endif
	}

	// debug only
	diameter_entry_t get_column_pivot(MatrixColumn* column,
	                    const std::vector<diameter_index_t>& columns_to_reduce,
	                    const size_t index, const coefficient_t factor, const size_t& dim) {
		WorkingColumn tmp_working_reduction_column, tmp_working_coboundary;
		add_coboundary(column, columns_to_reduce, index,
					   1, dim, tmp_working_reduction_column, tmp_working_coboundary);
		return get_pivot(tmp_working_coboundary);
	}

	void compute_pairs(const std::vector<diameter_index_t>& columns_to_reduce,
	                   entry_hash_map& pivot_column_index, const index_t dim) {

#if defined(PRINT_PERSISTENCE_PAIRS)
		std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
#endif

		Matrix reduction_matrix(columns_to_reduce.size());

#if defined(PRINT_PERSISTENCE_PAIRS) && !defined(USING_SERIAL)
		// extra vector is a work-around inability to store floats in the hash_map
		typedef hash_map<entry_t, size_t, entry_hash, equal_index> entry_diameter_index_map;
		std::atomic<size_t> last_diameter_index { 0 };
		std::vector<value_t> diameters(columns_to_reduce.size());
		entry_diameter_index_map deaths;
		deaths.reserve(columns_to_reduce.size());
#endif

#if defined(INDICATE_PROGRESS) && defined(USING_SERIAL)
		std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + time_step;
#endif
		foreach(columns_to_reduce, [&](index_t index_column_to_reduce, bool first, mrzv::MemoryManager<MatrixColumn>& memory_manager) {
			diameter_entry_t column_to_reduce(columns_to_reduce[index_column_to_reduce], 1);
			value_t diameter = get_diameter(column_to_reduce);

			WorkingColumn working_reduction_column, working_coboundary;

			diameter_entry_t pivot;
			if (first) {
				bool emergent;
				std::tie(pivot,emergent) = init_coboundary_and_get_pivot(
					column_to_reduce, working_coboundary, dim, pivot_column_index, reduction_matrix, index_column_to_reduce);
				if (emergent)
					return index_column_to_reduce;
			} else {
				MatrixColumn* reduction_column_to_reduce = reduction_matrix.column(index_column_to_reduce);
				add_coboundary(reduction_column_to_reduce, columns_to_reduce, index_column_to_reduce,
							   1, dim, working_reduction_column, working_coboundary, false);
				pivot = get_pivot(working_coboundary);
			}

			while (true) {
#if defined(INDICATE_PROGRESS) && defined(USING_SERIAL)
				if (std::chrono::steady_clock::now() > next) {
					std::cerr << clear_line << "reducing column " << index_column_to_reduce + 1
					          << "/" << columns_to_reduce.size() << " (diameter " << diameter << ")"
					          << std::flush;
					next = std::chrono::steady_clock::now() + time_step;
				}
#endif
				if (get_index(pivot) != -1) {
					auto pair = pivot_column_index.find(get_entry(pivot));
					if (pair != pivot_column_index.end()) {
						entry_t old_entry_column_to_add;
						index_t index_column_to_add;
						MatrixColumn* reduction_column_to_add;
						entry_t entry_column_to_add = pivot_column_index.value(pair);
						do
						{
							old_entry_column_to_add = entry_column_to_add;

							index_column_to_add = get_index(entry_column_to_add);

							reduction_column_to_add = reduction_matrix.column(index_column_to_add);

							// this is a weaker check than in the original lockfree
							// persistence paper (it would suffice that the pivot
							// in reduction_column_to_add) hasn't changed, but
							// given that matrix V is stored, rather than matrix R,
							// it's easier to check that pivot_column_index entry
							// we read hasn't changed
							// TODO: think through memory orders, and whether we need to adjust anything
							entry_column_to_add = pivot_column_index.value(pair);
						} while (old_entry_column_to_add != entry_column_to_add);

						if (index_column_to_add < index_column_to_reduce)
						{
							// pivot to the left; usual reduction
							coefficient_t factor =
								modulus - get_coefficient(pivot) *
											  multiplicative_inverse[get_coefficient(entry_column_to_add)] %
											  modulus;

							add_coboundary(reduction_column_to_add, columns_to_reduce, index_column_to_add,
										   factor, dim, working_reduction_column, working_coboundary);

							pivot = get_pivot(working_coboundary);
						} else {
							// pivot to the right
							MatrixColumn* new_column = generate_column(std::move(working_reduction_column));
							MatrixColumn* previous = reduction_matrix.exchange(index_column_to_reduce, new_column);
							assert(get_index(get_column_pivot(new_column, columns_to_reduce, index_column_to_reduce, 1, dim)) == get_index(pivot));
							memory_manager.retire(previous);

							if (pivot_column_index.update(pair, entry_column_to_add, make_entry(index_column_to_reduce, get_coefficient(pivot)))) {
								return index_column_to_add;
							} else {
								continue; // re-read the pair
							}
						}
					} else {
#if defined(PRINT_PERSISTENCE_PAIRS) && defined(USING_SERIAL)
						value_t death = get_diameter(pivot);
						if (death > diameter * ratio) {
#ifdef INDICATE_PROGRESS
							std::cerr << clear_line << std::flush;
#endif
							std::cout << " [" << diameter << "," << death << ")" << (first ? "" : " <- correction") << std::endl;
						}
#endif
#if defined(PRINT_PERSISTENCE_PAIRS) && !defined(USING_SERIAL)
						size_t location = last_diameter_index++;
						diameters[location] = get_diameter(pivot);
						deaths.insert({get_entry(pivot), location});
#endif
						MatrixColumn* new_column = generate_column(std::move(working_reduction_column));
						MatrixColumn* previous = reduction_matrix.exchange(index_column_to_reduce, new_column);
						memory_manager.retire(previous);

						assert(get_index(get_column_pivot(new_column, columns_to_reduce, index_column_to_reduce, 1, dim)) == get_index(pivot));

						// equivalent to CAS in the original algorithm
						auto insertion_result = pivot_column_index.insert({get_entry(pivot), make_entry(index_column_to_reduce, get_coefficient(pivot))});
						if (!insertion_result.second)		// failed to insert, somebody got there before us, continue reduction
							continue;						// TODO: insertion_result.first is the new pair; could elide and extra atomic load

						break;
					}
				} else {
					// TODO: these will need special attention, if output happens after the reduction, not during
#if defined(PRINT_PERSISTENCE_PAIRS) && defined(USING_SERIAL)
#ifdef INDICATE_PROGRESS
					std::cerr << clear_line << std::flush;
#endif
					std::cout << " [" << diameter << ", )" << (first ? "" : " <- correction") << std::endl;
#endif
					break;
				}
			}
			return index_column_to_reduce;
		});
#if defined(INDICATE_PROGRESS)
		std::cerr << clear_line << std::flush;
#endif
#if defined(PRINT_PERSISTENCE_PAIRS) && !defined(USING_SERIAL)
		pivot_column_index.foreach([&](const typename entry_hash_map::value_type& x) {
			auto it = deaths.find(x.first);
			if (it == deaths.end()) return;
			value_t death = diameters[it->second];
			value_t birth = get_diameter(columns_to_reduce[get_index(x.second)]);
			if (death > birth * ratio)
				std::cout << " [" << birth << "," << death << ")" << std::endl;
		});
		// TODO: this doesn't print unpaired values
#endif
	}

	std::vector<diameter_index_t> get_edges();

	void compute_barcodes() {
		std::vector<diameter_index_t> simplices, columns_to_reduce;

		compute_dim_0_pairs(simplices, columns_to_reduce);

		for (index_t dim = 1; dim <= dim_max; ++dim) {
			entry_hash_map pivot_column_index;
			pivot_column_index.reserve(columns_to_reduce.size());

			compute_pairs(columns_to_reduce, pivot_column_index, dim);

			if (dim < dim_max)
				assemble_columns_to_reduce(simplices, columns_to_reduce, pivot_column_index,
				                           dim + 1);
		}
	}
};

template <> class ripser<compressed_lower_distance_matrix>::simplex_coboundary_enumerator {
	index_t idx_below, idx_above, v, k;
	std::vector<index_t> vertices;
	const diameter_entry_t simplex;
	const coefficient_t modulus;
	const compressed_lower_distance_matrix& dist;
	const binomial_coeff_table& binomial_coeff;

public:
	simplex_coboundary_enumerator(const diameter_entry_t _simplex, const index_t _dim,
	                              const ripser& parent)
	    : idx_below(get_index(_simplex)), idx_above(0), v(parent.n - 1), k(_dim + 1),
	      vertices(_dim + 1), simplex(_simplex), modulus(parent.modulus), dist(parent.dist),
	      binomial_coeff(parent.binomial_coeff) {
		parent.get_simplex_vertices(get_index(_simplex), _dim, parent.n, vertices.rbegin());
	}

	bool has_next(bool all_cofacets = true) {
		return (v >= k && (all_cofacets || binomial_coeff(v, k) > idx_below));
	}

	diameter_entry_t next() {
		while ((binomial_coeff(v, k) <= idx_below)) {
			idx_below -= binomial_coeff(v, k);
			idx_above += binomial_coeff(v, k + 1);
			--v;
			--k;
			assert(k != -1);
		}
		value_t cofacet_diameter = get_diameter(simplex);
		for (index_t w : vertices) cofacet_diameter = std::max(cofacet_diameter, dist(v, w));
		index_t cofacet_index = idx_above + binomial_coeff(v--, k + 1) + idx_below;
		coefficient_t cofacet_coefficient =
		    (k & 1 ? modulus - 1 : 1) * get_coefficient(simplex) % modulus;
		return diameter_entry_t(cofacet_diameter, cofacet_index, cofacet_coefficient);
	}
};

template <> class ripser<sparse_distance_matrix>::simplex_coboundary_enumerator {
	const ripser& parent;
	index_t idx_below, idx_above, k;
	std::vector<index_t> vertices;
	const diameter_entry_t simplex;
	const coefficient_t modulus;
	const sparse_distance_matrix& dist;
	const binomial_coeff_table& binomial_coeff;
	static thread_local std::vector<std::vector<index_diameter_t>::const_reverse_iterator> neighbor_it;
	static thread_local std::vector<std::vector<index_diameter_t>::const_reverse_iterator> neighbor_end;
	index_diameter_t neighbor;

public:
	simplex_coboundary_enumerator(const diameter_entry_t _simplex, const index_t _dim,
	                              const ripser& _parent)
	    : parent(_parent), idx_below(get_index(_simplex)), idx_above(0), k(_dim + 1),
	      vertices(_dim + 1), simplex(_simplex), modulus(parent.modulus), dist(parent.dist),
	      binomial_coeff(parent.binomial_coeff) {
		neighbor_it.clear();
		neighbor_end.clear();

		parent.get_simplex_vertices(idx_below, _dim, parent.n, vertices.rbegin());
		for (auto v : vertices) {
			neighbor_it.push_back(dist.neighbors[v].rbegin());
			neighbor_end.push_back(dist.neighbors[v].rend());
		}
	}

	bool has_next(bool all_cofacets = true) {
		for (auto &it0 = neighbor_it[0], &end0 = neighbor_end[0]; it0 != end0; ++it0) {
			neighbor = *it0;
			for (size_t idx = 1; idx < neighbor_it.size(); ++idx) {
				auto &it = neighbor_it[idx], end = neighbor_end[idx];
				while (get_index(*it) > get_index(neighbor))
					if (++it == end) return false;
				if (get_index(*it) != get_index(neighbor))
					goto continue_outer;
				else
					neighbor = std::max(neighbor, *it);
			}
			while (k > 0 && vertices[k - 1] > get_index(neighbor)) {
				if (!all_cofacets) return false;
				idx_below -= binomial_coeff(vertices[k - 1], k);
				idx_above += binomial_coeff(vertices[k - 1], k + 1);
				--k;
			}
			return true;
		continue_outer:;
		}
		return false;
	}

	diameter_entry_t next() {
		++neighbor_it[0];
		value_t cofacet_diameter = std::max(get_diameter(simplex), get_diameter(neighbor));
		index_t cofacet_index = idx_above + binomial_coeff(get_index(neighbor), k + 1) + idx_below;
		coefficient_t cofacet_coefficient =
		    (k & 1 ? modulus - 1 : 1) * get_coefficient(simplex) % modulus;
		return diameter_entry_t(cofacet_diameter, cofacet_index, cofacet_coefficient);
	}
};

thread_local std::vector<std::vector<index_diameter_t>::const_reverse_iterator> ripser<sparse_distance_matrix>::simplex_coboundary_enumerator::neighbor_it;
thread_local std::vector<std::vector<index_diameter_t>::const_reverse_iterator> ripser<sparse_distance_matrix>::simplex_coboundary_enumerator::neighbor_end;

template <> std::vector<diameter_index_t> ripser<compressed_lower_distance_matrix>::get_edges() {
	std::vector<diameter_index_t> edges;
	std::vector<index_t> vertices(2);
	for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
		get_simplex_vertices(index, 1, dist.size(), vertices.rbegin());
		value_t length = dist(vertices[0], vertices[1]);
		if (length <= threshold) edges.push_back({length, index});
	}
	return edges;
}

template <> std::vector<diameter_index_t> ripser<sparse_distance_matrix>::get_edges() {
	std::vector<diameter_index_t> edges;
	for (index_t i = 0; i < n; ++i)
		for (auto n : dist.neighbors[i]) {
			index_t j = get_index(n);
			if (i > j) edges.push_back({get_diameter(n), get_edge_index(i, j)});
		}
	return edges;
}

enum file_format {
	LOWER_DISTANCE_MATRIX,
	UPPER_DISTANCE_MATRIX,
	DISTANCE_MATRIX,
	POINT_CLOUD,
	DIPHA,
	SPARSE,
	BINARY
};

static const uint16_t endian_check(0xff00);
static const bool is_big_endian = *reinterpret_cast<const uint8_t*>(&endian_check);

template <typename T> T read(std::istream& input_stream) {
	T result;
	char* p = reinterpret_cast<char*>(&result);
	if (input_stream.read(p, sizeof(T)).gcount() != sizeof(T)) return T();
	if (is_big_endian) std::reverse(p, p + sizeof(T));
	return result;
}

compressed_lower_distance_matrix read_point_cloud(std::istream& input_stream) {
	std::vector<std::vector<value_t>> points;

	std::string line;
	value_t value;
	while (std::getline(input_stream, line)) {
		std::vector<value_t> point;
		std::istringstream s(line);
		while (s >> value) {
			point.push_back(value);
			s.ignore();
		}
		if (!point.empty()) points.push_back(point);
		assert(point.size() == points.front().size());
	}

	euclidean_distance_matrix eucl_dist(std::move(points));
	index_t n = eucl_dist.size();
	std::cout << "point cloud with " << n << " points in dimension "
	          << eucl_dist.points.front().size() << std::endl;

	std::vector<value_t> distances;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < i; ++j) distances.push_back(eucl_dist(i, j));

	return compressed_lower_distance_matrix(std::move(distances));
}

sparse_distance_matrix read_sparse_distance_matrix(std::istream& input_stream) {
	std::vector<std::vector<index_diameter_t>> neighbors;
	index_t num_edges = 0;

	std::string line;
	while (std::getline(input_stream, line)) {
		std::istringstream s(line);
		size_t i, j;
		value_t value;
		s >> i;
		s >> j;
		s >> value;
		if (i != j) {
			neighbors.resize(std::max({neighbors.size(), i + 1, j + 1}));
			neighbors[i].push_back({j, value});
			neighbors[j].push_back({i, value});
			++num_edges;
		}
	}

	for (size_t i = 0; i < neighbors.size(); ++i)
		std::sort(neighbors[i].begin(), neighbors[i].end());

	return sparse_distance_matrix(std::move(neighbors), num_edges);
}

compressed_lower_distance_matrix read_lower_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_upper_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;
	value_t value;
	while (input_stream >> value) {
		distances.push_back(value);
		input_stream.ignore();
	}

	return compressed_lower_distance_matrix(compressed_upper_distance_matrix(std::move(distances)));
}

compressed_lower_distance_matrix read_distance_matrix(std::istream& input_stream) {
	std::vector<value_t> distances;

	std::string line;
	value_t value;
	for (int i = 0; std::getline(input_stream, line); ++i) {
		std::istringstream s(line);
		for (int j = 0; j < i && s >> value; ++j) {
			distances.push_back(value);
			s.ignore();
		}
	}

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_dipha(std::istream& input_stream) {
	if (read<int64_t>(input_stream) != 8067171840) {
		std::cerr << "input is not a Dipha file (magic number: 8067171840)" << std::endl;
		exit(-1);
	}

	if (read<int64_t>(input_stream) != 7) {
		std::cerr << "input is not a Dipha distance matrix (file type: 7)" << std::endl;
		exit(-1);
	}

	index_t n = read<int64_t>(input_stream);

	std::vector<value_t> distances;

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (i > j)
				distances.push_back(read<double>(input_stream));
			else
				read<double>(input_stream);

	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_binary(std::istream& input_stream) {
	std::vector<value_t> distances;
	while (!input_stream.eof()) distances.push_back(read<value_t>(input_stream));
	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_file(std::istream& input_stream, const file_format format) {
	switch (format) {
	case LOWER_DISTANCE_MATRIX:
		return read_lower_distance_matrix(input_stream);
	case UPPER_DISTANCE_MATRIX:
		return read_upper_distance_matrix(input_stream);
	case DISTANCE_MATRIX:
		return read_distance_matrix(input_stream);
	case POINT_CLOUD:
		return read_point_cloud(input_stream);
	case DIPHA:
		return read_dipha(input_stream);
	default:
		return read_binary(input_stream);
	}
}

void print_usage_and_exit(int exit_code) {
	std::cerr
	    << "Usage: "
	    << "ripser "
	    << "[options] [filename]" << std::endl
	    << std::endl
	    << "Options:" << std::endl
	    << std::endl
	    << "  --help           print this screen" << std::endl
	    << "  --format         use the specified file format for the input. Options are:"
	    << std::endl
	    << "                     lower-distance (lower triangular distance matrix; default)"
	    << std::endl
	    << "                     upper-distance (upper triangular distance matrix)" << std::endl
	    << "                     distance       (full distance matrix)" << std::endl
	    << "                     point-cloud    (point cloud in Euclidean space)" << std::endl
	    << "                     dipha          (distance matrix in DIPHA file format)" << std::endl
	    << "                     sparse         (sparse distance matrix in sparse triplet format)"
	    << std::endl
	    << "                     binary         (lower triangular distance matrix in binary format)"
	    << std::endl
	    << "  --dim <k>        compute persistent homology up to dimension k" << std::endl
	    << "  --threshold <t>  compute Rips complexes up to diameter t" << std::endl
#ifdef USE_COEFFICIENTS
	    << "  --modulus <p>    compute homology with coefficients in the prime field Z/pZ"
	    << std::endl
#endif
#ifndef USING_SERIAL
	    << "  --threads <t>    number of threads to use"
	    << std::endl
#endif
	    << "  --ratio <r>      only show persistence pairs with death/birth ratio > r" << std::endl
	    << std::endl;
	exit(exit_code);
}

int main(int argc, char** argv) {
	const char* filename = nullptr;

	file_format format = DISTANCE_MATRIX;

	index_t dim_max = 1;
	value_t threshold = std::numeric_limits<value_t>::max();
	float ratio = 1;
	coefficient_t modulus = 2;
	unsigned num_threads = 0;

	for (index_t i = 1; i < argc; ++i) {
		const std::string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--dim") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			dim_max = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--threshold") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			threshold = std::stof(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--ratio") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			ratio = std::stof(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--format") {
			std::string parameter = std::string(argv[++i]);
			if (parameter.rfind("lower", 0) == 0)
				format = LOWER_DISTANCE_MATRIX;
			else if (parameter.rfind("upper", 0) == 0)
				format = UPPER_DISTANCE_MATRIX;
			else if (parameter.rfind("dist", 0) == 0)
				format = DISTANCE_MATRIX;
			else if (parameter.rfind("point", 0) == 0)
				format = POINT_CLOUD;
			else if (parameter == "dipha")
				format = DIPHA;
			else if (parameter == "sparse")
				format = SPARSE;
			else if (parameter == "binary")
				format = BINARY;
			else
				print_usage_and_exit(-1);
#ifdef USE_COEFFICIENTS
		} else if (arg == "--modulus") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			modulus = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size() || !is_prime(modulus)) print_usage_and_exit(-1);
#endif
#ifndef USING_SERIAL
		} else if (arg == "--threads") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			num_threads = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size() || !is_prime(modulus)) print_usage_and_exit(-1);
#endif
		} else {
			if (filename) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}

#ifdef USE_TBB
	tbb::task_scheduler_init init(num_threads);
#endif

	std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "couldn't open file " << filename << std::endl;
		exit(-1);
	}

	if (format == SPARSE) {
		sparse_distance_matrix dist =
		    read_sparse_distance_matrix(filename ? file_stream : std::cin);
		std::cout << "sparse distance matrix with " << dist.size() << " points and "
		          << dist.num_edges << "/" << (dist.size() * (dist.size() - 1)) / 2 << " entries"
		          << std::endl;

		ripser<sparse_distance_matrix>(std::move(dist), dim_max, threshold, ratio, modulus, num_threads)
		    .compute_barcodes();
	} else {
		compressed_lower_distance_matrix dist =
		    read_file(filename ? file_stream : std::cin, format);

		value_t min = std::numeric_limits<value_t>::infinity(),
		        max = -std::numeric_limits<value_t>::infinity(), max_finite = max;
		int num_edges = 0;

		if (threshold == std::numeric_limits<value_t>::max()) {
			value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
			for (size_t i = 0; i < dist.size(); ++i) {
				value_t r_i = -std::numeric_limits<value_t>::infinity();
				for (size_t j = 0; j < dist.size(); ++j) r_i = std::max(r_i, dist(i, j));
				enclosing_radius = std::min(enclosing_radius, r_i);
			}
			threshold = enclosing_radius;
		}

		for (auto d : dist.distances) {
			min = std::min(min, d);
			max = std::max(max, d);
			max_finite =
			    d != std::numeric_limits<value_t>::infinity() ? std::max(max, d) : max_finite;
			if (d <= threshold) ++num_edges;
		}
		std::cout << "value range: [" << min << "," << max_finite << "]" << std::endl;

		if (threshold >= max) {
			std::cout << "distance matrix with " << dist.size() << " points" << std::endl;
			ripser<compressed_lower_distance_matrix>(std::move(dist), dim_max, threshold, ratio,
			                                         modulus, num_threads)
			    .compute_barcodes();
		} else {
			std::cout << "sparse distance matrix with " << dist.size() << " points and "
			          << num_edges << "/" << (dist.size() * dist.size() - 1) / 2 << " entries"
			          << std::endl;

			ripser<sparse_distance_matrix>(sparse_distance_matrix(std::move(dist), threshold),
			                               dim_max, threshold, ratio, modulus, num_threads)
			    .compute_barcodes();
		}
		exit(0);
	}
}
