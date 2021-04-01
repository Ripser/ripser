/*

 Infiltrator: a lean C++ code for computation of simplicial persistence barcodes

 MIT License

 Copyright (c) 2015–2021 Ulrich Bauer

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
#define PRINT_PERSISTENCE_PAIRS

//#define USE_GOOGLE_HASHMAP

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
#include <ctime>

#ifdef USE_GOOGLE_HASHMAP
#include <sparsehash/dense_hash_map>
template <class Key, class T, class H, class E>
class hash_map : public google::dense_hash_map<Key, T, H, E> {
public:
	explicit hash_map() : google::dense_hash_map<Key, T, H, E>() { this->set_empty_key(-1); }
	inline void reserve(size_t hint) { this->resize(hint); }
};
#else
template <class Key, class T, class H, class E>
class hash_map : public std::unordered_map<Key, T, H, E> {};
#endif

typedef float value_t;
typedef __int128_t index_t;
typedef uint16_t coefficient_t;

std::ostream& operator<<(std::ostream& stream, const __int128_t& i) {
	stream << (size_t)(i);
	return stream;
}

std::istream& operator>>(std::istream& stream, __int128_t& i) {
	size_t s;
	stream >> s;
	i = (__int128_t)(s);
	return stream;
}

#ifdef INDICATE_PROGRESS
static const std::chrono::milliseconds time_step(40);
#endif

static const std::string clear_line("\r\033[K");

#ifdef USE_COEFFICIENTS
static const size_t num_coefficient_bits = 8;

static const index_t max_simplex_index =
    (1l << (8 * sizeof(index_t) - 1 - num_coefficient_bits)) - 1;
#endif

void check_overflow(index_t i) {
	if
#ifdef USE_COEFFICIENTS
	    (i > max_simplex_index)
#else
	    (i < 0)
#endif
		throw std::overflow_error("simplex index " + std::to_string((uint64_t)i) +
		                          " in filtration is larger than maximum index"// + std::to_string(max_simplex_index)
								  );
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

struct __attribute__((packed)) entry_t {
	index_t index : 8 * sizeof(index_t) - num_coefficient_bits;
	coefficient_t coefficient : num_coefficient_bits;
	entry_t(index_t _index, coefficient_t _coefficient)
	    : index(_index), coefficient(_coefficient) {}
	entry_t(index_t _index) : index(_index), coefficient(0) {}
	entry_t() : index(0), coefficient(0) {}
};

static_assert(sizeof(entry_t) == sizeof(index_t), "size of entry_t is not the same as index_t");

entry_t make_entry(index_t i, coefficient_t c) { return entry_t(i, c); }
index_t get_index(const entry_t& e) { return e.index; }
index_t get_coefficient(const entry_t& e) { return e.coefficient; }
void set_coefficient(entry_t& e, const coefficient_t c) { e.coefficient = c; }

std::ostream& operator<<(std::ostream& stream, const entry_t& e) {
	stream << get_index(e) << ":" << get_coefficient(e);
	return stream;
}

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

template <typename Entry> struct smaller_diameter_or_greater_index {
	bool operator()(const Entry& a, const Entry& b) {
		return (get_diameter(a) < get_diameter(b)) ||
		((get_diameter(a) == get_diameter(b)) && (get_index(a) > get_index(b)));
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

	mutable std::vector<std::vector<index_diameter_t>::const_reverse_iterator> neighbor_it;
	mutable std::vector<std::vector<index_diameter_t>::const_reverse_iterator> neighbor_end;

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

index_t compute_index(const std::vector<index_t> vertices, const binomial_coeff_table& B) {

	index_t index = 0, j = vertices.size() - 1;
	for (index_t i : vertices) {
		index += B(i, j + 1);
		j--;
	}
	return index;
}

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
	std::vector<size_t> bounds;
	std::vector<ValueType> entries;

	typedef typename std::vector<ValueType>::iterator iterator;
	typedef std::pair<iterator, iterator> iterator_pair;

public:
	size_t size() const { return bounds.size(); }

	iterator_pair subrange(const index_t index) {
		return {entries.begin() + (index == 0 ? 0 : bounds[index - 1]),
		        entries.begin() + bounds[index]};
	}

	void append_column() { bounds.push_back(entries.size()); }

	void push_back(const ValueType e) {
		assert(0 < size());
		entries.push_back(e);
		++bounds.back();
	}
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

class infiltrator {
	std::vector<std::unordered_map<index_t, value_t>> filtration;
    index_t n, dim_max;
	value_t threshold;
	float ratio;
	coefficient_t modulus;
	const binomial_coeff_table binomial_coeff;
	const std::vector<coefficient_t> multiplicative_inverse;
	mutable std::vector<diameter_entry_t> cofacet_entries;

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

	typedef hash_map<entry_t, size_t, entry_hash, equal_index> entry_hash_map;

public:
    infiltrator(std::vector<std::unordered_map<index_t, value_t>>&& _filtration, index_t _n, index_t _dim_max, value_t _threshold,
	       float _ratio, coefficient_t _modulus)
	    : filtration(std::move(_filtration)), n(_n), dim_max(_dim_max), threshold(_threshold),
	      ratio(_ratio), modulus(_modulus), binomial_coeff(n, dim_max + 2),
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

	value_t compute_diameter(const index_t index, index_t dim) const {
		auto pair = filtration[dim].find(index);
		if (pair != filtration[dim].end()) return pair->second;
		else return std::numeric_limits<value_t>::infinity();
	}

	class simplex_boundary_enumerator {
	private:
		index_t idx_below, idx_above, v, k, face_dim;
		std::vector<index_t> vertices;
		const diameter_entry_t simplex;
		const coefficient_t modulus;
		const binomial_coeff_table& binomial_coeff;
		const infiltrator& parent;
		
	public:
		simplex_boundary_enumerator(const diameter_entry_t _simplex, index_t _dim,
									const infiltrator& _parent)
		: idx_below(get_index(_simplex)), idx_above(0), v(_parent.n - 1), k(_dim + 1),
		face_dim(_dim - 1), vertices(_dim + 1),  simplex(_simplex), modulus(_parent.modulus),
		binomial_coeff(_parent.binomial_coeff), parent(_parent) {
			parent.get_simplex_vertices(get_index(_simplex), _dim, parent.n, vertices.begin());
		}
		
		bool has_next() {
			v = (v == -1) ? -1 : parent.get_max_vertex(idx_below, k, v);
			return (v != -1) && (binomial_coeff(v, k) <= idx_below);
		}
		
		diameter_entry_t next() {
			index_t face_index = idx_above - binomial_coeff(v, k) + idx_below;
			
			value_t face_diameter = parent.compute_diameter(face_index, face_dim);
			
			coefficient_t face_coefficient = (k & 1 ? -1 + modulus : 1) * get_coefficient(simplex) % modulus;
			
			idx_below -= binomial_coeff(v, k);
			idx_above += binomial_coeff(v, k - 1);
			
			--v;
			--k;
			
			return diameter_entry_t(face_diameter, face_index, face_coefficient);
			;
		}
	};
	
	std::vector<index_t> vertices_of_simplex(const index_t simplex_index, const index_t dim) const {
		std::vector<index_t> vertices(dim + 1);
		get_simplex_vertices(simplex_index, dim, n, vertices.rbegin());
		return vertices;
	}
	


	void assemble_columns_to_reduce(std::vector<diameter_index_t>& columns_to_reduce,
                                    entry_hash_map& pivot_column_index, index_t dim);

	void compute_dim_0_pairs(std::vector<diameter_index_t>& columns_to_reduce) {
#ifdef PRINT_PERSISTENCE_PAIRS
		std::cout << "persistence intervals in dim 0:" << std::endl;
#endif

		union_find dset(n);

		std::vector<diameter_index_t> edges = get_edges();

		std::sort(edges.begin(), edges.end(),
		          smaller_diameter_or_greater_index<diameter_index_t>());
		std::vector<index_t> vertices_of_edge(2);
		for (auto e : edges) {
			get_simplex_vertices(get_index(e), 1, n, vertices_of_edge.begin());
			index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

			if (u != v) {
#ifdef PRINT_PERSISTENCE_PAIRS
				if (get_diameter(e) > 0)
					std::cout << " ["
					<< std::max(filtration[0].find(vertices_of_edge[0])->second,
								filtration[0].find(vertices_of_edge[1])->second)
					<< "," << get_diameter(e) << ")" << std::endl;
#endif
				dset.link(u, v);
			} else
				columns_to_reduce.push_back(e);
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

#ifdef PRINT_PERSISTENCE_PAIRS
		for (index_t i = 0; i < n; ++i)
            if (dset.find(i) == i) {
                std::cout << " [" << filtration[0].find(i)->second << ", )" << std::endl;
            }
#endif
    }

	template <typename Column> diameter_entry_t pop_pivot(Column& column) {
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

	template <typename Column> diameter_entry_t get_pivot(Column& column) {
		diameter_entry_t result = pop_pivot(column);
		if (get_index(result) != -1) column.push(result);
		return result;
	}

	template <typename Column>
	diameter_entry_t init_coboundary_and_get_pivot(const diameter_entry_t simplex,
	                                               Column& working_coboundary, const index_t& dim,
	                                               entry_hash_map& pivot_column_index) {
		bool check_for_emergent_pair = true;
		cofacet_entries.clear();
		simplex_boundary_enumerator cofacets(simplex, dim, *this);
		while (cofacets.has_next()) {
			diameter_entry_t cofacet = cofacets.next();
			if (get_diameter(cofacet) <= threshold) {
				cofacet_entries.push_back(cofacet);
				if (check_for_emergent_pair && (get_diameter(simplex) == get_diameter(cofacet))) {
					if (pivot_column_index.find(get_entry(cofacet)) == pivot_column_index.end())
						return cofacet;
					check_for_emergent_pair = false;
				}
			}
		}
		for (auto cofacet : cofacet_entries) working_coboundary.push(cofacet);
		return get_pivot(working_coboundary);
	}

	template <typename Column>
	void add_simplex_coboundary(const diameter_entry_t simplex, const index_t& dim,
	                            Column& working_reduction_column, Column& working_coboundary) {
		working_reduction_column.push(simplex);
		simplex_boundary_enumerator cofacets(simplex, dim, *this);
		while (cofacets.has_next()) {
			diameter_entry_t cofacet = cofacets.next();
			if (get_diameter(cofacet) <= threshold) working_coboundary.push(cofacet);
		}
	}

	template <typename Column>
	void add_coboundary(compressed_sparse_matrix<diameter_entry_t>& reduction_matrix,
	                    const std::vector<diameter_index_t>& columns_to_reduce,
	                    const size_t index_column_to_add, const coefficient_t factor, const size_t& dim,
	                    Column& working_reduction_column, Column& working_coboundary) {
		diameter_entry_t column_to_add(columns_to_reduce[index_column_to_add], factor);
		add_simplex_coboundary(column_to_add, dim, working_reduction_column, working_coboundary);

		for (diameter_entry_t simplex : reduction_matrix.subrange(index_column_to_add)) {
			set_coefficient(simplex, get_coefficient(simplex) * factor % modulus);
			add_simplex_coboundary(simplex, dim, working_reduction_column, working_coboundary);
		}
	}

	void compute_pairs(const std::vector<diameter_index_t>& columns_to_reduce,
	                   entry_hash_map& pivot_column_index, const index_t dim) {

#ifdef PRINT_PERSISTENCE_PAIRS
		std::cout << "persistence intervals in dim " << dim - 1 << " (" << columns_to_reduce.size() << " columns):" << std::endl;
#endif

		compressed_sparse_matrix<diameter_entry_t> reduction_matrix;

#ifdef INDICATE_PROGRESS
		std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + time_step;
#endif
		for (size_t index_column_to_reduce = 0; index_column_to_reduce < columns_to_reduce.size();
		     ++index_column_to_reduce) {

			diameter_entry_t column_to_reduce(columns_to_reduce[index_column_to_reduce], 1);
			value_t diameter = get_diameter(column_to_reduce);

			reduction_matrix.append_column();

			std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>,
			                    smaller_diameter_or_greater_index<diameter_entry_t>>
			    working_reduction_column, working_coboundary;

			diameter_entry_t pivot = init_coboundary_and_get_pivot(
			    column_to_reduce, working_coboundary, dim, pivot_column_index);

			while (true) {
#ifdef INDICATE_PROGRESS
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
						entry_t other_pivot = pair->first;
						index_t index_column_to_add = pair->second;
						coefficient_t factor =
						    modulus - get_coefficient(pivot) *
						                  multiplicative_inverse[get_coefficient(other_pivot)] %
						                  modulus;

						add_coboundary(reduction_matrix, columns_to_reduce, index_column_to_add,
						               factor, dim, working_reduction_column, working_coboundary);

						pivot = get_pivot(working_coboundary);
					} else {
#ifdef PRINT_PERSISTENCE_PAIRS
						value_t birth = get_diameter(pivot);
						if (birth * ratio < diameter) {
#ifdef INDICATE_PROGRESS
							std::cerr << clear_line << std::flush;
#endif
							std::cout << " [" << birth << "," << diameter << ")" << std::endl;
						}
#endif
						pivot_column_index.insert({get_entry(pivot), index_column_to_reduce});

						while (true) {
							diameter_entry_t e = pop_pivot(working_reduction_column);
							if (get_index(e) == -1) break;
							assert(get_coefficient(e) > 0);
							reduction_matrix.push_back(e);
						}
						break;
					}
				} else {
#ifdef PRINT_PERSISTENCE_PAIRS
#ifdef INDICATE_PROGRESS
					std::cerr << clear_line << std::flush;
#endif
					std::cout << "+[" << diameter << ", )" << std::endl;
#endif
					break;
				}
			}
		}
#ifdef INDICATE_PROGRESS
		std::cerr << clear_line << std::flush;
#endif
	}

	std::vector<diameter_index_t> get_edges();

	void compute_barcodes() {
		std::vector<diameter_index_t> simplices, columns_to_reduce;

		entry_hash_map pivot_column_index;


		for (index_t dim = dim_max; dim > 0; --dim) {
			assemble_columns_to_reduce(columns_to_reduce, pivot_column_index, dim);
			
			pivot_column_index = entry_hash_map();
			pivot_column_index.reserve(columns_to_reduce.size());

			compute_pairs(columns_to_reduce, pivot_column_index, dim);
		}

		compute_dim_0_pairs(columns_to_reduce);

	}
};

std::vector<diameter_index_t> infiltrator::get_edges() {
	std::vector<diameter_index_t> edges;
	for (auto entry: filtration[1]) {
		index_t index = entry.first;
		value_t diameter = entry.second;
		if (diameter <= threshold) edges.push_back(std::make_pair(diameter, index));
	}
	return edges;
}

void infiltrator::assemble_columns_to_reduce(std::vector<diameter_index_t>& columns_to_reduce,
                                        entry_hash_map& pivot_column_index,
                                        index_t dim) {
	columns_to_reduce.clear();
	
#ifdef INDICATE_PROGRESS
	std::cout << "\033[K"
	<< "assembling " << filtration[dim].size() << " columns" << std::flush << "\n";
#endif
	
	for (auto entry: filtration[dim]) {
		index_t index = entry.first;
		value_t diameter = entry.second;
		if (pivot_column_index.find(index) == pivot_column_index.end()) {
			if (diameter <= threshold)
				columns_to_reduce.push_back(std::make_pair(diameter, index));
#ifdef INDICATE_PROGRESS
			if ((index + 1) % 1000000 == 0)
				std::cout << "\033[K"
				<< "assembled " << columns_to_reduce.size() << " out of "
				<< (index + 1) << "/" << filtration[dim].size() << " columns" << std::flush
				<< "\n";
#endif
		}
	}
	
#ifdef INDICATE_PROGRESS
	std::cout << "\033[K"
	<< "sorting " << columns_to_reduce.size() << " columns" << std::flush << "\n";
#endif
	
	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
			  smaller_diameter_or_greater_index<diameter_index_t>());
#ifdef INDICATE_PROGRESS
	std::cout << "\033[K";
#endif
}

template <typename T> T read(std::istream& s) {
	T result;
	s.read(reinterpret_cast<char*>(&result), sizeof(T));
	return result; // on little endian: boost::endian::little_to_native(result);
}

void print_usage_and_exit(int exit_code) {
	std::cerr
	    << "Usage: "
	    << "infiltrator "
	    << "[options] [filename]" << std::endl
	    << std::endl
	    << "Options:" << std::endl
	    << std::endl
	    << "  --help           print this screen" << std::endl
	    << "  --dim <k>        compute persistent homology up to dimension k" << std::endl
	    << "  --threshold <t>  compute homology up to filtration value <t>" << std::endl
#ifdef USE_COEFFICIENTS
	    << "  --modulus <p>    compute homology with coefficients in the prime field Z/pZ"
	    << std::endl
#endif
	    << "  --ratio <r>      only show persistence pairs with death/birth ratio > r" << std::endl
	    << std::endl;
	exit(exit_code);
}

std::vector<std::unordered_map<index_t, value_t>> read_file(std::istream& input_stream, index_t& n,
                                                            index_t& dim_max) {
	std::string line;
	std::string delimiter = "]";
	std::getline(input_stream, line);
	std::stringstream stream(line);
	index_t num;
	std::vector<index_t> new_simplex;
	std::size_t string_end;

	stream >> n;
	stream >> dim_max;

	const binomial_coeff_table B(n, dim_max + 2);

	std::vector<std::unordered_map<index_t, value_t>> filtration(dim_max + 2);

	while (std::getline(input_stream, line)) {
		string_end = line.find(delimiter);
		std::stringstream stream(line.substr(1, string_end - 1));
		new_simplex.clear();
		while (stream >> num) { new_simplex.push_back(num); }
		std::sort(new_simplex.rbegin(), new_simplex.rend());

		value_t filtration_value;
		(std::stringstream(line.substr(string_end + 1, line.length()))) >> filtration_value;
		index_t index = compute_index(new_simplex, B);
		filtration[new_simplex.size() - 1].insert(std::make_pair(index, filtration_value));
	}
	return filtration;
}

#include <iostream>

int main(int argc, const char* argv[]) {

	const char* filename = nullptr;

	index_t n, dim_max = std::numeric_limits<index_t>::max(), dim;
	value_t threshold = std::numeric_limits<value_t>::max();
	float ratio = 1;
	coefficient_t modulus = 2;

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
#ifdef USE_COEFFICIENTS
		} else if (arg == "--modulus") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			modulus = std::stol(parameter, &next_pos);
			if (next_pos != parameter.size() || !is_prime(modulus)) print_usage_and_exit(-1);
#endif
		} else {
			if (filename) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}

	std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "couldn't open file " << filename << std::endl;
		exit(-1);
	}

	std::vector<std::unordered_map<index_t, value_t>> filtration =
	    read_file(filename ? file_stream : std::cin, n, dim);

	std::cout << "complex of dimension " << dim << " with " << n << " vertices" << std::endl;

	double clock_start = std::clock();

	infiltrator(std::move(filtration), n, std::min(dim, dim_max), threshold, ratio, modulus).compute_barcodes();

	std::cout << "Computed persistent homology in " << (std::clock()-clock_start) / CLOCKS_PER_SEC << " s\n";

}
