/*

Ripser: a lean C++ code for computation of Vietoris-Rips persistence barcodes

Copyright 2015-2016 Ulrich Bauer.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

//#define ASSEMBLE_REDUCTION_MATRIX
//#define USE_COEFFICIENTS

//#define INDICATE_PROGRESS
//#define PRINT_PERSISTENCE_PAIRS

// #define PYTHON_EXTENSION

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>


template <class Key, class T> class hash_map : public std::unordered_map<Key, T> {};
typedef float value_t;

typedef int64_t index_t;
typedef int16_t coefficient_t;

class binomial_coeff_table {
	std::vector<std::vector<index_t>> B;
	index_t n_max, k_max;

public:
	binomial_coeff_table(index_t n, index_t k) {
		n_max = n;
		k_max = k;

		B.resize(n + 1);
		for (index_t i = 0; i <= n; i++) {
			B[i].resize(k + 1);
			for (index_t j = 0; j <= std::min(i, k); j++) {
				if (j == 0 || j == i)
					B[i][j] = 1;
				else
					B[i][j] = B[i - 1][j - 1] + B[i - 1][j];
			}
		}
	}

	index_t operator()(index_t n, index_t k) const {
		assert(n <= n_max);
		assert(k <= k_max);
		return B[n][k];
	}
};

bool is_prime(const coefficient_t n) {
	if (!(n & 1) || n < 2) return n == 2;
	for (coefficient_t p = 3, q = n / p, r = n % p; p <= q; p += 2, q = n / p, r = n % p)
		if (!r) return false;
	return true;
}

coefficient_t normalize(const coefficient_t n, const coefficient_t modulus) {
    return n > modulus/2 ? n - modulus : n;
}

std::vector<coefficient_t> multiplicative_inverse_vector(const coefficient_t m) {
	std::vector<coefficient_t> inverse(m);
	inverse[1] = 1;
	// m = a * (m / a) + m % a
	// Multipying with inverse(a) * inverse(m % a):
	// 0 = inverse(m % a) * (m / a)  + inverse(a)  (mod m)
	for (coefficient_t a = 2; a < m; ++a) inverse[a] = m - (inverse[m % a] * (m / a)) % m;
	return inverse;
}

index_t get_next_vertex(index_t& v, const index_t idx, const index_t k, const binomial_coeff_table& binomial_coeff) {
	if (binomial_coeff(v, k) > idx) {
		index_t count = v;
		while (count > 0) {
			index_t i = v;
			index_t step = count >> 1;
			i -= step;
			if (binomial_coeff(i, k) > idx) {
				v = --i;
				count -= step + 1;
			} else
				count = step;
		}
	}
	assert(binomial_coeff(v, k) <= idx);
	assert(binomial_coeff(v + 1, k) > idx);
	return v;
}

template <typename OutputIterator>
OutputIterator get_simplex_vertices(index_t idx, const index_t dim, index_t v,
                                    const binomial_coeff_table& binomial_coeff, OutputIterator out) {
	--v;
	for (index_t k = dim + 1; k > 0; --k) {
		get_next_vertex(v, idx, k, binomial_coeff);
		*out++ = v;
		idx -= binomial_coeff(v, k);
	}
	return out;
}

std::vector<index_t> vertices_of_simplex(const index_t simplex_index, const index_t dim, const index_t n,
                                         const binomial_coeff_table& binomial_coeff) {
	std::vector<index_t> vertices;
	get_simplex_vertices(simplex_index, dim, n, binomial_coeff, std::back_inserter(vertices));
	return vertices;
}

#ifdef USE_COEFFICIENTS
struct entry_t {
	index_t index : 8 * (sizeof(index_t) - sizeof(coefficient_t));
	coefficient_t coefficient;
	entry_t(index_t _index, coefficient_t _coefficient) : index(_index), coefficient(_coefficient) {}
	entry_t(index_t _index) : index(_index), coefficient(1) {}
	entry_t() : index(0), coefficient(1) {}
} __attribute__((packed));

//static_assert(sizeof(entry_t) == sizeof(index_t), "size of entry_t is not the same as index_t");

entry_t make_entry(index_t _index, coefficient_t _coefficient) { return entry_t(_index, _coefficient); }
index_t get_index(entry_t e) { return e.index; }
index_t get_coefficient(entry_t e) { return e.coefficient; }
void set_coefficient(entry_t& e, const coefficient_t c) { e.coefficient = c; }

bool operator==(const entry_t& e1, const entry_t& e2) {
	return get_index(e1) == get_index(e2) && get_coefficient(e1) == get_coefficient(e2);
}

std::ostream& operator<<(std::ostream& stream, const entry_t& e) {
	stream << get_index(e) << ":" << get_coefficient(e);
	return stream;
}

#else

typedef index_t entry_t;
const index_t get_index(entry_t i) { return i; }
index_t get_coefficient(entry_t i) { return 1; }
entry_t make_entry(index_t _index, coefficient_t _value) { return entry_t(_index); }
void set_coefficient(index_t& e, const coefficient_t c) {}

#endif

const entry_t& get_entry(const entry_t& e) { return e; }

template <typename Entry> struct smaller_index {
	bool operator()(const Entry& a, const Entry& b) { return get_index(a) < get_index(b); }
};

class diameter_index_t : public std::pair<value_t, index_t> {
public:
	diameter_index_t() : std::pair<value_t, index_t>() {}
	diameter_index_t(std::pair<value_t, index_t> p) : std::pair<value_t, index_t>(p) {}
};
value_t get_diameter(diameter_index_t i) { return i.first; }
index_t get_index(diameter_index_t i) { return i.second; }

class diameter_entry_t : public std::pair<value_t, entry_t> {
public:
	diameter_entry_t(std::pair<value_t, entry_t> p) : std::pair<value_t, entry_t>(p) {}
	diameter_entry_t(entry_t e) : std::pair<value_t, entry_t>(0, e) {}
	diameter_entry_t() : diameter_entry_t(0) {}
	diameter_entry_t(value_t _diameter, index_t _index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(_diameter, make_entry(_index, _coefficient)) {}
	diameter_entry_t(diameter_index_t _diameter_index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(get_diameter(_diameter_index),
	                                  make_entry(get_index(_diameter_index), _coefficient)) {}
	diameter_entry_t(diameter_index_t _diameter_index) : diameter_entry_t(_diameter_index, 1) {}
};

const entry_t& get_entry(const diameter_entry_t& p) { return p.second; }
entry_t& get_entry(diameter_entry_t& p) { return p.second; }
const index_t get_index(const diameter_entry_t& p) { return get_index(get_entry(p)); }
const coefficient_t get_coefficient(const diameter_entry_t& p) { return get_coefficient(get_entry(p)); }
const value_t& get_diameter(const diameter_entry_t& p) { return p.first; }
void set_coefficient(diameter_entry_t& p, const coefficient_t c) { set_coefficient(get_entry(p), c); }

template <typename Entry> struct greater_diameter_or_smaller_index {
	bool operator()(const Entry& a, const Entry& b) {
		return (get_diameter(a) > get_diameter(b)) ||
		       ((get_diameter(a) == get_diameter(b)) && (get_index(a) < get_index(b)));
	}
};

template <typename DistanceMatrix> class rips_filtration_comparator {
public:
	const DistanceMatrix& dist;
	const index_t dim;

private:
	mutable std::vector<index_t> vertices;
	const binomial_coeff_table& binomial_coeff;

public:
	rips_filtration_comparator(const DistanceMatrix& _dist, const index_t _dim,
	                           const binomial_coeff_table& _binomial_coeff)
	    : dist(_dist), dim(_dim), vertices(_dim + 1), binomial_coeff(_binomial_coeff){};

	value_t diameter(const index_t index) const {
		value_t diam = 0;
		get_simplex_vertices(index, dim, dist.size(), binomial_coeff, vertices.begin());

		for (index_t i = 0; i <= dim; ++i)
			for (index_t j = 0; j < i; ++j) { diam = std::max(diam, dist(vertices[i], vertices[j])); }
		return diam;
	}

	bool operator()(const index_t a, const index_t b) const {
		assert(a < binomial_coeff(dist.size(), dim + 1));
		assert(b < binomial_coeff(dist.size(), dim + 1));

		return greater_diameter_or_smaller_index<diameter_index_t>()(diameter_index_t(diameter(a), a),
		                                                             diameter_index_t(diameter(b), b));
	}

	template <typename Entry> bool operator()(const Entry& a, const Entry& b) const {
		return operator()(get_index(a), get_index(b));
	}
};

template <class DistanceMatrix> class simplex_coboundary_enumerator {
private:
	index_t idx_below, idx_above, v, k;
	std::vector<index_t> vertices;
	const diameter_entry_t simplex;
	const coefficient_t modulus;
	const DistanceMatrix& dist;
	const binomial_coeff_table& binomial_coeff;

public:
	simplex_coboundary_enumerator(const diameter_entry_t _simplex, index_t _dim, index_t _n,
	                              const coefficient_t _modulus, const DistanceMatrix& _dist,
	                              const binomial_coeff_table& _binomial_coeff)
	    : simplex(_simplex), idx_below(get_index(_simplex)), idx_above(0), v(_n - 1), k(_dim + 1), modulus(_modulus),
	      binomial_coeff(_binomial_coeff), dist(_dist), vertices(_dim + 1) {
		get_simplex_vertices(get_index(_simplex), _dim, _n, binomial_coeff, vertices.begin());
	}

	bool has_next() {
		while ((v != -1) && (binomial_coeff(v, k) <= idx_below)) {
			idx_below -= binomial_coeff(v, k);
			idx_above += binomial_coeff(v, k + 1);

			--v;
			--k;
			assert(k != -1);
		}
		return v != -1;
	}

	index_t next_index() { return idx_above + binomial_coeff(v--, k + 1) + idx_below; }

	diameter_entry_t next() {
		value_t coface_diameter = get_diameter(simplex);
		for (index_t w : vertices) coface_diameter = std::max(coface_diameter, dist(v, w));
		coefficient_t coface_coefficient = (k & 1 ? -1 + modulus : 1) * get_coefficient(simplex) % modulus;
		return diameter_entry_t(coface_diameter, idx_above + binomial_coeff(v--, k + 1) + idx_below,
		                        coface_coefficient);
	}
};

enum compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

template <compressed_matrix_layout Layout> class compressed_distance_matrix {
public:
	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	void init_rows();

	compressed_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(_distances), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
		assert(distances.size() == size() * (size() - 1) / 2);
		init_rows();
	}

	template <typename DistanceMatrix>
	compressed_distance_matrix(const DistanceMatrix& mat)
	    : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size()) {
		init_rows();

		for (index_t i = 1; i < size(); ++i)
			for (index_t j = 0; j < i; ++j) rows[i][j] = mat(i, j);
	}

	value_t operator()(const index_t i, const index_t j) const;

	size_t size() const { return rows.size(); }
};

template <> void compressed_distance_matrix<LOWER_TRIANGULAR>::init_rows() {
	value_t* pointer = &distances[0];
	for (index_t i = 1; i < size(); ++i) {
		rows[i] = pointer;
		pointer += i;
	}
}

template <> void compressed_distance_matrix<UPPER_TRIANGULAR>::init_rows() {
	value_t* pointer = &distances[0] - 1;
	for (index_t i = 0; i < size() - 1; ++i) {
		rows[i] = pointer;
		pointer += size() - i - 2;
	}
}

template <> value_t compressed_distance_matrix<UPPER_TRIANGULAR>::operator()(index_t i, index_t j) const {
	if (i > j) std::swap(i, j);
	return i == j ? 0 : rows[i][j];
}

template <> value_t compressed_distance_matrix<LOWER_TRIANGULAR>::operator()(index_t i, index_t j) const {
	if (i > j) std::swap(i, j);
	return i == j ? 0 : rows[j][i];
}

typedef compressed_distance_matrix<LOWER_TRIANGULAR> compressed_lower_distance_matrix;
typedef compressed_distance_matrix<UPPER_TRIANGULAR> compressed_upper_distance_matrix;

class euclidean_distance_matrix {
public:
	std::vector<std::vector<value_t>> points;

	euclidean_distance_matrix(std::vector<std::vector<value_t>>&& _points) : points(_points) {}

	value_t operator()(const index_t i, const index_t j) const {
		return std::sqrt(std::inner_product(points[i].begin(), points[i].end(), points[j].begin(), value_t(),
		                                    std::plus<value_t>(),
		                                    [](value_t u, value_t v) { return (u - v) * (u - v); }));
	}

	size_t size() const { return points.size(); }
};

class union_find {
	std::vector<index_t> parent;
	std::vector<uint8_t> rank;

public:
	union_find(index_t n) : parent(n), rank(n, 0) {
		for (index_t i = 0; i < n; ++i) parent[i] = i;
	}

	index_t find(index_t x) {
		index_t y = x, z = parent[y];
		while (z != y) {
			y = z;
			z = parent[y];
		}
		y = parent[x];
		while (z != y) {
			parent[x] = z;
			x = y;
			y = parent[x];
		}
		return z;
	}
	void link(index_t x, index_t y) {
		x = find(x);
		y = find(y);
		if (x == y) return;
		if (rank[x] > rank[y])
			parent[y] = x;
		else {
			parent[x] = y;
			if (rank[x] == rank[y]) ++rank[y];
		}
	}
};

template <typename Heap> diameter_entry_t pop_pivot(Heap& column, coefficient_t modulus) {
	if (column.empty())
		return diameter_entry_t(-1);
	else {
		auto pivot = column.top();

#ifdef USE_COEFFICIENTS
		coefficient_t coefficient = 0;
		do {
			coefficient = (coefficient + get_coefficient(column.top())) % modulus;
			column.pop();

			if (coefficient == 0) {
				if (column.empty())
					return diameter_entry_t(-1);
				else
					pivot = column.top();
			}
		} while (!column.empty() && get_index(column.top()) == get_index(pivot));
		if (get_index(pivot) != -1) { set_coefficient(pivot, coefficient); }
#else
		column.pop();
		while (!column.empty() && get_index(column.top()) == get_index(pivot)) {
			column.pop();
			if (column.empty())
				return diameter_entry_t(-1);
			else {
				pivot = column.top();
				column.pop();
			}
		}
#endif
		return pivot;
	}
}

template <typename Heap> diameter_entry_t get_pivot(Heap& column, coefficient_t modulus) {
	diameter_entry_t result = pop_pivot(column, modulus);
	if (get_index(result) != -1) column.push(result);
	return result;
}

template <typename ValueType> class compressed_sparse_matrix {
	std::vector<size_t> bounds;
	std::vector<ValueType> entries;

public:
	size_t size() const { return bounds.size(); }

	typename std::vector<ValueType>::const_iterator cbegin(size_t index) const {
		assert(index < size());
		return index == 0 ? entries.cbegin() : entries.cbegin() + bounds[index - 1];
	}

	typename std::vector<ValueType>::const_iterator cend(size_t index) const {
		assert(index < size());
		return entries.cbegin() + bounds[index];
	}

	template <typename Iterator> void append_column(Iterator begin, Iterator end) {
		for (Iterator it = begin; it != end; ++it) { entries.push_back(*it); }
		bounds.push_back(entries.size());
	}

	void append_column() { bounds.push_back(entries.size()); }

	void push_back(ValueType e) {
		assert(0 < size());
		entries.push_back(e);
		++bounds.back();
	}

	void pop_back() {
		assert(0 < size());
		entries.pop_back();
		--bounds.back();
	}

	template <typename Collection> void append_column(const Collection collection) {
		append_column(collection.cbegin(), collection.cend());
	}
};

template <typename Heap> void push_entry(Heap& column, index_t i, coefficient_t c, value_t diameter) {
	entry_t e = make_entry(i, c);
	column.push(std::make_pair(diameter, e));
}

template <typename Comparator>
void assemble_columns_to_reduce(std::vector<diameter_index_t>& columns_to_reduce,
                                hash_map<index_t, index_t>& pivot_column_index, 
								const Comparator& comp, 
								index_t dim,
                                index_t n, 
								value_t threshold, 
								const binomial_coeff_table& binomial_coeff) {
	
	index_t num_simplices = binomial_coeff(n, dim + 2);

	columns_to_reduce.clear();

	for (index_t index = 0; index < num_simplices; ++index) {
		if (pivot_column_index.find(index) == pivot_column_index.end()) {
			value_t diameter = comp.diameter(index);
			if (diameter <= threshold) columns_to_reduce.push_back(std::make_pair(diameter, index));
		}
	}

	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
	          greater_diameter_or_smaller_index<diameter_index_t>());
}

template <typename DistanceMatrix, typename ComparatorCofaces, typename Comparator>
void compute_pairs(std::vector<diameter_index_t>& columns_to_reduce, 
				   hash_map<index_t, index_t>& pivot_column_index,
				   std::vector<value_t>& birthsanddeaths,
				   int do_cocycles, std::vector<value_t>& allcocycles,
				   index_t dim, index_t n, value_t threshold, coefficient_t modulus,
                   const std::vector<coefficient_t>& multiplicative_inverse, 
				   const DistanceMatrix& dist,
                   const ComparatorCofaces& comp, 
				   const Comparator& comp_prev,
                   const binomial_coeff_table& binomial_coeff) {

#ifdef ASSEMBLE_REDUCTION_MATRIX
	compressed_sparse_matrix<diameter_entry_t> reduction_coefficients;
#else
#ifdef USE_COEFFICIENTS
	std::vector<diameter_entry_t> reduction_coefficients;
#endif
#endif

	std::vector<diameter_entry_t> coface_entries;

	for (index_t i = 0; i < columns_to_reduce.size(); ++i) {
		auto column_to_reduce = columns_to_reduce[i];

#ifdef ASSEMBLE_REDUCTION_MATRIX
		std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>, smaller_index<diameter_entry_t>>
		    reduction_column;
#endif

		std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>,
		                    greater_diameter_or_smaller_index<diameter_entry_t>>
		    working_coboundary;

		value_t diameter = get_diameter(column_to_reduce);

		index_t j = i;

		// start with a dummy pivot entry with coefficient -1 in order to initialize
		// working_coboundary with the coboundary of the simplex with index column_to_reduce
		diameter_entry_t pivot(0, -1, -1 + modulus);

#ifdef ASSEMBLE_REDUCTION_MATRIX
		// initialize reduction_coefficients as identity matrix
		reduction_coefficients.append_column();
		reduction_coefficients.push_back(diameter_entry_t(column_to_reduce, 1));
#else
#ifdef USE_COEFFICIENTS
		reduction_coefficients.push_back(diameter_entry_t(column_to_reduce, 1));
#endif
#endif

		bool might_be_apparent_pair = (i == j);

		do {
			const coefficient_t factor = modulus - get_coefficient(pivot);

#ifdef ASSEMBLE_REDUCTION_MATRIX
			auto coeffs_begin = reduction_coefficients.cbegin(j), coeffs_end = reduction_coefficients.cend(j);
#else
#ifdef USE_COEFFICIENTS
			auto coeffs_begin = &reduction_coefficients[j], coeffs_end = &reduction_coefficients[j] + 1;
#else
			auto coeffs_begin = &columns_to_reduce[j], coeffs_end = &columns_to_reduce[j] + 1;
#endif
#endif

			for (auto it = coeffs_begin; it != coeffs_end; ++it) {
				diameter_entry_t simplex = *it;
				set_coefficient(simplex, get_coefficient(simplex) * factor % modulus);

#ifdef ASSEMBLE_REDUCTION_MATRIX
				reduction_column.push(simplex);
#endif

				coface_entries.clear();
				simplex_coboundary_enumerator<decltype(dist)> cofaces(simplex, dim, n, modulus, dist, binomial_coeff);
				while (cofaces.has_next()) {
					diameter_entry_t coface = cofaces.next();
					if (get_diameter(coface) <= threshold) {
						coface_entries.push_back(coface);
						if (might_be_apparent_pair && (get_diameter(simplex) == get_diameter(coface))) {
							if (pivot_column_index.find(get_index(coface)) == pivot_column_index.end()) {
								pivot = coface;
								goto found_persistence_pair;
							}
							might_be_apparent_pair = false;
						}
					}
				}
				for (auto e : coface_entries) {
					working_coboundary.push(e);
				}
			}

			pivot = get_pivot(working_coboundary, modulus);

			if (get_index(pivot) != -1) {
				auto pair = pivot_column_index.find(get_index(pivot));

				if (pair != pivot_column_index.end()) {
					j = pair->second;
					continue;
				}
			} else {
				//Infinite class; make the death time the diameter of the point cloud
				birthsanddeaths.push_back(diameter);
				birthsanddeaths.push_back(INFINITY);

#ifdef ASSEMBLE_REDUCTION_MATRIX
				if (do_cocycles) {
					//Representative cocycle
					auto cocycle = reduction_column;
					diameter_entry_t e;
					std::vector<index_t> simplex;
					std::vector<value_t> thiscocycle;
					index_t lencocycle = 0;
					while (get_index(e = get_pivot(cocycle, modulus)) != -1) {
						simplex = vertices_of_simplex(get_index(e), dim, n, binomial_coeff);
						for (size_t k = 0; k < simplex.size(); k++) {
							thiscocycle.push_back((value_t)simplex[k]);
						}
						thiscocycle.push_back(normalize(get_coefficient(e), modulus));
						cocycle.pop();
						lencocycle++;
					}
					allcocycles.push_back(lencocycle);
					allcocycles.insert(allcocycles.end(), thiscocycle.begin(), thiscocycle.end());
				}
#endif
				break;
			}

			found_persistence_pair:
			value_t death = get_diameter(pivot);

			if (diameter != death) {
				birthsanddeaths.push_back(diameter);
				birthsanddeaths.push_back(death);
#ifdef ASSEMBLE_REDUCTION_MATRIX
				if (do_cocycles) {
					//Representative cocycle
					auto cocycle = reduction_column;
					diameter_entry_t e;
					std::vector<index_t> simplex;
					std::vector<value_t> thiscocycle;
					index_t lencocycle = 0;
					while (get_index(e = get_pivot(cocycle, modulus)) != -1) {
						simplex = vertices_of_simplex(get_index(e), dim, n, binomial_coeff);
						for (size_t k = 0; k < simplex.size(); k++) {
							thiscocycle.push_back((value_t)simplex[k]);
						}
						thiscocycle.push_back(normalize(get_coefficient(e), modulus));
						cocycle.pop();
						lencocycle++;
					}
					allcocycles.push_back(lencocycle);
					allcocycles.insert(allcocycles.end(), thiscocycle.begin(), thiscocycle.end());
				}
#endif
			}

			pivot_column_index.insert(std::make_pair(get_index(pivot), i));

#ifdef USE_COEFFICIENTS
			const coefficient_t inverse = multiplicative_inverse[get_coefficient(pivot)];
#endif

#ifdef ASSEMBLE_REDUCTION_MATRIX
			// replace current column of reduction_coefficients (with a single diagonal 1 entry)
			// by reduction_column (possibly with a different entry on the diagonal)
			reduction_coefficients.pop_back();
			while (true) {
				diameter_entry_t e = pop_pivot(reduction_column, modulus);
				if (get_index(e) == -1) break;
#ifdef USE_COEFFICIENTS
				set_coefficient(e, inverse * get_coefficient(e) % modulus);
				assert(get_coefficient(e) > 0);
#endif
				reduction_coefficients.push_back(e);
			}
#else
#ifdef USE_COEFFICIENTS
			reduction_coefficients.pop_back();
			reduction_coefficients.push_back(diameter_entry_t(column_to_reduce, inverse));
#endif
#endif
			break;
		} while (true);
	}

}

enum file_format { LOWER_DISTANCE_MATRIX, UPPER_DISTANCE_MATRIX, DISTANCE_MATRIX, POINT_CLOUD, DIPHA };

template <typename T> T read(std::istream& s) {
	T result;
	s.read(reinterpret_cast<char*>(&result), sizeof(T));
	return result; // on little endian: boost::endian::little_to_native(result);
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

	std::cout << "point cloud with " << n << " points in dimension " << eucl_dist.points.front().size() << std::endl;

	std::vector<value_t> distances;

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < i; ++j) distances.push_back(eucl_dist(i, j));

	return compressed_lower_distance_matrix(std::move(distances));
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

compressed_lower_distance_matrix read_file(std::istream& input_stream, file_format format) {
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
	}
}

void print_usage_and_exit(int exit_code) {
	std::cerr << "Usage: "
	          << "ripser "
	          << "[options] [filename]" << std::endl
	          << std::endl
	          << "Options:" << std::endl
	          << std::endl
	          << "  --help           print this screen" << std::endl
	          << "  --format         use the specified file format for the input. Options are:" << std::endl
	          << "                     lower-distance (lower triangular distance matrix; default)" << std::endl
	          << "                     upper-distance (upper triangular distance matrix)" << std::endl
	          << "                     distance       (full distance matrix)" << std::endl
	          << "                     point-cloud    (point cloud in Euclidean space)" << std::endl
	          << "                     dipha          (distance matrix in DIPHA file format)" << std::endl
	          << "  --dim <k>        compute persistent homology up to dimension <k>" << std::endl
	          << "  --threshold <t>  compute Rips complexes up to diameter <t>" << std::endl
	          << "  --output <f>     output persistence pairs to file <f>" << std::endl
#ifdef USE_COEFFICIENTS
	          << "  --modulus <p>    compute homology with coefficients in the prime field Z/<p>Z"
#endif
	          << std::endl;

	exit(exit_code);
}

std::vector<value_t> pythondm(float* D, int N, int modulus, int dim_max, float threshold, int do_cocycles) {
	//Setup distance matrix, coefficient tables, etc
	std::vector<value_t> distances(D, D+N);
	compressed_lower_distance_matrix dist = compressed_lower_distance_matrix(compressed_upper_distance_matrix(std::move(distances)));
	index_t n = dist.size();
	dim_max = std::min((index_t)(dim_max), n-2);
	if (threshold == -1) {
		threshold = INFINITY;
	}

	binomial_coeff_table binomial_coeff(n, dim_max + 2);
	std::vector<coefficient_t> multiplicative_inverse(multiplicative_inverse_vector(modulus));
	std::vector<value_t> birthsanddeaths;
	std::vector<value_t> retvec;
	/*Copy everything into a return array with the following format
	 [NumperClass0, PD0, 
	  NumPerClass1, PD1, Cocycle1_0Len, Cocycle1_0, Cocycle1_1Len, Cocycle1_1, ..., Cocycle1_NumPerClass1Len, Cocycle1_NumPerClass1
	  NumPerClass2, PD2, Cocycle1_0Len, Cocycle1_0, Cocycle1_1Len, Cocycle1_1, ..., Cocycle2_NumPerClass1Len, Cocycle2_NumPerClass1
	  ...
	  ]
	
	where PD0 is b0, d0, b1, d1, ...
	and each cocycle is simplexidx0, simplexidx1, ..., simplexidxd, val mod p
	*/


	//Do 0D persistence first
	std::vector<diameter_index_t> columns_to_reduce;

	{
		union_find dset(n);
		std::vector<diameter_index_t> edges;
		rips_filtration_comparator<decltype(dist)> comp(dist, 1, binomial_coeff);
		for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
			value_t diameter = comp.diameter(index);
			if (diameter <= threshold) edges.push_back(std::make_pair(diameter, index));
		}
		std::sort(edges.rbegin(), edges.rend(), greater_diameter_or_smaller_index<diameter_index_t>());

		std::vector<index_t> vertices_of_edge(2);
		for (auto e : edges) {
			vertices_of_edge.clear();
			get_simplex_vertices(get_index(e), 1, n, binomial_coeff, std::back_inserter(vertices_of_edge));
			index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

			if (u != v) {
				if (get_diameter(e) > 0) {
					birthsanddeaths.push_back(0);
					birthsanddeaths.push_back(get_diameter(e));
				}
				dset.link(u, v);
			} else
				columns_to_reduce.push_back(e);
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

		for (index_t i = 0; i < n; ++i)
			if (dset.find(i) == i) {
				//Have one element representing the infinite class
				birthsanddeaths.push_back(0);
				birthsanddeaths.push_back(INFINITY);
			}
		retvec.push_back(birthsanddeaths.size()/2);
		retvec.insert(retvec.end(), birthsanddeaths.begin(), birthsanddeaths.end());
	}

	for (index_t dim = 1; dim <= dim_max; ++dim) {
		rips_filtration_comparator<decltype(dist)> comp(dist, dim + 1, binomial_coeff);
		rips_filtration_comparator<decltype(dist)> comp_prev(dist, dim, binomial_coeff);

		hash_map<index_t, index_t> pivot_column_index;
		pivot_column_index.reserve(columns_to_reduce.size());

		std::vector<value_t> allcocycles;
		birthsanddeaths.clear();
		
		compute_pairs(columns_to_reduce, pivot_column_index, birthsanddeaths, 
						do_cocycles, allcocycles, dim, n, threshold, modulus, 
						multiplicative_inverse, dist, comp, comp_prev, binomial_coeff);
		
		//Copy over persistence diagram
		retvec.push_back(birthsanddeaths.size()/2);
		retvec.insert(retvec.end(), birthsanddeaths.begin(), birthsanddeaths.end());
		//Copy over cocycles
		retvec.insert(retvec.end(), allcocycles.begin(), allcocycles.end());

		if (dim < dim_max) {
			assemble_columns_to_reduce(columns_to_reduce, pivot_column_index, comp, dim, n, threshold, binomial_coeff);
		}
	}
	return retvec;
}

