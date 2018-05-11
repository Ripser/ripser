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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
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

public:
	binomial_coeff_table(index_t n, index_t k) : B(n + 1) {
		for (index_t i = 0; i <= n; i++) {
			B[i].resize(k + 1);
			for (index_t j = 0; j <= std::min(i, k); j++)
				if (j == 0 || j == i)
					B[i][j] = 1;
				else
					B[i][j] = B[i - 1][j - 1] + B[i - 1][j];
		}
	}

	index_t operator()(index_t n, index_t k) const {
		assert(n < B.size() && k < B[n].size());
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

#ifdef USE_COEFFICIENTS
struct __attribute__((packed)) entry_t {
	index_t index : 8 * (sizeof(index_t) - sizeof(coefficient_t));
	coefficient_t coefficient;
	entry_t(index_t _index, coefficient_t _coefficient)
	    : index(_index), coefficient(_coefficient) {}
	entry_t(index_t _index) : index(_index), coefficient(1) {}
	entry_t() : index(0), coefficient(1) {}
};

//static_assert(sizeof(entry_t) == sizeof(index_t), "size of entry_t is not the same as index_t");

entry_t make_entry(index_t _index, coefficient_t _coefficient) {
	return entry_t(_index, _coefficient);
}
index_t get_index(const entry_t& e) { return e.index; }
index_t get_coefficient(const entry_t& e) { return e.coefficient; }
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
const index_t get_index(const entry_t& i) { return i; }
index_t get_coefficient(const entry_t& i) { return 1; }
entry_t make_entry(index_t _index, coefficient_t _value) { return entry_t(_index); }
void set_coefficient(entry_t& e, const coefficient_t c) {}

#endif

const entry_t& get_entry(const entry_t& e) { return e; }

typedef std::pair<value_t, index_t> diameter_index_t;
value_t get_diameter(const diameter_index_t& i) { return i.first; }
index_t get_index(const diameter_index_t& i) { return i.second; }

class diameter_entry_t : public std::pair<value_t, entry_t> {
public:
	diameter_entry_t() {}
	diameter_entry_t(const entry_t& e) : std::pair<value_t, entry_t>(0, e) {}
	diameter_entry_t(value_t _diameter, index_t _index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(_diameter, make_entry(_index, _coefficient)) {}
	diameter_entry_t(const diameter_index_t& _diameter_index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(get_diameter(_diameter_index),
	                                  make_entry(get_index(_diameter_index), _coefficient)) {}
	diameter_entry_t(const diameter_index_t& _diameter_index)
	    : diameter_entry_t(_diameter_index, 1) {}
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
	bool operator()(const Entry& a, const Entry& b) {
		return (get_diameter(a) > get_diameter(b)) ||
		       ((get_diameter(a) == get_diameter(b)) && (get_index(a) < get_index(b)));
	}
};

enum compressed_matrix_layout { LOWER_TRIANGULAR, UPPER_TRIANGULAR };

template <compressed_matrix_layout Layout> class compressed_distance_matrix {
public:
	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	void init_rows();

	compressed_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(std::move(_distances)), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
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

class sparse_distance_matrix {
public:
	std::vector<std::vector<diameter_index_t>> neighbors;

	//Initialize from thresholded dense distance matrix
	template <typename DistanceMatrix>
	sparse_distance_matrix(const DistanceMatrix& mat, value_t threshold) : neighbors(mat.size()) {
		for (index_t i = 0; i < size(); ++i) {
			for (index_t j = 0; j < size(); ++j) {
				if (i != j && mat(i, j) <= threshold) {
					neighbors[i].push_back(std::make_pair(mat(i, j), j));
				}
			}
		}
	}
	//Initialize from COO format
	sparse_distance_matrix(int* I, int* J, float* V, int NEdges, int N, float threshold) : neighbors(N) {
		int i, j;
		value_t val;
		for (int idx = 0; idx < NEdges; idx++) {
			i = I[idx];
			j = J[idx];
			val = V[idx];
			if (i < j && val <= threshold) {
				neighbors[i].push_back(std::make_pair(val, j));
				neighbors[j].push_back(std::make_pair(val, i));
			}
		}
	}

	size_t size() const { return neighbors.size(); }
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

template <>
value_t compressed_distance_matrix<UPPER_TRIANGULAR>::operator()(const index_t i,
                                                                 const index_t j) const {
	return i == j ? 0 : i > j ? rows[j][i] : rows[i][j];
}

template <>
value_t compressed_distance_matrix<LOWER_TRIANGULAR>::operator()(const index_t i,
                                                                 const index_t j) const {
	return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
}

typedef compressed_distance_matrix<LOWER_TRIANGULAR> compressed_lower_distance_matrix;
typedef compressed_distance_matrix<UPPER_TRIANGULAR> compressed_upper_distance_matrix;


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

template <typename Heap>
void push_entry(Heap& column, index_t i, coefficient_t c, value_t diameter) {
	entry_t e = make_entry(i, c);
	column.push(std::make_pair(diameter, e));
}

template <typename DistanceMatrix> class ripser {
	DistanceMatrix dist;
	index_t n, dim_max;
	value_t threshold;
	int do_cocycles; //If this flag is off, don't extract the representative cocycles to save time
	float ratio;
	coefficient_t modulus;
	const binomial_coeff_table binomial_coeff;
	std::vector<coefficient_t> multiplicative_inverse;
	mutable std::vector<index_t> vertices;
	mutable std::vector<std::vector<diameter_index_t>::const_reverse_iterator> neighbor_it;
	mutable std::vector<std::vector<diameter_index_t>::const_reverse_iterator> neighbor_end;
	mutable std::vector<diameter_entry_t> coface_entries;


public:
	/*Store all of the results into an array with the following format
	 [NumperClass0, PD0,
	  NumPerClass1, PD1, Cocycle1_0Len, Cocycle1_0, Cocycle1_1Len, Cocycle1_1, ..., Cocycle1_NumPerClass1Len, Cocycle1_NumPerClass1
	  NumPerClass2, PD2, Cocycle1_0Len, Cocycle1_0, Cocycle1_1Len, Cocycle1_1, ..., Cocycle2_NumPerClass1Len, Cocycle2_NumPerClass1
	  ...
	  , n_edges]

	where PD0 is b0, d0, b1, d1, ...
	and each cocycle is simplexidx0, simplexidx1, ..., simplexidxd, val mod p
	*/
	mutable std::vector<value_t> retvec;

	ripser(DistanceMatrix&& _dist, index_t _dim_max, value_t _threshold, float _ratio,
	       coefficient_t _modulus, int _do_cocycles)
	    : dist(std::move(_dist)), n(dist.size()),
	      dim_max(std::min(_dim_max, index_t(dist.size() - 2))), threshold(_threshold),
	      ratio(_ratio), modulus(_modulus), do_cocycles(_do_cocycles),
		  binomial_coeff(n, dim_max + 2),
	      multiplicative_inverse(multiplicative_inverse_vector(_modulus)) {}

	index_t get_next_vertex(index_t& v, const index_t idx, const index_t k) const {
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
		assert(binomial_coeff(v, k) <= idx && binomial_coeff(v + 1, k) > idx);
		return v;
	}

	index_t get_edge_index(const index_t i, const index_t j) const {
		return binomial_coeff(i, 2) + j;
	}

	template <typename OutputIterator>
	OutputIterator get_simplex_vertices(index_t idx, const index_t dim, index_t v,
	                                    OutputIterator out) const {
		--v;
		for (index_t k = dim + 1; k > 0; --k) {
			get_next_vertex(v, idx, k);
			*out++ = v;
			idx -= binomial_coeff(v, k);
		}
		return out;
	}

	value_t compute_diameter(const index_t index, index_t dim) const {
		value_t diam = -std::numeric_limits<value_t>::infinity();

		vertices.clear();
		get_simplex_vertices(index, dim, dist.size(), std::back_inserter(vertices));

		for (index_t i = 0; i <= dim; ++i)
			for (index_t j = 0; j < i; ++j) {
				diam = std::max(diam, dist(vertices[i], vertices[j]));
			}
		return diam;
	}

	class simplex_coboundary_enumerator;

	void assemble_columns_to_reduce(std::vector<diameter_index_t>& simplices,
	                                std::vector<diameter_index_t>& columns_to_reduce,
	                                hash_map<index_t, index_t>& pivot_column_index, index_t dim);

	void compute_dim_0_pairs(std::vector<diameter_index_t>& edges,
	                         std::vector<diameter_index_t>& columns_to_reduce) {
		std::vector<value_t> birthsanddeaths;
		//TODO: Get correct birth times if the edges are negative (required for lower star)
		union_find dset(n);

		edges = get_edges();

		std::sort(edges.rbegin(), edges.rend(),
		          greater_diameter_or_smaller_index<diameter_index_t>());

		std::vector<index_t> vertices_of_edge(2);
		for (auto e : edges) {
			vertices_of_edge.clear();
			get_simplex_vertices(get_index(e), 1, n, std::back_inserter(vertices_of_edge));
			index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

			if (u != v) {
				if (get_diameter(e) != 0) {
					birthsanddeaths.push_back(0);
					birthsanddeaths.push_back((value_t)get_diameter(e));
				}
				dset.link(u, v);
			} else
				columns_to_reduce.push_back(e);
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

		for (index_t i = 0; i < n; ++i) {
			if (dset.find(i) == i) {
				birthsanddeaths.push_back(0);
				birthsanddeaths.push_back(std::numeric_limits<value_t>::infinity());
			}
		}
		retvec.push_back(birthsanddeaths.size()/2);
		retvec.insert(retvec.end(), birthsanddeaths.begin(), birthsanddeaths.end());
	}

	template <typename Column, typename Iterator>
	diameter_entry_t add_coboundary_and_get_pivot(Iterator column_begin, Iterator column_end,
	                                              coefficient_t factor_column_to_add,
#ifdef ASSEMBLE_REDUCTION_MATRIX
	                                              Column& working_reduction_column,
#endif
	                                              Column& working_coboundary, const index_t& dim,
	                                              hash_map<index_t, index_t>& pivot_column_index,
	                                              bool& might_be_apparent_pair) {
		for (auto it = column_begin; it != column_end; ++it) {
			diameter_entry_t simplex = *it;
			set_coefficient(simplex, get_coefficient(simplex) * factor_column_to_add % modulus);

#ifdef ASSEMBLE_REDUCTION_MATRIX
			working_reduction_column.push(simplex);
#endif

			coface_entries.clear();
			simplex_coboundary_enumerator cofaces(simplex, dim, *this);
			while (cofaces.has_next()) {
				diameter_entry_t coface = cofaces.next();
				if (get_diameter(coface) <= threshold) {
					coface_entries.push_back(coface);
					if (might_be_apparent_pair && (get_diameter(simplex) == get_diameter(coface))) {
						if (pivot_column_index.find(get_index(coface)) ==
						    pivot_column_index.end()) {
							return coface;
						}
						might_be_apparent_pair = false;
					}
				}
			}
			for (auto coface : coface_entries) working_coboundary.push(coface);
		}

		return get_pivot(working_coboundary, modulus);
	}

	void compute_pairs(std::vector<diameter_index_t>& columns_to_reduce,
	                   hash_map<index_t, index_t>& pivot_column_index, index_t dim) {
#ifdef ASSEMBLE_REDUCTION_MATRIX
		compressed_sparse_matrix<diameter_entry_t> reduction_matrix;
#else
#ifdef USE_COEFFICIENTS
		std::vector<diameter_entry_t> reduction_matrix;
#endif
#endif

		std::vector<diameter_entry_t> coface_entries;
		std::vector<value_t> birthsanddeaths;
		std::vector<value_t> allcocycles;

		for (index_t index_column_to_reduce = 0; index_column_to_reduce < columns_to_reduce.size();
		     ++index_column_to_reduce) {
			auto column_to_reduce = columns_to_reduce[index_column_to_reduce];

#ifdef ASSEMBLE_REDUCTION_MATRIX
			std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>,
			                    greater_diameter_or_smaller_index<diameter_entry_t>>
			    working_reduction_column;
#endif

			std::priority_queue<diameter_entry_t, std::vector<diameter_entry_t>,
			                    greater_diameter_or_smaller_index<diameter_entry_t>>
			    working_coboundary;

			value_t diameter = get_diameter(column_to_reduce);

			index_t index_column_to_add = index_column_to_reduce;

			diameter_entry_t pivot;

			// start with factor 1 in order to initialize working_coboundary
			// with the coboundary of the simplex with index column_to_reduce
			coefficient_t factor_column_to_add = 1;

#ifdef ASSEMBLE_REDUCTION_MATRIX
			// initialize reduction_matrix as identity matrix
			reduction_matrix.append_column();
#endif
#ifdef USE_COEFFICIENTS
			reduction_matrix.push_back(diameter_entry_t(column_to_reduce, 1));
#endif

			bool might_be_apparent_pair = (index_column_to_reduce == index_column_to_add);

			while (true) {
#ifdef ASSEMBLE_REDUCTION_MATRIX
#ifdef USE_COEFFICIENTS
				auto reduction_column_begin = reduction_matrix.cbegin(index_column_to_add),
				     reduction_column_end = reduction_matrix.cend(index_column_to_add);
#else
				std::vector<diameter_entry_t> coeffs;
				coeffs.push_back(columns_to_reduce[index_column_to_add]);
				for (auto it = reduction_matrix.cbegin(index_column_to_add);
				     it != reduction_matrix.cend(index_column_to_add); ++it)
					coeffs.push_back(*it);
				auto reduction_column_begin = coeffs.begin(), reduction_column_end = coeffs.end();
#endif
#else
#ifdef USE_COEFFICIENTS
				auto reduction_column_begin = &reduction_matrix[index_column_to_add],
				     reduction_column_end = &reduction_matrix[index_column_to_add] + 1;
#else
				auto reduction_column_begin = &columns_to_reduce[index_column_to_add],
				     reduction_column_end = &columns_to_reduce[index_column_to_add] + 1;
#endif
#endif

				pivot = add_coboundary_and_get_pivot(
				    reduction_column_begin, reduction_column_end, factor_column_to_add,
#ifdef ASSEMBLE_REDUCTION_MATRIX
				    working_reduction_column,
#endif
				    working_coboundary, dim, pivot_column_index, might_be_apparent_pair);

				if (get_index(pivot) != -1) {
					auto pair = pivot_column_index.find(get_index(pivot));

					if (pair != pivot_column_index.end()) {
						index_column_to_add = pair->second;
						factor_column_to_add = modulus - get_coefficient(pivot);
					} else {
						value_t death = get_diameter(pivot);
						if (diameter != death) {
							birthsanddeaths.push_back(diameter);
							birthsanddeaths.push_back(death);

#ifdef ASSEMBLE_REDUCTION_MATRIX
							if (do_cocycles) {
								std::vector<value_t> thiscocycle;
<<<<<<< HEAD
								index_t lencocycle = 0;
								while (get_index(e = get_pivot(cocycle, modulus)) != -1) {
									get_simplex_vertices(get_index(e), dim, n, std::back_inserter(simplex));
									for (size_t k = 0; k < simplex.size(); k++) {
										thiscocycle.push_back((value_t)simplex[k]);
									}
									thiscocycle.push_back(normalize(get_coefficient(e), modulus));
									cocycle.pop();
									lencocycle++;
								}
								allcocycles.push_back(lencocycle);
=======
								index_t lencocycle = extract_cocycle(thiscocycle, working_reduction_column, dim);
								allcocycles.push_back((value_t)lencocycle);
>>>>>>> 631fe3378581c5cf8e033a5b57146dc43d4743ae
								allcocycles.insert(allcocycles.end(), thiscocycle.begin(), thiscocycle.end());
							}
#endif

						}
						pivot_column_index.insert(
						    std::make_pair(get_index(pivot), index_column_to_reduce));

#ifdef USE_COEFFICIENTS
						const coefficient_t inverse =
						    multiplicative_inverse[get_coefficient(pivot)];
#endif

#ifdef ASSEMBLE_REDUCTION_MATRIX
						// replace current column of reduction_matrix (with a single diagonal 1
						// entry) by reduction_column (possibly with a different entry on the
						// diagonal)
#ifdef USE_COEFFICIENTS
						reduction_matrix.pop_back();
#else
						pop_pivot(working_reduction_column, modulus);
#endif

						while (true) {
							diameter_entry_t e = pop_pivot(working_reduction_column, modulus);
							if (get_index(e) == -1) break;
#ifdef USE_COEFFICIENTS
							set_coefficient(e, inverse * get_coefficient(e) % modulus);
							assert(get_coefficient(e) > 0);
#endif
							reduction_matrix.push_back(e);
						}
#else
#ifdef USE_COEFFICIENTS
						reduction_matrix.pop_back();
						reduction_matrix.push_back(diameter_entry_t(column_to_reduce, inverse));
#endif
#endif
						break;
					}
				} else {
					birthsanddeaths.push_back(diameter);
					birthsanddeaths.push_back(std::numeric_limits<value_t>::infinity());


#ifdef ASSEMBLE_REDUCTION_MATRIX
					if (do_cocycles) {
						//Representative cocycle
						std::vector<value_t> thiscocycle;
<<<<<<< HEAD
						index_t lencocycle = 0;
						while (get_index(e = get_pivot(cocycle, modulus)) != -1) {
							get_simplex_vertices(get_index(e), dim, n, std::back_inserter(simplex));
							for (size_t k = 0; k < simplex.size(); k++) {
								thiscocycle.push_back((value_t)simplex[k]);
							}
							thiscocycle.push_back(normalize(get_coefficient(e), modulus));
							cocycle.pop();
							lencocycle++;
						}
						allcocycles.push_back(lencocycle);
=======
						index_t lencocycle = extract_cocycle(thiscocycle, working_reduction_column, dim);
						allcocycles.push_back((value_t)lencocycle);
>>>>>>> 631fe3378581c5cf8e033a5b57146dc43d4743ae
						allcocycles.insert(allcocycles.end(), thiscocycle.begin(), thiscocycle.end());
					}
#endif

					break;
				}
			}
		}

		//Now copy everything into the result array
		//Copy over persistence diagram
		retvec.push_back(birthsanddeaths.size()/2);
		retvec.insert(retvec.end(), birthsanddeaths.begin(), birthsanddeaths.end());
		//Copy over cocycles
		retvec.insert(retvec.end(), allcocycles.begin(), allcocycles.end());
	}

	std::vector<diameter_index_t> get_edges();

	void compute_barcodes() {

		std::vector<diameter_index_t> simplices, columns_to_reduce;

		compute_dim_0_pairs(simplices, columns_to_reduce);

		for (index_t dim = 1; dim <= dim_max; ++dim) {
			hash_map<index_t, index_t> pivot_column_index;
			pivot_column_index.reserve(columns_to_reduce.size());

			compute_pairs(columns_to_reduce, pivot_column_index, dim);

			if (dim < dim_max) {
				assemble_columns_to_reduce(simplices, columns_to_reduce, pivot_column_index,
				                           dim + 1);
			}
		}
	}

	template<typename Chain>
	value_t extract_cocycle(std::vector<value_t>& thiscocycle, Chain& cycle, index_t dim) {
		diameter_entry_t e;
		value_t lencocycle = 0;
		while (get_index(e = get_pivot(cycle, modulus)) != -1) {
			vertices.resize(dim + 1);
			get_simplex_vertices(get_index(e), dim, n,
								 vertices.rbegin());
			
			auto it = vertices.begin();
			if (it != vertices.end()) {
				thiscocycle.push_back((value_t)(*it++));
				while (it != vertices.end()) thiscocycle.push_back((value_t)(*it++));
			}
			thiscocycle.push_back((value_t)(normalize(get_coefficient(e), modulus)));
			cycle.pop();
			lencocycle++;
		}
		return lencocycle;
	}
};

template <> class ripser<compressed_lower_distance_matrix>::simplex_coboundary_enumerator {
private:
	index_t idx_below, idx_above, v, k;
	std::vector<index_t> vertices;
	const diameter_entry_t simplex;
	const coefficient_t modulus;
	const compressed_lower_distance_matrix& dist;
	const binomial_coeff_table& binomial_coeff;

public:
	simplex_coboundary_enumerator(const diameter_entry_t _simplex, index_t _dim,
	                              const ripser& parent)
	    : idx_below(get_index(_simplex)), idx_above(0), v(parent.n - 1), k(_dim + 1),
	      vertices(_dim + 1), simplex(_simplex), modulus(parent.modulus), dist(parent.dist),
	      binomial_coeff(parent.binomial_coeff) {
		parent.get_simplex_vertices(get_index(_simplex), _dim, parent.n, vertices.begin());
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

	diameter_entry_t next() {
		value_t coface_diameter = get_diameter(simplex);
		for (index_t w : vertices) coface_diameter = std::max(coface_diameter, dist(v, w));
		index_t coface_index = idx_above + binomial_coeff(v--, k + 1) + idx_below;
		coefficient_t coface_coefficient =
		    (k & 1 ? -1 + modulus : 1) * get_coefficient(simplex) % modulus;
		return diameter_entry_t(coface_diameter, coface_index, coface_coefficient);
	}
};

template <> class ripser<sparse_distance_matrix>::simplex_coboundary_enumerator {
private:
	const ripser& parent;

	index_t idx_below, idx_above, v, k, max_vertex_below;
	const diameter_entry_t simplex;
	const coefficient_t modulus;
	const sparse_distance_matrix& dist;
	const binomial_coeff_table& binomial_coeff;

	std::vector<index_t>& vertices;
	std::vector<std::vector<diameter_index_t>::const_reverse_iterator>& neighbor_it;
	std::vector<std::vector<diameter_index_t>::const_reverse_iterator>& neighbor_end;
	diameter_index_t x;

public:
	simplex_coboundary_enumerator(const diameter_entry_t _simplex, index_t _dim,
	                              const ripser& _parent)
	    : parent(_parent), idx_below(get_index(_simplex)), idx_above(0), v(parent.n - 1),
	      k(_dim + 1), max_vertex_below(parent.n - 1), simplex(_simplex), modulus(parent.modulus),
	      dist(parent.dist), binomial_coeff(parent.binomial_coeff), vertices(parent.vertices),
	      neighbor_it(parent.neighbor_it), neighbor_end(parent.neighbor_end) {

		neighbor_it.clear();
		neighbor_end.clear();
		vertices.clear();

		parent.get_simplex_vertices(idx_below, _dim, parent.n, std::back_inserter(vertices));

		for (auto v : vertices) {
			neighbor_it.push_back(dist.neighbors[v].rbegin());
			neighbor_end.push_back(dist.neighbors[v].rend());
		}
	}

	bool has_next(bool all_cofaces = true) {
		for (auto &it0 = neighbor_it[0], &end0 = neighbor_end[0]; it0 != end0; ++it0) {
			x = *it0;
			for (size_t idx = 1; idx < neighbor_it.size(); ++idx) {
				auto &it = neighbor_it[idx], end = neighbor_end[idx];
				while (get_index(*it) > get_index(x))
					if (++it == end) return false;
				auto y = *it;
				if (get_index(y) != get_index(x))
					goto continue_outer;
				else
					x = std::max(x, y);
			}
			return all_cofaces || !(k > 0 && parent.get_next_vertex(max_vertex_below, idx_below,
			                                                        k) > get_index(x));
		continue_outer:;
		}
		return false;
	}

	diameter_entry_t next() {
		++neighbor_it[0];

		while (k > 0 && parent.get_next_vertex(max_vertex_below, idx_below, k) > get_index(x)) {
			idx_below -= binomial_coeff(max_vertex_below, k);
			idx_above += binomial_coeff(max_vertex_below, k + 1);
			--k;
		}

		value_t coface_diameter = std::max(get_diameter(simplex), get_diameter(x));

		coefficient_t coface_coefficient =
		    (k & 1 ? -1 + modulus : 1) * get_coefficient(simplex) % modulus;

		return diameter_entry_t(coface_diameter,
		                        idx_above + binomial_coeff(get_index(x), k + 1) + idx_below,
		                        coface_coefficient);
	}
};

template <> std::vector<diameter_index_t> ripser<compressed_lower_distance_matrix>::get_edges() {
	std::vector<diameter_index_t> edges;
	for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
		value_t diameter = compute_diameter(index, 1);
		if (diameter <= threshold) edges.push_back(std::make_pair(diameter, index));
	}
	return edges;
}

template <> std::vector<diameter_index_t> ripser<sparse_distance_matrix>::get_edges() {
	std::vector<diameter_index_t> edges;
	for (index_t i = 0; i < n; ++i)
		for (auto n : dist.neighbors[i]) {
			index_t j = get_index(n);
			if (i > j) edges.push_back(std::make_pair(get_diameter(n), get_edge_index(i, j)));
		}
	return edges;
}

template <>
void ripser<compressed_lower_distance_matrix>::assemble_columns_to_reduce(
    std::vector<diameter_index_t>& simplices, std::vector<diameter_index_t>& columns_to_reduce,
    hash_map<index_t, index_t>& pivot_column_index, index_t dim) {
	index_t num_simplices = binomial_coeff(n, dim + 1);

	columns_to_reduce.clear();

	for (index_t index = 0; index < num_simplices; ++index) {
		if (pivot_column_index.find(index) == pivot_column_index.end()) {
			value_t diameter = compute_diameter(index, dim);
			if (diameter <= threshold) columns_to_reduce.push_back(std::make_pair(diameter, index));
		}
	}

	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
	          greater_diameter_or_smaller_index<diameter_index_t>());
}

template <>
void ripser<sparse_distance_matrix>::assemble_columns_to_reduce(
    std::vector<diameter_index_t>& simplices, std::vector<diameter_index_t>& columns_to_reduce,
    hash_map<index_t, index_t>& pivot_column_index, index_t dim) {

	--dim;
	columns_to_reduce.clear();

	std::vector<diameter_index_t> next_simplices;

	for (diameter_index_t simplex : simplices) {
		simplex_coboundary_enumerator cofaces(simplex, dim, *this);

		while (cofaces.has_next(false)) {
			auto coface = cofaces.next();

			next_simplices.push_back(std::make_pair(get_diameter(coface), get_index(coface)));

			if (pivot_column_index.find(get_index(coface)) == pivot_column_index.end())
				columns_to_reduce.push_back(
				    std::make_pair(get_diameter(coface), get_index(coface)));
		}
	}

	simplices.swap(next_simplices);
}

std::vector<value_t> pythondm(float* D, int N, int modulus, int dim_max, float threshold, int do_cocycles) {
	//Setup distance matrix and figure out threshold
	std::vector<value_t> retvec;
	std::vector<value_t> distances(D, D+N);
	compressed_lower_distance_matrix dist =
		compressed_lower_distance_matrix(compressed_upper_distance_matrix(std::move(distances)));
	index_t n = dist.size();
	dim_max = std::min((index_t)(dim_max), n-2);
	float ratio = 1.0; //TODO: This seems like a dummy parameter at the moment

	value_t min = std::numeric_limits<value_t>::infinity(),
	        max = -std::numeric_limits<value_t>::infinity(), max_finite = max;
	int num_edges = 0;

	value_t enclosing_radius = std::numeric_limits<value_t>::infinity();
	for (index_t i = 0; i < dist.size(); ++i) {
		value_t r_i = -std::numeric_limits<value_t>::infinity();
		for (index_t j = 0; j < dist.size(); ++j) r_i = std::max(r_i, dist(i, j));
		enclosing_radius = std::min(enclosing_radius, r_i);
	}

	if (threshold == std::numeric_limits<value_t>::max()) threshold = enclosing_radius;

	for (auto d : dist.distances) {
		min = std::min(min, d);
		max = std::max(max, d);
		max_finite = d != std::numeric_limits<value_t>::infinity() ? std::max(max, d) : max_finite;
		if (d <= threshold) ++num_edges;
	}


	if (threshold >= max) {
		ripser<compressed_lower_distance_matrix> r(std::move(dist), dim_max, threshold, ratio,
		                                         modulus, do_cocycles);
		r.compute_barcodes();
		retvec.insert(retvec.end(), r.retvec.begin(), r.retvec.end());
	} else {
		ripser<sparse_distance_matrix> r(sparse_distance_matrix(std::move(dist), threshold), dim_max,
		                               threshold, ratio, modulus, do_cocycles);
		r.compute_barcodes();
		retvec.insert(retvec.end(), r.retvec.begin(), r.retvec.end());
	}
	//Report the number of edges that were added
	retvec.push_back((value_t)num_edges);
	return retvec;
}


std::vector<value_t> pythondmsparse(int* I, int* J, float* V, int NEdges,
								 int N, int modulus, int dim_max, float threshold, int do_cocycles) {
	//Setup distance matrix and figure out threshold
	std::vector<value_t> retvec;
	float ratio = 1.0; //TODO: This seems like a dummy parameter at the moment
	ripser<sparse_distance_matrix> r(sparse_distance_matrix(I, J, V, NEdges, N, threshold),
									dim_max, threshold, ratio, modulus, do_cocycles);
	r.compute_barcodes();
	retvec.insert(retvec.end(), r.retvec.begin(), r.retvec.end());
	//Report the number of edges that were added
	value_t num_edges;
	for (int idx = 0; idx < NEdges; idx++) {
		if (I[idx] < J[idx] && V[idx] <= threshold) {
			num_edges++;
		}
	}
	retvec.push_back(num_edges);
	return retvec;
}

int unwrapvector(std::vector<value_t> vec, float** out) {
    int length = vec.size();
    float* arr = (float*) malloc(length * sizeof(value_t));
    std::memcpy(arr, &vec[0], length * sizeof(value_t));
    *out = arr;
    return length;
}

extern "C" {
    int cripser(float** out, float* D, int N,
                int modulus, int dim_max, float threshold, int do_cocycles) {
        std::vector<value_t> resvec = pythondm(D, N, modulus, dim_max, threshold, do_cocycles);
        return unwrapvector(resvec, out);
    }

    int cripsersparse(float** out, int* I, int* J, float* V, int NEdges, int N,
                      int modulus, int dim_max, float threshold, int do_cocycles) {
        std::vector<value_t> resvec = pythondmsparse(I, J, V, NEdges, N, modulus, dim_max, threshold, do_cocycles);
        return unwrapvector(resvec, out);
    }
}
