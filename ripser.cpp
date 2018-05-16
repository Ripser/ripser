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

#define ASSEMBLE_REDUCTION_MATRIX
//#define USE_COEFFICIENTS

//#define INDICATE_PROGRESS
#define PRINT_PERSISTENCE_PAIRS

//#define USE_GOOGLE_HASHMAP

//#define COFACE_BASED_FILTRATION
#define FACE_BASED_FILTRATION

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>

#ifdef USE_GOOGLE_HASHMAP
#include <sparsehash/sparse_hash_map>
template <class Key, class T> class hash_map : public google::sparse_hash_map<Key, T> {
public:
	inline void reserve(size_t hint) { this->resize(hint); }
};
#else
template <class Key, class T> class hash_map : public std::unordered_map<Key, T> {};
#endif

typedef float value_t;
// typedef uint16_t value_t;

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

template <typename Entry> struct smaller_index2 {
	bool operator()(Entry& a, Entry& b) { return a.get_index() < b.get_index(); }
};

class filtered_simplex {
	index_t index;
	coefficient_t coefficient;
	value_t filtration_value;

public:
	filtered_simplex(const filtered_simplex& s) {
		index = s.index;
		filtration_value = s.filtration_value;
	}
	filtered_simplex(index_t ind, value_t fval) {
		index = ind;
		filtration_value = fval;
	}
	filtered_simplex(index_t ind) {
		index = ind;
		filtration_value = 0;
	}
	value_t get_filtration_value() { return this->filtration_value; }
	index_t get_index() { return this->index; }
};

class filtered_simplex_coeff {
	index_t index;
	coefficient_t coefficient;
	value_t filtration_value;

public:
	filtered_simplex_coeff(const filtered_simplex_coeff& s) {
		index = s.index;
		coefficient = s.coefficient;
		filtration_value = s.filtration_value;
	}
	filtered_simplex_coeff(index_t ind, coefficient_t coeff, value_t fval) {
		index = ind;
		coefficient = coeff;
		filtration_value = fval;
	}
	filtered_simplex_coeff(index_t ind) {
		index = ind;
		coefficient = 1;
		filtration_value = 0;
	}
	filtered_simplex_coeff(filtered_simplex s, coefficient_t coeff) {
		this->index = s.get_index();
		this->filtration_value = s.get_filtration_value();
		this->coefficient = coeff;
	}
	coefficient_t get_coefficient() { return this->coefficient; }
	value_t get_filtration_value() { return this->filtration_value; }
	index_t get_index() { return this->index; }
	void set_coefficient(coefficient_t coeff) { this->coefficient = coeff; }
	void set_filtration(value_t fval) { this->filtration_value = fval; }
};

template <typename Entry> struct greater_filtration_or_smaller_index {
	bool operator()(Entry& a, Entry& b) {
		return (a.get_filtration_value() > b.get_filtration_value()) ||
		       ((a.get_filtration_value() == b.get_filtration_value()) && (a.get_index() < b.get_index()));
	}
};

template <typename Iterator> bool next_face(const Iterator first, Iterator k, const Iterator last);

bool descend_compare(index_t i, index_t j) { return (i > j); }

class simplex {
	std::vector<index_t> vertices;
	value_t filtration_value;
	index_t CNS_index;

public:
	void add_vertex(index_t num) { this->vertices.push_back(num); }

	void assign_filtration(std::string str) {
		std::stringstream stream(str);
		stream >> this->filtration_value;
	}

	void assign_vertices(std::vector<index_t>::iterator begin, std::vector<index_t>::iterator end) {
		vertices.assign(begin, end);
	}

	void assign_highest_facet(std::unordered_map<index_t, simplex>* simplicialFiltration,
	                          const binomial_coeff_table& B) {

		index_t facet_dim = this->vertices.size() - 2;
		value_t max = 0;

		if (facet_dim < 0) {
			filtration_value = 0;
			return;
		}
		do {
			simplex facet;
			facet.clear();
			facet.assign_vertices(vertices.begin(), vertices.begin() + facet_dim + 1);
			facet.sort();
			// face.assign_filtration(line.substr(string_end+1,line.length()));
			index_t index = facet.compute_index(B);

			if (simplicialFiltration[facet_dim].find(index) != simplicialFiltration[facet_dim].end()) {
				std::unordered_map<index_t, simplex>::const_iterator S = simplicialFiltration[facet_dim].find(index);
				value_t fval = S->second.get_filtration();
				max = max < fval ? fval : max;
			} else {
				std::cout << "Unexpected outcome: facet unassigned" << std::endl;
				exit(-1);
			}

		} while (next_face(vertices.begin(), vertices.begin() + facet_dim + 1, vertices.end()));
		filtration_value = max;
	}

	void assign_lowest_cofacet(std::unordered_map<index_t, simplex>* simplicialFiltration,
	                           const binomial_coeff_table& B, index_t n) {

		index_t min = std::numeric_limits<value_t>::max();
		std::vector<index_t> cofacet_vertices;
		simplex cofacet;
		index_t cofacet_dim = this->vertices.size();

		for (index_t i = 0; i < n; ++i) {
			if (!(std::binary_search(vertices.begin(), vertices.end(), i, descend_compare))) {
				cofacet_vertices.assign(vertices.begin(), vertices.end());
				cofacet_vertices.push_back(i);
				cofacet.clear();
				cofacet.assign_vertices(cofacet_vertices.begin(), cofacet_vertices.end());
				index_t index = cofacet.compute_index(B);
				if (simplicialFiltration[cofacet_dim].find(index) != simplicialFiltration[cofacet_dim].end()) {
					std::unordered_map<index_t, simplex>::const_iterator S =
					    simplicialFiltration[cofacet_dim].find(index);
					value_t fval = S->second.get_filtration();
					min = min > fval ? fval : min;
				}
			}
		}
		filtration_value = min;
	}

	void clear() {
		this->vertices.clear();
		this->filtration_value = 0;
	}
	index_t compute_index(const binomial_coeff_table& B) {

		index_t index = 0, j = this->vertices.size() - 1;
		for (index_t i : vertices) {
			index += B(i, j + 1);
			j--;
		}
		this->CNS_index = index;
		return index;
	}
	void sort() { std::sort(this->vertices.begin(), this->vertices.end(), descend_compare); }
	value_t get_filtration() const { return filtration_value; }
	index_t get_index() const { return CNS_index; }
	std::vector<index_t>::iterator vertices_begin() { return vertices.begin(); }
	std::vector<index_t>::iterator vertices_end() { return vertices.end(); }
	void set_index(index_t ind) { CNS_index = ind; }
	void set_filtration(value_t fval) { filtration_value = fval; }
};

class simplex_coboundary_enumerator2 {
private:
	index_t idx_below, idx_above, v, k;
	std::vector<index_t> vertices;
	filtered_simplex_coeff simplex;
	const coefficient_t modulus;
	const binomial_coeff_table& binomial_coeff;

public:
	simplex_coboundary_enumerator2(filtered_simplex_coeff _simplex, index_t _dim, index_t _n,
	                               const coefficient_t _modulus, const binomial_coeff_table& _binomial_coeff)
	    : simplex(_simplex), idx_below(_simplex.get_index()), idx_above(0), v(_n - 1), k(_dim + 1), modulus(_modulus),
	      binomial_coeff(_binomial_coeff), vertices(_dim + 1) {
		get_simplex_vertices(_simplex.get_index(), _dim, _n, binomial_coeff, vertices.begin());
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

	filtered_simplex_coeff next() {
		index_t coface_index = idx_above + binomial_coeff(v--, k + 1) + idx_below;
		value_t coface_fvalue = 1; // dummy
		coefficient_t coface_coefficient = (k & 1 ? -1 + modulus : 1) * simplex.get_coefficient() % modulus;
		return filtered_simplex_coeff(coface_index, coface_coefficient, coface_fvalue);
	}
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

template <typename Heap> filtered_simplex_coeff pop_pivot2(Heap& column, coefficient_t modulus) {
	if (column.empty())
		return filtered_simplex_coeff(-1);
	else {
		auto pivot = column.top();

#ifdef USE_COEFFICIENTS
		coefficient_t coefficient = 0;
		do {
			coefficient = (coefficient + get_coefficient(column.top())) % modulus;
			column.pop();

			if (coefficient == 0) {
				if (column.empty())
					return filtered_simplex_coeff(-1);
				else
					pivot = column.top();
			}
		} while (!column.empty() && get_index(column.top()) == get_index(pivot));
		if (get_index(pivot) != -1) { pivot.set_coefficient(coefficient); }
#else
		column.pop();
		auto column_top = column.top();
		while (!column.empty() && column_top.get_index() == pivot.get_index()) {
			column.pop();
			if (column.empty())
				return filtered_simplex_coeff(-1);
			else {
				pivot = column.top();
				column.pop();
			}
			column_top = column.top();
		}
#endif
		return pivot;
	}
}

template <typename Heap> filtered_simplex_coeff get_pivot2(Heap& column, coefficient_t modulus) {
	filtered_simplex_coeff result = pop_pivot2(column, modulus);
	if (result.get_index() != -1) column.push(result);
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

template <typename Heap> void push_entry2(Heap& column, index_t i, coefficient_t c, value_t fvalue) {
	filtered_simplex_coeff S(i, c, fvalue);
	column.push(S);
}

void assemble_columns_to_reduce2(std::vector<filtered_simplex>& columns_to_reduce,
                                 hash_map<index_t, index_t>& pivot_column_index,
                                 std::unordered_map<index_t, simplex>* simplicialFiltration, index_t dim, index_t n,
                                 value_t threshold, const binomial_coeff_table& binomial_coeff) {
	index_t num_simplices = binomial_coeff(n, dim + 2);

	columns_to_reduce.clear();

#ifdef INDICATE_PROGRESS
	std::cout << "\033[K"
	          << "assembling " << num_simplices << " columns" << std::flush << "\r";
#endif

	for (index_t index = 0; index < num_simplices; ++index) {
		if (pivot_column_index.find(index) == pivot_column_index.end() &&
		    simplicialFiltration[dim + 1].find(index) != simplicialFiltration[dim + 1].end()) {
			std::unordered_map<index_t, simplex>::const_iterator S = simplicialFiltration[dim + 1].find(index);
			value_t fval = S->second.get_filtration();
			index_t ind = index;
			filtered_simplex F(ind, fval);
			if (fval <= threshold) columns_to_reduce.push_back(F);
#ifdef INDICATE_PROGRESS
			if ((index + 1) % 1000 == 0)
				std::cout << "\033[K"
				          << "assembled " << columns_to_reduce.size() << " out of " << (index + 1) << "/"
				          << num_simplices << " columns" << std::flush << "\r";
#endif
		}
	}

#ifdef INDICATE_PROGRESS
	std::cout << "\033[K"
	          << "sorting " << num_simplices << " columns" << std::flush << "\r";
#endif

	std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
	          greater_filtration_or_smaller_index<filtered_simplex>());
#ifdef INDICATE_PROGRESS
	std::cout << "\033[K";
#endif
}

void compute_pairs2(std::vector<filtered_simplex>& columns_to_reduce, hash_map<index_t, index_t>& pivot_column_index,
                    index_t dim, index_t n, value_t threshold, coefficient_t modulus,
                    const std::vector<coefficient_t>& multiplicative_inverse,
                    std::unordered_map<index_t, simplex>* simplicialFiltration,
                    const binomial_coeff_table& binomial_coeff) {

#ifdef PRINT_PERSISTENCE_PAIRS
	std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
#endif

#ifdef ASSEMBLE_REDUCTION_MATRIX
	compressed_sparse_matrix<filtered_simplex_coeff> reduction_coefficients;
#else
#ifdef USE_COEFFICIENTS
	std::vector<filtered_simplex_coeff> reduction_coefficients;
#endif
#endif

	std::vector<filtered_simplex_coeff> coface_entries;

	for (index_t i = 0; i < columns_to_reduce.size(); ++i) {
		auto column_to_reduce = columns_to_reduce[i];

#ifdef ASSEMBLE_REDUCTION_MATRIX
		std::priority_queue<filtered_simplex_coeff, std::vector<filtered_simplex_coeff>,
		                    smaller_index2<filtered_simplex_coeff>>
		    reduction_column;
#endif

		std::priority_queue<filtered_simplex_coeff, std::vector<filtered_simplex_coeff>,
		                    greater_filtration_or_smaller_index<filtered_simplex_coeff>>
		    working_coboundary;

		value_t filtration_val = column_to_reduce.get_filtration_value();

#ifdef INDICATE_PROGRESS
		if ((i + 1) % 1000 == 0)
			std::cout << "\033[K"
			          << "reducing column " << i + 1 << "/" << columns_to_reduce.size()
			          << " (Filtration value:  " << filtration_val << ")" << std::flush << "\r";
#endif

		index_t j = i;

		// start with a dummy pivot entry with coefficient -1 in order to initialize
		// working_coboundary with the coboundary of the simplex with index column_to_reduce
		filtered_simplex_coeff pivot(-1, -1 + modulus, 0);

#ifdef ASSEMBLE_REDUCTION_MATRIX
		// initialize reduction_coefficients as identity matrix
		reduction_coefficients.append_column();
		reduction_coefficients.push_back(filtered_simplex_coeff(column_to_reduce, 1));
#else
#ifdef USE_COEFFICIENTS
		reduction_coefficients.push_back(filtered_simplex_coeff(column_to_reduce, 1));
#endif
#endif

		bool might_be_apparent_pair = (i == j);

		do {
			const coefficient_t factor = modulus - pivot.get_coefficient();

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
				filtered_simplex_coeff current_simplex = *it;
				current_simplex.set_coefficient(current_simplex.get_coefficient() * factor % modulus);

#ifdef ASSEMBLE_REDUCTION_MATRIX
				reduction_column.push(current_simplex);
#endif

				coface_entries.clear();
				simplex_coboundary_enumerator2 cofaces(current_simplex, dim, n, modulus, binomial_coeff);

				while (cofaces.has_next()) {
					filtered_simplex_coeff coface = cofaces.next();
					index_t ind = coface.get_index();
					if (simplicialFiltration[dim + 1].find(ind) != simplicialFiltration[dim + 1].end()) {
						std::unordered_map<index_t, simplex>::const_iterator S =
						    simplicialFiltration[dim + 1].find(ind);
						value_t fval = S->second.get_filtration();
						coface.set_filtration(fval);
						if (fval <= threshold) {
							coface_entries.push_back(coface);
							if (might_be_apparent_pair &&
							    (current_simplex.get_filtration_value() == coface.get_filtration_value())) {
								if (pivot_column_index.find(coface.get_index()) == pivot_column_index.end()) {
									pivot = coface;
									goto found_persistence_pair;
								}
								might_be_apparent_pair = false;
							}
						}
					}
				}
				for (auto e : coface_entries) working_coboundary.push(e);
			}

			pivot = get_pivot2(working_coboundary, modulus);

			if (pivot.get_index() != -1) {
				auto pair = pivot_column_index.find(pivot.get_index());

				if (pair != pivot_column_index.end()) {
					j = pair->second;
					continue;
				}
			} else {
#ifdef PRINT_PERSISTENCE_PAIRS
#ifdef INDICATE_PROGRESS
				std::cout << "\033[K";
#endif
				// std::cout << "Index "  << " [" << column_to_reduce.get_index() << ", )";
				std::cout << " [" << filtration_val << ", )" << std::endl << std::flush;
#endif
				break;
			}

		found_persistence_pair:
#ifdef PRINT_PERSISTENCE_PAIRS
			value_t death = pivot.get_filtration_value();
			if (filtration_val != death) {
#ifdef INDICATE_PROGRESS
				std::cout << "\033[K";
#endif
				//  std::cout << "Index "  << " [" << pivot.get_index() << "," << i << ")";
				std::cout << " [" << filtration_val << "," << death << ")" << std::endl << std::flush;
			}
#endif

			pivot_column_index.insert(std::make_pair(pivot.get_index(), i));

#ifdef USE_COEFFICIENTS
			const coefficient_t inverse = multiplicative_inverse[pivot.get_coefficient()];
#endif

#ifdef ASSEMBLE_REDUCTION_MATRIX
			// replace current column of reduction_coefficients (with a single diagonal 1 entry)
			// by reduction_column (possibly with a different entry on the diagonal)
			reduction_coefficients.pop_back();
			while (true) {
				filtered_simplex_coeff e = pop_pivot2(reduction_column, modulus);
				if (e.get_index() == -1) break;
#ifdef USE_COEFFICIENTS
				e.set_coefficient(inverse * e.get_coefficient() % modulus);
				assert(e.get_coefficient() > 0);
#endif
				reduction_coefficients.push_back(e);
			}
#else
#ifdef USE_COEFFICIENTS
			reduction_coefficients.pop_back();
			reduction_coefficients.push_back(filtered_simplex_coeff(column_to_reduce, inverse));
#endif
#endif
			break;
		} while (true);
	}

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
#ifdef USE_COEFFICIENTS
	          << "  --modulus <p>    compute homology with coefficients in the prime field Z/<p>Z"
#endif
	          << std::endl;

	exit(exit_code);
}

std::unordered_map<index_t, simplex>* read_file2(std::istream& input_stream, index_t& n, index_t& dim_max) {
	std::string line;
	std::string delimiter = "]";
	std::getline(input_stream, line);
	std::stringstream stream(line);
	index_t num;
	std::vector<index_t> new_simplex;
	std::size_t string_end;
	simplex current_simplex;
	stream >> n;
	stream >> dim_max;
	index_t dim;

	std::cout << "n " << n << std::endl;
	std::cout << "dim_max " << dim_max << std::endl;
	const binomial_coeff_table B(n, dim_max + 2);

	std::unordered_map<index_t, simplex>* simplicialFiltration = new std::unordered_map<index_t, simplex>[dim_max + 1];

	while (std::getline(input_stream, line)) {
		std::cout << line << std::endl;
		string_end = line.find(delimiter);
		//  std::cout << "substring " << line.substr(1,string_end-1) << std::endl;
		stream.str("");
		stream.clear();
		stream << line.substr(1, string_end - 1);
		// std::cout << line[line.length()] << "gotcha ";

		// std::cout << stream.str() << std::endl;
		current_simplex.clear();
		dim = 0;
		while (stream >> num) {
			current_simplex.add_vertex(num);
			// std::cout << num << std::endl;
			dim++;
		}
		dim--;
		current_simplex.sort();
		current_simplex.assign_filtration(line.substr(string_end + 1, line.length()));
		index_t index = current_simplex.compute_index(B);
		// std::cout << "index " << index << std::endl;
		std::pair<index_t, simplex> indexed_simplex;
		indexed_simplex.first = index;
		indexed_simplex.second = current_simplex;
		// std::cout << "dim " << dim << std::endl;
		simplicialFiltration[dim].insert(indexed_simplex);
	}
	return simplicialFiltration;
}

void missing_faces_helper_routine(std::unordered_map<index_t, simplex>* simplicialFiltration,
                                  std::vector<index_t> simplex_vertices, const binomial_coeff_table& B, index_t n,
                                  index_t face_dim) {
	simplex face;
	face.clear();
	face.assign_vertices(simplex_vertices.begin(), simplex_vertices.begin() + face_dim + 1);
	face.sort();
	// face.assign_filtration(line.substr(string_end+1,line.length()));

	index_t index = face.compute_index(B);

	if (simplicialFiltration[face_dim].find(index) == simplicialFiltration[face_dim].end()) {
		std::cout << "New face added" << std::endl;
		for (auto i = simplex_vertices.begin(); i != simplex_vertices.begin() + face_dim + 1; ++i)
			std::cout << *i << ' ';
		std::cout << std::endl;
		std::pair<index_t, simplex> indexed_simplex;
		indexed_simplex.second.clear();
#ifdef FACE_BASED_FILTRATION
		face.assign_highest_facet(simplicialFiltration, B);
#else
		face.assign_lowest_cofacet(simplicialFiltration, B, n);
#endif

		std::cout << "Filtration value: " << face.get_filtration() << std::endl;
		indexed_simplex.first = index;
		indexed_simplex.second = face;
		simplicialFiltration[face_dim].insert(indexed_simplex);
	}
}

template <typename Iterator> bool next_face(const Iterator first, Iterator k, const Iterator last) {
	if ((first == last) || (first == k) || (last == k)) {
		// std::cout << "First false " << std::endl;
		return false;
	}
	Iterator i1 = first;
	Iterator i2 = last;
	++i1;
	if (last == i1) {
		//     std::cout << "Second false " << std::endl;
		return false;
	}
	i1 = last;
	--i1;
	i1 = k;
	--i2;
	while (first != i1) {
		if (*--i1 < *i2) {
			Iterator j = k;
			while (!(*i1 < *j)) ++j;
			std::iter_swap(i1, j);
			++i1;
			++j;
			i2 = k;
			std::rotate(i1, j, last);
			while (last != j) {
				++j;
				++i2;
			}
			std::rotate(k, i2, last);
			return true;
		}
	}
	std::rotate(first, k, last);
	// std::cout << "Third false " << std::endl;

	return false;
}

void add_missing_faces(std::unordered_map<index_t, simplex>* simplicialFiltration, index_t dim_max, index_t n,
                       const binomial_coeff_table& B) {
/*

   hash_map<index_t, simplex>::iterator S = simplicialFiltration[0].begin();
  std::vector<index_t>::iterator a, b;
  a = S->second.vertices_begin(), b = S->second.vertices_end();
  std::vector<index_t> simplex_vertices ( a, b );

*/

#ifdef FACE_BASED_FILTRATION
	for (index_t dim = 1; dim <= dim_max; ++dim) {
		for (auto simplicialIterator = simplicialFiltration[dim].begin();
		     simplicialIterator != simplicialFiltration[dim].end(); ++simplicialIterator) {
			std::vector<index_t>::iterator a, b;
			a = simplicialIterator->second.vertices_begin(), b = simplicialIterator->second.vertices_end();
			std::vector<index_t> simplex_vertices(a, b);
			/*   for ( auto i = simplex_vertices.begin(); i != simplex_vertices.end() ; ++i)
			       std::cout << *i << ' ';
			   std::cout << std::endl; */

			for (index_t face_dim = 0; face_dim < dim; ++face_dim) {
				do {
					/*  std::cout << "face is "  << std::endl;
					  for ( auto i = simplex_vertices.begin(); i != simplex_vertices.begin()+ face_dim + 1 ; ++i)
					      std::cout << *i << ' ';
					  std::cout << std::endl; */
					missing_faces_helper_routine(simplicialFiltration, simplex_vertices, B, n, face_dim);

				} while (next_face(simplex_vertices.begin(), simplex_vertices.begin() + face_dim + 1,
				                   simplex_vertices.end()));
			}
		}
	}

#else
	for (index_t face_dim = dim_max - 1; face_dim >= 0; --face_dim) {
		for (index_t dim = dim_max; dim >= face_dim + 1; --dim) {
			for (auto simplicialIterator = simplicialFiltration[dim].begin();
			     simplicialIterator != simplicialFiltration[dim].end(); ++simplicialIterator) {
				std::vector<index_t>::iterator a, b;
				a = simplicialIterator->second.vertices_begin(), b = simplicialIterator->second.vertices_end();
				std::vector<index_t> simplex_vertices(a, b);
				/*  for ( auto i = simplex_vertices.begin(); i != simplex_vertices.begin()+ dim + 1 ; ++i)
				      std::cout << *i << ' ';
				  std::cout << std::endl; */

				do {
					/*     std::cout << "face is "  << std::endl;
					     for ( auto i = simplex_vertices.begin(); i != simplex_vertices.begin()+ face_dim + 1 ; ++i)
					         std::cout << *i << ' ';
					     std::cout << std::endl; */
					missing_faces_helper_routine(simplicialFiltration, simplex_vertices, B, n, face_dim);

				} while (next_face(simplex_vertices.begin(), simplex_vertices.begin() + face_dim + 1,
				                   simplex_vertices.end()));
			}
		}
	}

#endif
}

#include <iostream>

int main(int argc, const char* argv[]) {
	// insert code here...
	// std::cout << "We are here" << std::endl;

	const char* filename = nullptr;

	index_t n, dim_max = 1;
	value_t threshold = std::numeric_limits<value_t>::max();

#ifdef USE_COEFFICIENTS
	coefficient_t modulus = 2;
#else
	const coefficient_t modulus = 2;
#endif

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

	std::unordered_map<index_t, simplex>* simplicialFiltration;

	simplicialFiltration = read_file2(filename ? file_stream : std::cin, n, dim_max);

	std::cout << "Complex of dimension " << dim_max << " with " << n << " points" << std::endl;
	dim_max = std::min(dim_max, n - 2);

	binomial_coeff_table binomial_coeff(n, dim_max + 2);

	add_missing_faces(simplicialFiltration, dim_max, n, binomial_coeff);

	std::vector<coefficient_t> multiplicative_inverse(multiplicative_inverse_vector(modulus));

	std::vector<filtered_simplex> columns_to_reduce;

	{
		union_find dset(n);
		std::vector<filtered_simplex> edges;
		for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
			if (simplicialFiltration[1].find(index) != simplicialFiltration[1].end()) {
				hash_map<index_t, simplex>::const_iterator S = simplicialFiltration[1].find(index);
				value_t f_val = S->second.get_filtration();
				if (f_val <= threshold) edges.push_back(filtered_simplex(index, f_val));
			}
		}
		std::sort(edges.rbegin(), edges.rend(), greater_filtration_or_smaller_index<filtered_simplex>());

#ifdef PRINT_PERSISTENCE_PAIRS
		std::cout << "persistence intervals in dim 0:" << std::endl;
#endif

		std::vector<index_t> vertices_of_edge(2);
		for (auto e : edges) {
			vertices_of_edge.clear();
			get_simplex_vertices(e.get_index(), 1, n, binomial_coeff, std::back_inserter(vertices_of_edge));
			index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

			if (u != v) {
#ifdef PRINT_PERSISTENCE_PAIRS
				hash_map<index_t, simplex>::const_iterator v0 = simplicialFiltration[0].find(vertices_of_edge[0]);
				value_t f_val0 = v0->second.get_filtration();
				hash_map<index_t, simplex>::const_iterator v1 = simplicialFiltration[0].find(vertices_of_edge[1]);
				value_t f_val1 = v1->second.get_filtration();
				value_t paired_fval = f_val0 > f_val1 ? f_val0 : f_val1;
				if (e.get_filtration_value() > 0)
					std::cout << " [" << paired_fval << "," << e.get_filtration_value() << ")" << std::endl;
#endif
				dset.link(u, v);
			} else
				columns_to_reduce.push_back(e);
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

#ifdef PRINT_PERSISTENCE_PAIRS
		for (index_t i = 0; i < n; ++i)
			if (dset.find(i) == i) {
				// std::cout << " [0, )" << std::endl << std::flush;
				hash_map<index_t, simplex>::const_iterator unpaired_vertex = simplicialFiltration[0].find(i);
				value_t f_val = unpaired_vertex->second.get_filtration();
				std::cout << " [" << f_val << ","
				          << ")" << std::endl;
			}
#endif
	}

	for (index_t dim = 1; dim < dim_max; ++dim) {

		hash_map<index_t, index_t> pivot_column_index;
		pivot_column_index.reserve(columns_to_reduce.size());

		compute_pairs2(columns_to_reduce, pivot_column_index, dim, n, threshold, modulus, multiplicative_inverse,
		               simplicialFiltration, binomial_coeff);

		if (dim < dim_max - 1) {
			assemble_columns_to_reduce2(columns_to_reduce, pivot_column_index, simplicialFiltration, dim, n, threshold,
			                            binomial_coeff);
		}
	}

	std::cout << "Hello, World!\n";
	return 0;
}
