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

#ifdef USE_COEFFICIENTS
struct entry_t {
	index_t index : 8 * (sizeof(index_t) - sizeof(coefficient_t));
	coefficient_t coefficient;
	entry_t(index_t _index, coefficient_t _coefficient) : index(_index), coefficient(_coefficient) {}
	entry_t(index_t _index) : index(_index), coefficient(1) {}
	entry_t() : index(0), coefficient(1) {}
} __attribute__((packed));

static_assert(sizeof(entry_t) == sizeof(index_t), "size of entry_t is not the same as index_t");

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

index_t compute_index(const std::vector<index_t> vertices, const binomial_coeff_table& B) {
	
	index_t index = 0, j = vertices.size() - 1;
	for (index_t i : vertices) {
		index += B(i, j + 1);
		j--;
	}
	return index;
}

class simplex_coboundary_enumerator {
private:
	index_t idx_below, idx_above, v, k;
	std::vector<index_t> vertices;
	diameter_entry_t simplex;
	const coefficient_t modulus;
	const binomial_coeff_table& binomial_coeff;

public:
	simplex_coboundary_enumerator(diameter_entry_t _simplex, index_t _dim, index_t _n,
	                               const coefficient_t _modulus, const binomial_coeff_table& _binomial_coeff)
	    : simplex(_simplex), idx_below(get_index(_simplex)), idx_above(0), v(_n - 1), k(_dim + 1), modulus(_modulus),
	      binomial_coeff(_binomial_coeff), vertices(_dim + 1) {
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
		index_t coface_index = idx_above + binomial_coeff(v--, k + 1) + idx_below;
		value_t coface_fvalue = 1; // dummy
		coefficient_t coface_coefficient = (k & 1 ? -1 + modulus : 1) * get_coefficient(simplex) % modulus;
		return diameter_entry_t(coface_fvalue, coface_index, coface_coefficient);
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

void assemble_columns_to_reduce(std::vector<diameter_index_t>& columns_to_reduce,
                                 hash_map<index_t, index_t>& pivot_column_index,
                                 std::unordered_map<index_t, value_t>* simplicialFiltration, index_t dim, index_t n,
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
			auto S = simplicialFiltration[dim + 1].find(index);
			value_t fval = S->second;
			index_t ind = index;
			if (fval <= threshold) columns_to_reduce.push_back(std::make_pair(fval, index));
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
	          greater_diameter_or_smaller_index<diameter_index_t>());
#ifdef INDICATE_PROGRESS
	std::cout << "\033[K";
#endif
}

void compute_pairs(std::vector<diameter_index_t>& columns_to_reduce, hash_map<index_t, index_t>& pivot_column_index,
                    index_t dim, index_t n, value_t threshold, coefficient_t modulus,
                    const std::vector<coefficient_t>& multiplicative_inverse,
                    std::unordered_map<index_t, value_t>* simplicialFiltration,
                    const binomial_coeff_table& binomial_coeff) {

#ifdef PRINT_PERSISTENCE_PAIRS
	std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
#endif

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

#ifdef INDICATE_PROGRESS
		if ((i + 1) % 1000 == 0)
			std::cout << "\033[K"
			          << "reducing column " << i + 1 << "/" << columns_to_reduce.size()
			          << " (filtration value:  " << diameter << ")" << std::flush << "\r";
#endif

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
				simplex_coboundary_enumerator cofaces(simplex, dim, n, modulus, binomial_coeff);

				while (cofaces.has_next()) {
					diameter_entry_t coface = cofaces.next();
					auto S =
					simplicialFiltration[dim + 1].find(get_index(coface));
					if (S != simplicialFiltration[dim + 1].end()) {
						value_t fval = S->second;
						coface = std::make_pair(fval, get_entry(coface));
						if (fval <= threshold) {
							coface_entries.push_back(coface);
							if (might_be_apparent_pair &&
							    (get_diameter(simplex) == get_diameter(coface))) {
								if (pivot_column_index.find(get_index(coface)) == pivot_column_index.end()) {
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

			pivot = get_pivot(working_coboundary, modulus);

			if (get_index(pivot) != -1) {
				auto pair = pivot_column_index.find(get_index(pivot));

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
				std::cout << " [" << diameter << ", )" << std::endl << std::flush;
#endif
				break;
			}

		found_persistence_pair:
#ifdef PRINT_PERSISTENCE_PAIRS
			value_t death = get_diameter(pivot);
			if (diameter != death) {
#ifdef INDICATE_PROGRESS
				std::cout << "\033[K";
#endif
				//  std::cout << "Index "  << " [" << pivot.get_index() << "," << i << ")";
				std::cout << " [" << diameter << "," << death << ")" << std::endl << std::flush;
			}
#endif

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

std::unordered_map<index_t, value_t>* read_file(std::istream& input_stream, index_t& n, index_t& dim_max) {
	std::string line;
	std::string delimiter = "]";
	std::getline(input_stream, line);
	std::stringstream stream(line);
	index_t num;
	std::vector<index_t> new_simplex;
	std::size_t string_end;
	index_t dim;

	stream >> n;
	stream >> dim_max;

	std::cout << "n " << n << std::endl;
	std::cout << "dim_max " << dim_max << std::endl;
	const binomial_coeff_table B(n, dim_max + 2);

	std::unordered_map<index_t, value_t>* simplicialFiltration = new std::unordered_map<index_t, value_t>[dim_max + 1];

	while (std::getline(input_stream, line)) {
		std::cout << line << std::endl;
		string_end = line.find(delimiter);
		//  std::cout << "substring " << line.substr(1,string_end-1) << std::endl;
		stream.str("");
		stream.clear();
		stream << line.substr(1, string_end - 1);
		// std::cout << line[line.length()] << "gotcha ";

		// std::cout << stream.str() << std::endl;
		new_simplex.clear();
		dim = 0;
		while (stream >> num) {
			new_simplex.push_back(num);
			// std::cout << num << std::endl;
			dim++;
		}
		dim--;
		std::sort(new_simplex.rbegin(), new_simplex.rend());
		
		value_t filtration_value;
		(std::stringstream(line.substr(string_end + 1, line.length()))) >> filtration_value;
		index_t index = compute_index(new_simplex, B);
		// std::cout << "index " << index << std::endl;
		// std::cout << "dim " << dim << std::endl;
		simplicialFiltration[dim].insert(std::make_pair(index, filtration_value));
	}
	return simplicialFiltration;
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

	std::unordered_map<index_t, value_t>* simplicialFiltration;

	simplicialFiltration = read_file(filename ? file_stream : std::cin, n, dim_max);

	std::cout << "Complex of dimension " << dim_max << " with " << n << " points" << std::endl;
	dim_max = std::min(dim_max, n - 2);

	binomial_coeff_table binomial_coeff(n, dim_max + 2);
	std::vector<coefficient_t> multiplicative_inverse(multiplicative_inverse_vector(modulus));

	std::vector<diameter_index_t> columns_to_reduce;

	{
		union_find dset(n);
		std::vector<diameter_index_t> edges;
		for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
			if (simplicialFiltration[1].find(index) != simplicialFiltration[1].end()) {
				hash_map<index_t, value_t>::const_iterator S = simplicialFiltration[1].find(index);
				value_t diameter = S->second;
				if (diameter <= threshold) edges.push_back(std::make_pair(diameter, index));
			}
		}
		std::sort(edges.rbegin(), edges.rend(), greater_diameter_or_smaller_index<diameter_index_t>());

#ifdef PRINT_PERSISTENCE_PAIRS
		std::cout << "persistence intervals in dim 0:" << std::endl;
#endif

		std::vector<index_t> vertices_of_edge(2);
		for (auto e : edges) {
			vertices_of_edge.clear();
			get_simplex_vertices(get_index(e), 1, n, binomial_coeff, std::back_inserter(vertices_of_edge));
			index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

			if (u != v) {
#ifdef PRINT_PERSISTENCE_PAIRS
				auto v0 = simplicialFiltration[0].find(vertices_of_edge[0]);
				value_t f_val0 = v0->second;
				auto v1 = simplicialFiltration[0].find(vertices_of_edge[1]);
				value_t f_val1 = v1->second;
				value_t paired_fval = std::max(f_val0, f_val1);
				if (get_diameter(e) > 0)
					std::cout << " [" << paired_fval << "," << get_diameter(e) << ")" << std::endl;
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
				auto unpaired_vertex = simplicialFiltration[0].find(i);
				value_t f_val = unpaired_vertex->second;
				std::cout << " [" << f_val << ","
				          << ")" << std::endl;
			}
#endif
	}

	for (index_t dim = 1; dim < dim_max; ++dim) {

		hash_map<index_t, index_t> pivot_column_index;
		pivot_column_index.reserve(columns_to_reduce.size());

		compute_pairs(columns_to_reduce, pivot_column_index, dim, n, threshold, modulus, multiplicative_inverse,
		               simplicialFiltration, binomial_coeff);

		if (dim < dim_max - 1) {
			assemble_columns_to_reduce(columns_to_reduce, pivot_column_index, simplicialFiltration, dim, n, threshold,
			                            binomial_coeff);
		}
	}

	return 0;
}
