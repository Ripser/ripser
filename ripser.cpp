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

#define PRINT_PERSISTENCE_PAIRS

#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>

template <class Key, class T> class hash_map : public std::unordered_map<Key, T> {};

typedef float value_t;
typedef int64_t index_t;

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


struct diameter_index_t {
	value_t diameter, diameter_sub;
	index_t index;
	
	diameter_index_t(): diameter(0), diameter_sub(0), index(-1) {}
	diameter_index_t(const diameter_index_t& i): diameter(i.diameter_sub), diameter_sub(i.diameter_sub), index(i.index) {}
	diameter_index_t(const value_t _diameter, const value_t _diameter_sub, const index_t _index): diameter(_diameter), diameter_sub(_diameter_sub), index(_index) {}
};

value_t get_diameter(const diameter_index_t& i) { return i.diameter; }
value_t get_diameter_sub(const diameter_index_t& i) { return i.diameter_sub; }
index_t get_index(const diameter_index_t& i) { return i.index; }
	
template <typename Entry> struct greater_diameter_or_smaller_index {
	bool use_diameter_sub;
	greater_diameter_or_smaller_index(bool _use_diameter_sub) : use_diameter_sub(_use_diameter_sub) {}
	bool operator()(const Entry& a, const Entry& b) {
		return use_diameter_sub ?
		(get_diameter_sub(a) > get_diameter_sub(b)) ||
		       ((get_diameter_sub(a) == get_diameter_sub(b)) && (get_index(a) < get_index(b))):
		(get_diameter(a) > get_diameter(b)) ||
		       ((get_diameter(a) == get_diameter(b)) && (get_index(a) < get_index(b)));
	}
};

class compressed_lower_distance_matrix {
public:
	std::vector<value_t> distances;
	std::vector<value_t*> rows;

	void init_rows() {
		value_t* pointer = &distances[0];
		for (index_t i = 1; i < size(); ++i) {
			rows[i] = pointer;
			pointer += i;
		}
	}

	compressed_lower_distance_matrix(std::vector<value_t>&& _distances)
	    : distances(std::move(_distances)), rows((1 + std::sqrt(1 + 8 * distances.size())) / 2) {
		assert(distances.size() == size() * (size() - 1) / 2);
		init_rows();
	}

	template <typename DistanceMatrix>
	compressed_lower_distance_matrix(const DistanceMatrix& mat)
	    : distances(mat.size() * (mat.size() - 1) / 2), rows(mat.size()) {
		init_rows();

		for (index_t i = 1; i < size(); ++i)
			for (index_t j = 0; j < i; ++j) rows[i][j] = mat(i, j);
	}

	value_t operator()(const index_t i, const index_t j) const {
		return i == j ? 0 : i < j ? rows[j][i] : rows[i][j];
	}

	size_t size() const { return rows.size(); }
};

class euclidean_distance_matrix {
public:
	std::vector<std::vector<value_t>> points;

	euclidean_distance_matrix(std::vector<std::vector<value_t>>&& _points)
	    : points(std::move(_points)) {
			for (auto p: points) {
				assert(p.size() == points.front().size());
			}
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

template <typename Heap> diameter_index_t pop_pivot(Heap& column) {
	if (column.empty())
		return diameter_index_t();
	else {
		auto pivot = column.top();

		column.pop();
		while (!column.empty() && get_index(column.top()) == get_index(pivot)) {
			column.pop();
			if (column.empty())
				return diameter_index_t();
			else {
				pivot = column.top();
				column.pop();
			}
		}
		return pivot;
	}
}

template <typename Heap> diameter_index_t get_pivot(Heap& column) {
	diameter_index_t result = pop_pivot(column);
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

class ripser {
	compressed_lower_distance_matrix dist, dist_sub;
	index_t dim_max, n;
	value_t threshold;
	const binomial_coeff_table binomial_coeff;
	mutable std::vector<index_t> vertices;
	mutable std::vector<diameter_index_t> coface_entries;

public:
	ripser(compressed_lower_distance_matrix&& _dist, compressed_lower_distance_matrix&& _dist_sub, index_t _dim_max, value_t _threshold)
    : dist(std::move(_dist)), dist_sub(std::move(_dist_sub)), dim_max(std::min(_dim_max, index_t(dist.size() - 2))),
	      n(dist.size()), threshold(_threshold), binomial_coeff(n, dim_max + 2) {}

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
    
    value_t compute_diameter_sub(const index_t index, index_t dim) const {
        value_t diam = -std::numeric_limits<value_t>::infinity();
        
        vertices.clear();
        get_simplex_vertices(index, dim, dist_sub.size(), std::back_inserter(vertices));
        
        for (index_t i = 0; i <= dim; ++i)
            for (index_t j = 0; j < i; ++j) {
                diam = std::max(diam, dist_sub(vertices[i], vertices[j]));
            }
        return diam;
    }

	class simplex_coboundary_enumerator {
	private:
		index_t idx_below, idx_above, v, k;
		std::vector<index_t> vertices;
		const diameter_index_t simplex;
		const compressed_lower_distance_matrix& dist;
		const compressed_lower_distance_matrix& dist_sub;
		const binomial_coeff_table& binomial_coeff;

	public:
		simplex_coboundary_enumerator(const diameter_index_t _simplex, index_t _dim,
		                              const ripser& parent)
		    : idx_below(get_index(_simplex)), idx_above(0), v(parent.n - 1), k(_dim + 1),
              vertices(_dim + 1), simplex(_simplex), dist(parent.dist), dist_sub(parent.dist_sub),
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

		diameter_index_t next() {
			value_t coface_diameter = get_diameter(simplex);
			value_t coface_diameter_sub = get_diameter_sub(simplex);
			for (index_t w : vertices) {
				coface_diameter = std::max(coface_diameter, dist(v, w));
				coface_diameter_sub = std::max(coface_diameter_sub, dist_sub(v, w));
			}
			index_t coface_index = idx_above + binomial_coeff(v--, k + 1) + idx_below;
			return diameter_index_t(coface_diameter, coface_diameter_sub, coface_index);
		}
	};

	void compute_barcodes();

	void assemble_columns_to_reduce(std::vector<diameter_index_t>& columns_to_reduce,
	                                hash_map<index_t, index_t>& pivot_column_index, index_t dim) {
		index_t num_simplices = binomial_coeff(n, dim + 1);

		columns_to_reduce.clear();

		for (index_t index = 0; index < num_simplices; ++index) {
			if (pivot_column_index.find(index) == pivot_column_index.end()) {
				value_t diameter = compute_diameter(index, dim);
				value_t diameter_sub = compute_diameter_sub(index, dim);
                if (diameter_sub <= threshold) {
					columns_to_reduce.push_back(diameter_index_t(diameter, diameter_sub, index));
                }
			}
		}

		std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
		          greater_diameter_or_smaller_index<diameter_index_t>(true));
	}
	
	template <typename Column, typename Iterator>
	diameter_index_t add_coboundary_and_get_pivot(Iterator column_begin,
												  Iterator column_end,
												  Column& working_reduction_column,
												  Column& working_coboundary,
												  const index_t& dim,
												  hash_map<index_t, index_t>& pivot_column_index,
												  bool& might_be_apparent_pair,
												  bool image)
	{
		for (auto it = column_begin;
		it != column_end; ++it) {
			diameter_index_t simplex = *it;
			
			working_reduction_column.push(simplex);
			
			coface_entries.clear();
			simplex_coboundary_enumerator cofaces(simplex, dim, *this);
			while (cofaces.has_next()) {
				diameter_index_t coface = cofaces.next();
				value_t diam_coface = image ? get_diameter(coface) : get_diameter_sub(coface);
				
				if (diam_coface <= threshold) {
					coface_entries.push_back(coface);
					if (might_be_apparent_pair &&
						(get_diameter_sub(simplex) == diam_coface)) {
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
		
		return get_pivot(working_coboundary);
	}
	

	void compute_pairs(std::vector<diameter_index_t>& columns_to_reduce,
	                   hash_map<index_t, index_t>& pivot_column_index, index_t dim, bool image) {

		compressed_sparse_matrix<diameter_index_t> reduction_matrix;

		for (index_t index_column_to_reduce = 0; index_column_to_reduce < columns_to_reduce.size(); ++index_column_to_reduce) {
			auto column_to_reduce = columns_to_reduce[index_column_to_reduce];

			std::priority_queue<diameter_index_t, std::vector<diameter_index_t>,
			                    greater_diameter_or_smaller_index<diameter_index_t>>
			    working_reduction_column(!image), working_coboundary(!image);
            
			value_t diameter = get_diameter_sub(column_to_reduce);
			
			index_t index_column_to_add = index_column_to_reduce;
			
			diameter_index_t pivot;
            
			// initialize reduction_matrix as identity matrix
			reduction_matrix.append_column();
			reduction_matrix.push_back(diameter_index_t(column_to_reduce));
			
			bool might_be_apparent_pair = (index_column_to_reduce == index_column_to_add);

			while (true) {
                
                pivot = add_coboundary_and_get_pivot(reduction_matrix.cbegin(index_column_to_add),
                                                         reduction_matrix.cend(index_column_to_add),
                                                         working_reduction_column,
                                                         working_coboundary,
                                                         dim,
                                                         pivot_column_index,
                                                         might_be_apparent_pair,
                                                         image);
                
				if (get_index(pivot) != -1) {
					auto pair = pivot_column_index.find(get_index(pivot));

					if (pair != pivot_column_index.end()) {
						index_column_to_add = pair->second;
						continue;
					}
				} else {
#ifdef PRINT_PERSISTENCE_PAIRS
                    if (!image) {
                        std::cout << " [" << diameter << ", )" << std::endl << std::flush;
                    }
#endif
                    break;
                }
                
#ifdef PRINT_PERSISTENCE_PAIRS
                if (image) {
                    value_t death = get_diameter(pivot);
                    if (diameter != death) {
                        std::cout << " [" << diameter << "," << death << ")" << std::endl << std::flush;
                    }
                }
#endif

                pivot_column_index.insert(std::make_pair(get_index(pivot), index_column_to_reduce));
                
                
// replace current column of reduction_matrix (with a single diagonal 1 entry)
// by reduction_column (possibly with a different entry on the diagonal)

				reduction_matrix.pop_back();
				
				while (true) {
					diameter_index_t e = pop_pivot(working_reduction_column);
					if (get_index(e) == -1) break;
					reduction_matrix.push_back(e);
				}
				break;
			}
		}
	}
};

enum file_format {
	LOWER_DISTANCE_MATRIX,
	DISTANCE_MATRIX,
	POINT_CLOUD,
	RIPSER
};

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

	std::cout << "point cloud with " << n << " points in dimension "
	          << eucl_dist.points.front().size() << std::endl;

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

compressed_lower_distance_matrix read_ripser(std::istream& input_stream) {
	std::vector<value_t> distances;
	while (!input_stream.eof()) distances.push_back(read<value_t>(input_stream));
	return compressed_lower_distance_matrix(std::move(distances));
}

compressed_lower_distance_matrix read_file(std::istream& input_stream, file_format format) {
	switch (format) {
	case LOWER_DISTANCE_MATRIX:
		return read_lower_distance_matrix(input_stream);
	case DISTANCE_MATRIX:
		return read_distance_matrix(input_stream);
	case POINT_CLOUD:
		return read_point_cloud(input_stream);
	case RIPSER:
		return read_ripser(input_stream);
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
	    << "                     distance       (full distance matrix)" << std::endl
	    << "                     point-cloud    (point cloud in Euclidean space)" << std::endl
	    << "                     binary         (distance matrix in Ripser binary file format)"
	    << std::endl
        << "  --subfiltration <f>  use f as second filtration for image persistence" << std::endl
	    << "  --dim <k>        compute persistent homology up to dimension <k>" << std::endl
	    << "  --threshold <t>  compute Rips complexes up to diameter <t>" << std::endl
	    << std::endl;

	exit(exit_code);
}

int main(int argc, char** argv) {

	const char* filename = nullptr;
    const char* filename_sub = nullptr;

	file_format format = DISTANCE_MATRIX;

	index_t dim_max = 1;
	value_t threshold = std::numeric_limits<value_t>::infinity();

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
		} else if (arg == "--format") {
			std::string parameter = std::string(argv[++i]);
			if (parameter == "lower-distance")
				format = LOWER_DISTANCE_MATRIX;
			else if (parameter == "distance")
				format = DISTANCE_MATRIX;
			else if (parameter == "point-cloud")
				format = POINT_CLOUD;
			else if (parameter == "ripser")
				format = RIPSER;
			else
				print_usage_and_exit(-1);
        } else if (arg == "--subfiltration") {
            if (filename_sub) { print_usage_and_exit(-1); }
            filename_sub = argv[++i];
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

	compressed_lower_distance_matrix dist = read_file(filename ? file_stream : std::cin, format);

	std::cout << "distance matrix with " << dist.size() << " points" << std::endl;

	value_t min = std::numeric_limits<value_t>::infinity(), max = -std::numeric_limits<value_t>::infinity();
	
	for (auto d: dist.distances) {
		if (d != std::numeric_limits<value_t>::infinity() ) {
			min = std::min(min, d);
			max = std::max(max, d);
		} else {
			threshold = std::min(threshold, std::numeric_limits<value_t>::max());
		}
	}
	
	std::cout << "value range: [" << min << "," << max << "]"
	          << std::endl;

    file_stream.close();
    
    if (!filename_sub) filename_sub = filename;
    
    std::ifstream file_stream_sub(filename_sub);
    if (file_stream_sub.fail()) {
        std::cerr << "couldn't open file " << filename_sub << std::endl;
        exit(-1);
    }
    
    compressed_lower_distance_matrix dist_sub = read_file(file_stream_sub, format);
    
    std::cout << "subfiltration distance matrix with " << dist_sub.size() << " points" << std::endl;
    
    min = std::numeric_limits<value_t>::infinity();
    max = -std::numeric_limits<value_t>::infinity();
    
    for (auto d: dist_sub.distances) {
        if (d != std::numeric_limits<value_t>::infinity() ) {
            min = std::min(min, d);
            max = std::max(max, d);
        } else {
            threshold = std::min(threshold, std::numeric_limits<value_t>::max());
        }
    }
    
    std::cout << "subfiltration value range: [" << min << "," << max << "]"
    << std::endl;
    
	ripser(std::move(dist), std::move(dist_sub), dim_max, threshold).compute_barcodes();
}

void ripser::compute_barcodes() {

	std::vector<diameter_index_t> columns_to_reduce;
	std::vector<index_t> vertices_of_edge(2);


	{
		union_find dset(n), dset_sub(n);
		std::vector<diameter_index_t> edges;
		for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
			value_t diameter = compute_diameter(index, 1);
			value_t diameter_sub = compute_diameter_sub(index, 1);
			if (diameter_sub <= threshold) edges.push_back(diameter_index_t(diameter, diameter_sub, index));
		}


		std::sort(edges.rbegin(), edges.rend(),
				  greater_diameter_or_smaller_index<diameter_index_t>(false));

		//computing 0-dimensional barcode of image by computing that of upper filtration
#ifdef PRINT_PERSISTENCE_PAIRS
        std::cout << "persistence intervals in dim 0:" << std::endl;
#endif
        
        for (auto& e : edges) {
            vertices_of_edge.clear();
            get_simplex_vertices(get_index(e), 1, n, std::back_inserter(vertices_of_edge));
            index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);
            
            if (u != v) {
#ifdef PRINT_PERSISTENCE_PAIRS
                if (get_diameter(e) != 0)
                    std::cout << " [0," << get_diameter(e) << ")" << std::endl;
#endif
                dset.link(u, v);
            }
        }
        
		std::sort(edges.rbegin(), edges.rend(),
				  greater_diameter_or_smaller_index<diameter_index_t>(true));
		
		//computing 0-dimensional barcode of lower filtration to initialize matrix for higher dimensions
		for (auto& e : edges) {
			vertices_of_edge.clear();
			get_simplex_vertices(get_index(e), 1, n, std::back_inserter(vertices_of_edge));
			index_t u = dset_sub.find(vertices_of_edge[0]), v = dset_sub.find(vertices_of_edge[1]);

			if (u != v) {
				dset_sub.link(u, v);
			} else
				columns_to_reduce.push_back(e);
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

        
#ifdef PRINT_PERSISTENCE_PAIRS
        for (index_t i = 0; i < n; ++i)
            if (dset.find(i) == i) std::cout << " [0, )" << std::endl << std::flush;
#endif
	}

	for (index_t dim = 1; dim <= dim_max; ++dim) {
		hash_map<index_t, index_t> pivot_column_index;
        
#ifdef PRINT_PERSISTENCE_PAIRS
        std::cout << "persistence intervals in dim " << dim << ":" << std::endl;
#endif
        pivot_column_index.reserve(columns_to_reduce.size());

        //reducing matrix corresponding to image
        compute_pairs(columns_to_reduce, pivot_column_index, dim, true);

        //reducing matrix corresponding to lower filtration
        pivot_column_index.clear();
        compute_pairs(columns_to_reduce, pivot_column_index, dim, false);
        
		if (dim < dim_max) {
			assemble_columns_to_reduce(columns_to_reduce, pivot_column_index, dim + 1);
		}
	}
}
