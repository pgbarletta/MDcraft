#include <cmath>
#include <iostream>
#include <mdcraft/neibs/grid.h>

namespace mdcraft::neibs {

constexpr double max_double = std::numeric_limits<double>::max();
constexpr double min_double = std::numeric_limits<double>::lowest();


Grid::Grid(Domain& domain, Threads& threads)
  : domain(domain), pool(threads)
{
	initialize_full();
}

void Grid::initialize_full() {
	// Одна ячейка сетки на весь domain
	// One cell for the whole domain
	cell_sizes = {domain.xsize(), domain.ysize(), domain.zsize()};
	domain_num_cells = {1, 1, 1};
	grid_offset = {0, 0, 0};
	min_coords = {domain.xmin(), domain.ymin(), domain.zmin()};
	max_coords = {domain.xmax(), domain.ymax(), domain.zmax()};
	num_cells = {1, 1, 1};
}

void Grid::initialize_empty() {
	// Одна ячейка сетки на весь domain, но локальных ячеек нет
	// One grid cell for the whole domain, but no local cells
	cell_sizes = {domain.xsize(), domain.ysize(), domain.zsize()};
	domain_num_cells = {1, 1, 1};

	// Локальные ячейки "за границей"
	// Local cells "outside the border"
	grid_offset = {1, 1, 1};
	min_coords = {NAN, NAN, NAN};
	max_coords = {NAN, NAN, NAN};
	num_cells = {0, 0, 0};
}


void Grid::build(const Atoms& atoms) {
	if (atoms.empty()) {
		initialize_empty();
		return;
	}

	// bmin, bmax -- границы bounding box
	// cell_size  -- максимальный радиус поиска
	// bmin, bmax -- bounding box borders
	// cell_size  -- max search radius
	auto [bmin, bmax, cell_size] = find_bounds(atoms);

	// Установить число ячеек сетки на весь domain
	// Set number of grid cells for the whole domain
	domain_num_cells = {1, 1, 1};
	for (int i = 0; i < domain.dim(); ++i) {
		domain_num_cells[i] = std::max(1, static_cast<cell_id_t>(std::floor(domain.size(i) / cell_size)));
	}

	// Установить размеры ячеек сетки
	// Set grid cells' sizes
	cell_sizes = {domain.xsize(), domain.ysize(), domain.zsize()};
	for (int i = 0; i < domain.dim(); ++i) {
		cell_sizes[i] /= domain_num_cells[i];
	}

	// Смещение локальной сетки внутри domain
	// Local grid offset inside domain
	grid_offset = global_index3D(bmin);

	// Верхняя граница индексов для локальной сетки
	// Upper bound for indices of local grid
	auto upper_bounds = global_index3D(bmax);
	for (int i = 0; i < domain.dim(); ++i) {
		upper_bounds[i] = std::min(upper_bounds[i] + 1, domain_num_cells[i]);
	}

	// Число ячеек локальной сетки
	// Number of local grid cells
	num_cells = {1, 1, 1};
	for (int i = 0; i < domain.dim(); ++i) {
		num_cells[i] = upper_bounds[i] - grid_offset[i];
	}

	// Точные границы локальной сетки (зачем они нужны?)
	// Precise boundaries of local grid (what for?)
	min_coords = {domain.xmin(), domain.ymin(), domain.zmin()};
	max_coords = {domain.xmax(), domain.ymax(), domain.zmax()};
	for (int i = 0; i < domain.dim(); ++i) {
		min_coords[i] = domain.bmin(i) + cell_sizes[i] * grid_offset[i];
		max_coords[i] = domain.bmin(i) + cell_sizes[i] * upper_bounds[i];
	}
}

void Grid::set_indices(Atoms& atoms, std::vector<element_ids>& indices) const {
	indices.resize(atoms.size());
	auto find_one = [&](point_id i) {
		indices[i].cell_id = cell_index(local_index3D(atoms[i].r));
		indices[i].atom_id = i;
	};

	pool.parallel_for(0ul, atoms.size(), find_one);
}

Index3D Grid::global_index3D(const vector& v) const {
	Index3D c = {0, 0, 0};
	for (int i = 0; i < domain.dim(); ++i)
		c[i] = std::floor((v[i] - domain.bmin(i)) / cell_sizes[i]);
	return c;
}

Index3D Grid::local_index3D(const vector& v) const {
	Index3D c = global_index3D(v);
	for (int i = 0; i < domain.dim(); ++i)
		c[i] -= grid_offset[i];
	return c;
}

Index3D Grid::local_index3D(cell_id_t i) const {
	int n_x = num_cells[0];
	int n_y = num_cells[1];
	int n_xy = n_x * n_y;

	Index3D c;
	c[2] = domain.dim() > 2 ? i / n_xy : 0;
	c[1] = domain.dim() > 1 ? (i - c[2] * n_xy) / n_x : 0;
	c[0] = domain.dim() > 0 ? i - (n_x * (c[1] + n_y * c[2])) : 0;

	return c;
}

cell_id_t Grid::cell_index(const Index3D& ijk) const {
	return ijk[0] + num_cells[0] * (ijk[1] + num_cells[1] * ijk[2]);
}

std::tuple<std::array<cell_id_t, 27>, Flag3D>
	Grid::adjacent_cells(const vector& r) const {

	Index3D cell = global_index3D(r);

	auto n_x = domain_num_cells[0];
	auto n_y = domain_num_cells[1];
	auto n_z = domain_num_cells[2];

	int list_size = 0;
	std::array<cell_id_t, 27> list;

	for (cell_id_t x: {-1, 0, 1}) {
		if (domain.xperiodic()) {
			if (n_x == 1 && x != 0) continue;
			if (n_x == 2 && x == 1) continue;
		}

		for (cell_id_t y: {-1, 0, 1}) {
			if (domain.yperiodic()) {
				if (n_y == 1 && y != 0) continue;
				if (n_y == 2 && y == 1) continue;
			}

			for (cell_id_t z: {-1, 0, 1}) {
				if (domain.zperiodic()) {
					if (n_z == 1 && z != 0) continue;
					if (n_z == 2 && z == 1) continue;
				}

				// Индексы соседней ячейки на глобальной сетке
				// Neighvour cell indices on global grid
				Index3D neib = { cell[0] + x, cell[1] + y, cell[2] + z };

				// Учтем периодичность в глобальной индексации
				// Take into account periodicity in global indexing
				if (domain.xperiodic()) { neib[0] = (neib[0] + n_x) % n_x; }
				if (domain.yperiodic()) { neib[1] = (neib[1] + n_y) % n_y; }
				if (domain.zperiodic()) { neib[2] = (neib[2] + n_z) % n_z; }

				// Переведем в локальный индекс
				// Transform into local index
				neib[0] -= grid_offset[0];
				neib[1] -= grid_offset[1];
				neib[2] -= grid_offset[2];

				// Индекс внутри локальной сетки? Можем добавить ячейку.
				// Index inside local grid? Can add a cell.
				if (0 <= neib[0] && neib[0] < num_cells[0] &&
					0 <= neib[1] && neib[1] < num_cells[1] &&
					0 <= neib[2] && neib[2] < num_cells[2]) {

					list[list_size++] = cell_index(neib);
				}
			}
		}
	}

	for (int i = list_size; i < 27; ++i) {
		list[i] = empty_cell;
	}

	// Вроде как этот флаг сейчас не нужен
	// Now this flag doesn't seem to be required
	Flag3D is_edge = {false, false, false};
	for (int i = 0; i < domain.dim(); ++i) {
		is_edge[i] = cell[i] == 0 || cell[i] == num_cells[i] - 1;
	}
	return {list, is_edge};
}

std::tuple<vector, vector, double> Grid::find_bounds(const Atoms& atoms) const {
	// Вспомогательная структура
	// Auxillary structure
	struct Info {
		double max_rns   = 0.0;
		vector min_coord = {max_double, max_double, max_double};
		vector max_coord = {min_double, min_double, min_double};

		// Операция свёртки
		// Convolution operation
		Info &operator&=(const Info &foo) {
			max_rns   = std::max(max_rns, foo.max_rns);
			min_coord = min_coord.cwiseMin(foo.min_coord);
			max_coord = max_coord.cwiseMax(foo.max_coord);
			return *this;
		}
	};

	// Свёртка
	// Convoluion
	auto res = pool.reduce(
			atoms.begin(), atoms.end(), Info{},
			[](const Atom& atom) -> Info {
				return {atom.rns, atom.r, atom.r};
			});

	vector min_coords = vector::Zero();
	vector max_coords = vector::Zero();
	for (int i = 0; i < domain.dim(); ++i) {
		min_coords[i] = res.min_coord[i];
		max_coords[i] = res.max_coord[i];
	}

	return {min_coords, max_coords, res.max_rns};
}

void Grid::print() const {
	std::cout << "Domain: ["
			  << domain.xmin() << ", " << domain.xmax() << "] x ["
	          << domain.ymin() << ", " << domain.ymax() << "] x ["
	          << domain.zmin() << ", " << domain.zmax() << "];\n";
	std::cout << "  size:       " << domain.xsize() << " x " << domain.ysize() << " x " << domain.zsize() << "\n";
	std::cout << "  num cells:  " << domain_num_cells[0] << " x " << domain_num_cells[1] << " x " << domain_num_cells[2] << "\n";
	std::cout << "  cell_sizes: " << cell_sizes.transpose() << "\n";

	std::cout << "Subdomain: ["
			  << min_coords.x() << ", " << max_coords.x() << "] x ["
			  << min_coords.y() << ", " << max_coords.y() << "] x ["
			  << min_coords.z() << ", " << max_coords.z() << "];\n";
	std::cout << "  size:     " << (max_coords.x() - min_coords.x()) << " x "
	                            << (max_coords.y() - min_coords.y()) << " x "
	                            << (max_coords.z() - min_coords.z()) << "\n";
	std::cout << "  num cells:  " << num_cells[0] << " x " << num_cells[1] << " x " << num_cells[2] << "\n";
	std::cout << "  offsets:  " << grid_offset[0] << ", " << grid_offset[1] << ", " << grid_offset[2] << "\n";
}

} // namespace mdcraft::neibs