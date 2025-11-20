#pragma once
#include <functional>
#include <mdcraft/data/atom.h>
#include <mdcraft/tools/threads.h>
#include <mdcraft/lattice/domain.h>

namespace mdcraft::neibs {

using data::vector;
using data::Atom;
using data::Atoms;

using tools::Threads;
using lattice::Domain;

// Type for indexing elements (atoms)
using point_id = std::size_t;

// Type for the grid cells indexing
using cell_id_t = int;
using Index3D = std::array<cell_id_t, 3>;

static constexpr cell_id_t empty_cell = -1;

using Flag3D = std::array<bool, 3>;

struct element_ids {
	cell_id_t cell_id = 0;  ///< Owner cell index
	point_id  atom_id = 0;  ///< Atom index

	bool operator<(const element_ids& rhs) const {
		return cell_id < rhs.cell_id;
	}

	bool operator==(const element_ids& rhs) const {
		return cell_id == rhs.cell_id;
	}
};

/// @brief
///		@russian Виртуальная декартова сетка поверх расчетной области. Внутри
///		глобальной сетки выделена подсетка меньшего размера, которая содержит
///		только актуальную часть элементов.
///		Сетку для целой расчетной области будем называть глобальной сеткой, а
///		внутреннюю подсетку -- локальной. На обеих виртуальных сетках введена
///		трёхмерная индексация ячеек (глобальная и локальная).
///		Исключительно для локальной сетки также введена одномерная индексация.
///		@english Virual cartesian grid over the computational domain. Inside 
///		the global grid, a sub-grid of smaller size is selected. It contains 
///     only the relevant part of elements.
///     Let us call the whole domain grid a global grid, 
///     and the sub-grid a local grid. Over both grids a 3-dimentional cell 
///     indexing is introduced.

class Grid {
public:
	/// @brief
	///     @russian Конструктор класса
	///     @english Class constructor
	/// @param domain
	///		@russian Вычислительная область
	///		@english Computational domain
	/// @param threads
	///		@russian Реализация многопоточности
	///		@english Multithreading implementation
	explicit Grid(Domain& domain, Threads& threads = tools::dummy_pool);

	/// @brief
	///		@russian Построить сетку для множества элементов
	///		@english Build a grid for elements
	/// @param atoms
	///		@russian Массив элементов
	///		@english Array of elements
	void build(const Atoms& atoms);

	/// @brief
	///		@russian Получить 3D индекс ячейки для точки
	///		@english Get a 3D index of cell for a point
	Index3D local_index3D(const vector& v) const;

	/// @brief
	///		@russian Получить одномерный индекс ячейки локальной сетки
	///		@english Get plain 1D cell index for local grid
	/// @param cell
	///		@russian
	///		@english 3D cell index
	cell_id_t cell_index(const Index3D& cell) const;

	/// @brief
	///		@russian Установить одномерные локальные индексы элементов
	///		@english Set 1D indices of elements in a local grid
	/// @param atoms[in]
	///		@russian Массив элементов
	///		@english Array of elements
	///	@param indices[out]
	///		@russian Индексы элементов и ячеек сетки
	///		@english Indices of atoms and grid cells
	void set_indices(Atoms& atoms, std::vector<element_ids>& indices) const;

	/// @brief
	///		@russian Найти одномерные индексы соседних ячеек
	///		@english Get 1D indices of adjacent cells
	/// @param v
	///		@russian Координаты точки
	///		@english Point coordinates
	///	@return
	///		@russian Результирующий список ячеек. Индикатор того, находится ли
	///		ячейка сетки на границе сетки, по трем осям.
	///		@english Resulting indices list. Indicator of being boundary cell, by axes.
	std::tuple<std::array<cell_id_t, 27>, Flag3D>
		adjacent_cells(const vector& v) const;

	/// @brief
	///		@russian Максимальный линейный размер ячейки сетки
	///		@english Maximum linear size of the grid cell
	double cell_size() const { return cell_sizes.maxCoeff(); }

	/// @brief
	///		@russian Число ячеек локальной сетки
	///		@english Number of cells of the local grid
	cell_id_t cells_number() const {
		return num_cells[0] * num_cells[1] * num_cells[2];
	}

	/// @brief
	///		@russian Вывести информацию о сетке
	///		@english Print information about grid and subgrid
	void print() const;

private:
	/// @brief Unique cell in whole domain, local grid == global grid.
	void initialize_full();

	/// @brief Unique cell in whole domain, local grid is empty.
	void initialize_empty();

	/// @brief
	///		@russian Получить глобальный трехмерный индекс для точки
	///		@english Get a 3D index for a point
	Index3D global_index3D(const vector& v) const;

	/// @brief
	///		@russian Получить локальный 3D индекс ячейки по локальному 1D индексу
	///		@english Get a local 3D index of cell from a local 1D index.
	Index3D local_index3D(cell_id_t i) const;

	/// @return
	///		@russian Возвращает границы bounding box и максимальный радиус
	///		поиска соседей {bb_min, bb_max, max_search_rad}
	///		@english Return bounding box and maximum search radius
	std::tuple<vector, vector, double> find_bounds(const Atoms& atoms) const;


	/// @russian Ссылка на полную вычислительную область
	/// @english Reference to the full computational domain
	Domain& domain;

	/// @russian Число ячеек глобальной сетки по осям координат внутри domain
	/// @english Number of cells for the domain along Ox, Oy, Oz axes
	Index3D domain_num_cells{};

	/// @russian Размеры ячеек сетки по осям координат
	/// @english Grid cell sizes along Ox, Oy, Oz axes
	vector cell_sizes;

	/// @russian Смещение локальной сетки внутри сетки domain
	/// @english Offset of the local grid inside domain
	Index3D grid_offset{};

	/// @russian Число ячеек локальной сетки по трем осям
	/// @english Number of cells of the local grid by each axis
	Index3D num_cells{};

	/// @russian Минимум локальной сетки по Ox, Oy, Oz
	/// @english Local grid minimum by Ox, Oy, Oz
	vector min_coords;

	/// @russian Максимум локальной сетки по Ox, Oy, Oz
	/// @english Local grid maximum by Ox, Oy, Oz
	vector max_coords;

	/// @russian Ссылка на объект тредов
	/// @english Thread pool object reference
	Threads& pool;
};

} // namespace mdcraft::neibs