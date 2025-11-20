#pragma once
#include <array>
#include <functional>
#include <mdcraft/configuration.h>
#include <mdcraft/tools/threads.h>
#include <mdcraft/data/atom.h>
#include <mdcraft/neibs/verlet-list.h>

#include <mdcraft/solver/pair.h>
#include <mdcraft/solver/potential/base.h>
#include <mdcraft/solver/thermostat/base.h>
#include <mdcraft/solver/boundary/base.h>

namespace mdcraft {
namespace solver {

using Potential = potential::Base;
using Boundary  = boundary::Base;

/** \brief
	\~russian Класс решателя для расчета сил межатомного 
		взаимодействия с одним потенциалом.
	\~english Single solver class for interatomic forces
	          evaluation.
	\~
*/
class Multi : public Pair {
public:
	/** \brief
		\~russian Конструктор класса.
		\~english Class constructor.
		\~
		\param[in] potential
			\~russian Межатомный потенциал.
			\~english Interatomic potential.
			\~
		\param[in] threads
			\~russian пул тредов для параллельного счета.
			\~english thread pool for parallel calculation.
			\~
	*/
	Multi(
		Boundary&   boundary   = boundary::dummy_boundary,
		Thermostat& thermostat = thermostat::dummy_thermostat,
		std::size_t nkinds     = 1,
		Threads&    threads    = tools::dummy_pool
	);
	
	void add_potential(std::size_t i, std::size_t j, Potential& pot);

	/** \brief
		\~russian Подготовить поля хранилища 
			для одной частицы к решению.
			(Например, записать 0 в поля структуры, 
			где будет суммироваться решение.)
		\~english Prepare storage fields 
			of a single particle to solve.
			Ex., write 0 to fields which a used for result 
			accomodating.
		\~
		\param[in] atom
			\~russian итератор, обеспечивающий доступ 
				к частице.
			\~english iterator which gives access to 
				a particle.
			\~
	*/
	void prepare_one(
		Atoms::iterator atom, 
		Atoms&          neibs, 
		NeibsListOne&   nlist
	) override;

	/** \brief
		\~russian Рассчитать зарядовую плотность  
			для одной частицы.
		\~english Calculate charge density 
			of a single particle.
		\~
		\param[in] atom
			\~russian итератор, обеспечивающий доступ 
				к частице.
			\~english iterator which gives access to 
				a particle.
			\~
		\param[in] neibs
			\~russian хранилище с частицами, в котором
				есть соседние частицы с \p atom.
			\~english a Storage with particles 
				which contains neighboring particles
				with \p atom.
			\~
		\param[in] list
			\~russian список индексов частиц из \p neibs,
				которые являются соседними с \p atom.
			\~english a list of indices of particles 
				from \p neibs, which are neighbors of \p atom.
			\~
	*/
	void charge_one(
		Atoms::iterator atom, 
		Atoms&          neibs, 
		NeibsListOne&   nlist
	) override;

	/** \brief
		\~russian Рассчитать силу для одной частицы.
		\~english Solve for a particle.
		\~
		\param[in] atom
			\~russian итератор, обеспечивающий доступ 
				к частице.
			\~english iterator which gives access to 
				a particle.
			\~
		\param[in] neibs
			\~russian хранилище с частицами, в котором
				есть соседние частицы с \p atom.
			\~english a Storage with particles 
				which contains neighboring particles
				with \p atom.
			\~
		\param[in] list
			\~russian список индексов частиц из \p neibs,
				которые являются соседними с \p atom.
			\~english a list of indices of particles 
				from \p neibs, which are neighbors of \p atom.
			\~
	*/
	void force_one(
		Atoms::iterator atom, 
		Atoms&          neibs, 
		NeibsListOne&   nlist
	) override;

	/** \brief
		\~russian Рассчитать силу, вириал и потенциальную энергию для одной частицы.
		\~english Solve for a particle with virial and potential energy calculation.
		\~
		\param[in] atom
			\~russian итератор, обеспечивающий доступ 
				к частице.
			\~english iterator which gives access to 
				a particle.
			\~
		\param[in] neibs
			\~russian хранилище с частицами, в котором
				есть соседние частицы с \p me.
			\~english a Storage with particles 
				which contains neighboring particles
				with \p me.
			\~
		\param[in] list
			\~russian список индексов частиц из \p neibs,
				которые являются соседними с \p me.
			\~english a list of indices of particles 
				from \p neibs, which are neighbors of \p me.
			\~
	*/
	void virial_one(
		Atoms::iterator atom, 
		Atoms&          neibs, 
		NeibsListOne&   nlist
	) override;

private:
	
	std::vector<std::vector<Potential*>> potentials; /**<
		// \~russian Межатомные потенциалы.
		// \~english Interatomic potentials.
		// \~
	*/
	std::size_t nkinds;
};

} // solver
} // mdcraft
