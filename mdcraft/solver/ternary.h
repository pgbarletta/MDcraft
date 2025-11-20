#pragma once

#include <array>
#include <functional>

#include <mdcraft/solver/isolver.h>

#include <mdcraft/tools/threads.h>
#include <mdcraft/solver/thermostat/base.h>
#include <mdcraft/solver/boundary/base.h>

#include <mdcraft/solver/potential/base.h>

namespace mdcraft::solver {

using ::mdcraft::data::vector;
using ::mdcraft::data::matrix;

using ::mdcraft::tools::Threads;
using ::mdcraft::data::Atoms;

using NeibsList        = ::mdcraft::neibs::List;
using NeibsListOne     = ::mdcraft::neibs::ListOne;
using Thermostat       = thermostat::Base;
using BaseTernary      = potential::BaseTernary;


class Ternary : public ISolver {
public:
	/** \brief
		\~russian Конструктор класса.
		\~english Class constructor.
		\~
		\param[in] threads
			\~russian пул тредов для параллельного счета.
			\~english thread pool for parallel calculation.
			\~
	*/
	Ternary(
		BaseTernary& potential,
		Thermostat&  thermostat = thermostat::dummy_thermostat,
		Threads&     threads    = tools::dummy_pool
	);
	/** \brief
		\~russian Деструктор класса.
		\~english Class destructor.
		\~
	*/
	virtual ~Ternary() {};

	/** \brief
		\~russian Подготовить  к расчету.
		\~english Prepare atoms to solve.
		\~
		\param[in] atoms
			\~russian данные атомов.
			\~english Atoms data.
			\~
		\details
		\~russian Функция запускает другую
			функцию, \p prepare_one, для всех элементов 
			входного атомы. 
		\~english This function executes \p prepare_one
			for all the \p atoms.
		\~
	*/
	void prepare(
		Atoms&     atoms, 
		Atoms&     neibs,
		NeibsList& nlist
	) override;

	/** \brief
		\~russian Провести расчет сил.
		\~english Calculate forces actin on atoms.
		\~
		\param[in] atoms
			\~russian Данные атомов.
			\~english Atoms data.
			\~
		\details
		\~russian Функция запускает другую
			функцию, \p solve_one, для всех элементов 
			входного атомы. 
		\~english This function executes \p solve_one
			for all the \p atoms.
		\~
	*/
	void forces(
		Atoms&     atoms,
		Atoms&     neibs,
		const NeibsList& nlist
	);

	/** \brief
		\~russian Провести расчет вириалов.
		\~english Solve atoms with virial calculation.
		\~
		\param[in] atoms
			\~russian Данные атомов.
			\~english Atoms data.
			\~
		\details
		\~russian Функция запускает другую
			функцию, \p solve_one, для всех элементов 
			входного атомы. 
		\~english This function executes \p solve_one
			for all the \p atoms.
		\~
	*/
	void virials(
		Atoms&     atoms, 
		Atoms&     neibs, 
		const NeibsList& nlist
	);

protected:

	BaseTernary& potential; /**<
			\~russian Межатомный трехчастичный потенциал.
			\~english Interatomic ternary potential.
			\~
		*/

	// very sad that I cannot have a const version
	// because apply_one members changes Thermostat class state
	Thermostat& thermostat;

	Threads& pool;
		/**<\~russian ссылка на пул тредов для работы функции 
				в многопоточном режиме.
			\~english reference to a pool of threads for 
				multithreaded execution.
			\~
		*/
};

} // mdcraft