#pragma once

#include <vector>
#include <functional>
#include <future>

#include <mdcraft/configuration.h>
#include <mdcraft/data/atom.h>

#include <mdcraft/tools/threads.h>
#include <mdcraft/neibs/verlet-list.h>

#include <mdcraft/solver/isolver.h>

namespace mdcraft::solver::stepper {

using tools::Threads;
using data::Atoms;
using NeibsList = neibs::List;

/** \brief
	\~russian Базовый класс интегрирования уравнений движения.
	\~english Base class for equations of motion integration.
	\~
*/
class Base {
public:
	/** \brief
		\~russian Конструктор класса.
		\~english Class constructor.
		\~
		\param[in] solver
			\~russian ссылка на метод вычисления сил.
			\~english reference to a method forces evaluation.
			\~
		\param[in] threads
			\~russian пул тредов для параллельного счета.
			\~english thread pool for parallel calculation.
			\~
	*/	
	Base(
		ISolver& solver,
		Threads& threads = tools::dummy_pool
	);
	/** \brief
		\~russian Деструктор класса.
		\~english Class destructor.
		\~
	*/	
	virtual ~Base() = default;

	/** \brief
		\~russian Функция шага по времени. 
		\~english Make time step.
		\~
		\param[in] atoms
			\~russian хранилище с частицами.
			\~english Storage with particles data.
			\~
		\param[in] nlist
			\~russian объект списка соседей для частиц.
			\~english neighbors list object reference.
			\~
		\param[in] dt
			\~russian шаг по времени.
			\~english time step.
			\~
		\details
		\~russian Может включать в себя 
			несколько вызовов функций \p solver, если 
			схема интегрирования по времени многошаговая.
			Также вызывает применение граничных условий.
		\~english May contain several calls to \p solver
			if the integration scheme is a multistep method.
			Also calls boundary conditions application.
		\~
	*/	
	virtual void make_step(
		Atoms&      atoms, 
		NeibsList&  nlist, 
		double      dt
		) = 0;

#ifdef mdcraft_ENABLE_MPI
	/** \brief
		\~russian Функция шага по времени.
		\~english Make time step.
		\~
		\param[in] decomp
			\~russian хранилище с частицами.
			\~english Storage with particles data.
			\~
		\param[in] nlist1
			\~russian объект списка соседей для частиц.
			\~english neighbors list object reference.
			\~
		\param[in] dt
			\~russian шаг по времени.
			\~english time step.
			\~
		\details
		\~russian Может включать в себя
			несколько вызовов функций \p solver, если
			схема интегрирования по времени многошаговая.
			Также вызывает применение граничных условий.
		\~english May contain several calls to \p solver
			if the integration scheme is a multistep method.
			Also calls boundary conditions application.
		\~
	*/
	virtual void make_step(
		Decomp&     decomp,
		NeibsList&  nlist1,
		NeibsList&  nlist2,
		double      dt
	) = 0;
#endif

protected:
	/** \brief
		\~russian Интегрировать по времени одну частицу.
		\~english Time integrate one particle.
		\~
		\param[in] it
			\~russian итератор, обеспечивающий доступ к частице.
			\~english iterator which gives access to a particle.
			\~
		\param[in] dt
			\~russian шаг по времени.
			\~english time step.
			\~
	*/	
	// virtual void step_one(Storage::iterator& it, double dt) = 0;

	ISolver& solver;
	/**<\~russian Аппроксиматор производных по пространству.
		\~english Spatial derivatives approximator. \~ */

	Threads& pool;
		/**<\~russian ссылка на пул тредов для работы функции 
				в многопоточном режиме.
			\~english reference to a pool of threads for 
				multithreaded execution.
			\~
		*/
};

} // namespace mdcraft::solver::stepper

