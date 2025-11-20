#pragma once
#include <functional>
#include <mdcraft/solver/stepper/base.h>

namespace mdcraft::solver::stepper {

/** \brief
	\~russian Класс интегратора движения атомов 
		по времени методом Верле.
	\~english Verlet time integration scheme for atoms.
	\~
*/
class Verlet : public Base {
public:
	/** \brief
		\~russian Конструктор класса.
		\~english Class constructor.
		\~
		\param[in] solver
			\~russian ссылка на метод аппроксимации производных 
				по пространству.
			\~english reference to a method of 
				approximation of spatial derivatives.
			\~
		\param[in] boundaries
			\~russian объект граничных условий.
			\~english boundary conditions object.
			\~
		\param[in] threads
			\~russian пул тредов для параллельного счета.
			\~english thread pool for parallel calculation.
			\~
	*/	
	Verlet(
		ISolver& solver,
		Threads& threads = tools::dummy_pool
	);

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
		\~russian Также вызывает применение граничных условий.
		\~english Also calls boundary conditions application.
		\~
	*/	
	void make_step(
		Atoms&      atoms, 
		NeibsList&  nlist, 
		double      dt
		) override;

#ifdef mdcraft_ENABLE_MPI
	/** \brief
		\~russian Функция шага по времени.
		\~english Make time step.
		\~
		\param[in] decomp
			\~russian
			\~english
			\~
		\param[in] nlist1, nlist2
			\~russian объекты списка соседей для частиц.
			\~english neighbors list object reference.
			\~
		\param[in] dt
			\~russian шаг по времени.
			\~english time step.
			\~
		\details
		\~russian Также вызывает применение граничных условий.
		\~english Also calls boundary conditions application.
		\~
	*/
	void make_step(
		Decomp&     decomp,
		NeibsList&  nlist1,
		NeibsList&  nlist2,
		double      dt
	) override;
#endif

private:

	void step_coords(Atoms& atoms, double dt);
	void step_velocities(Atoms& atoms, double dt);

	/** \brief
		\~russian Интегрировать по времени одну частицу.
		\~english Time integrate one particle.
		\~
		\param[in] atom
			\~russian итератор, обеспечивающий доступ 
				к частице.
			\~english iterator which gives access to 
				a particle.
			\~
		\param[in] dt
			\~russian шаг по времени.
			\~english time step.
			\~
	*/	
	void step_one_coords(Atoms::iterator atom, double dt);

	/** \brief
		\~russian Интегрировать по времени одну частицу.
		\~english Time integrate one particle.
		\~
		\param[in] atom
			\~russian итератор, обеспечивающий доступ 
				к частице.
			\~english iterator which gives access to 
				a particle.
			\~
		\param[in] dt
			\~russian шаг по времени.
			\~english time step.
			\~
	*/	
	void step_one_velocities(Atoms::iterator atom, double dt);
};

} // namespace mdcraft::solver::stepper

