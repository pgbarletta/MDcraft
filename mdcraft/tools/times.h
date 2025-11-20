#pragma once
#include <ctime>
#include <chrono>

namespace mdcraft::tools {

using wall_clock = std::chrono::system_clock;
	/**<\~russian Системные часы.
		\~english System clock.
		\~
	*/	
using fseconds = std::chrono::duration<float>;
	/**<\~russian Секунды.
		\~english Seconds.
		\~
	*/	

using wall_t = std::chrono::time_point<wall_clock>;
	/**<\~russian Точка во времени, для системных часов.
		\~english Point in time, for system clock.
		\~
	*/	
using proc_t = std::clock_t;
	/**<\~russian Точка во времени, для процесорных часов.
		\~english Point in time, for processor clock ticks time.
		\~
	*/	

using wall_t_elapsed = std::chrono::duration<float>;
	/**<\~russian Время для системных часов измеряется числом с плавающей точкой.
		\~english System time is measured by float variable.
		\~
	*/	
using proc_t_elapsed = std::clock_t;
	/**<\~russian Процессорное время измеряется в тиках.
		\~english System time is measured by float variable.
		\~
	*/


/** \brief
	\~russian Системное время сейчас.
	\~english System time now.
	\~
*/
static wall_t now(wall_t&) { return wall_clock::now(); }

/** \brief
	\~russian Процессорное время сейчас.
	\~english Processor time now.
	\~
*/
static proc_t now(proc_t&) { return std::clock(); }

/** \brief
	\~russian Точка во времени: процессорное и системное время.
	\~english Time point: processor and system time.
	\~
*/
struct time_point_t { proc_t p; wall_t w; };

/** \brief
	\~russian Измерить время сейчас, во временной точке, 
		сначала процессорное, затем системное.
	\~english Find out time now; first \--- processor time, 
		second \--- system time.
	\~
*/
static time_point_t now_wall_and_proc() { time_point_t a; a.w = now(a.w); a.p = now(a.p); return a; }

/** \brief
	\~russian Измерить время сейчас, во временной точке, 
		сначала системное, затем процессорное.
	\~english Find out time now; first \--- system time, 
		second \--- processor time.
	\~
*/
static time_point_t now_proc_and_wall() { time_point_t a; a.p = now(a.p); a.w = now(a.w); return a; }

/** \brief
	\~russian Измерить время сейчас. 
	\~english Find out time now.
	\~
	\see now_wall_and_proc().
*/
static time_point_t now() { return now_wall_and_proc(); }


/** \brief
	\~russian Класс, измеряющий отрезки времени. 
	\~english Class to measure elapsed times (time intervals).
	\~
*/
class elapsed_time_t {
public:
	/** \brief
		\~russian Конструктор класса для нулевого отрезка времени.
		\~english Class constructor for zero elapsed time.
		\~
	*/	
	elapsed_time_t() : p(0), 
					   w(std::chrono::duration<float>::zero()) {}
	/** \brief
		\~russian Конструктор класса для ненулевого отрезка времени.
		\~english Class constructor for non-zero elapsed time.
		\~
		\param[in] p
			\~russian процессорное время.
			\~english process time.
			\~
		\param[in] w
			\~russian системное время.
			\~english system time.
			\~
	*/	
	elapsed_time_t(proc_t_elapsed p, wall_t_elapsed w) : p(p), w(w) {}

	/** \brief
		\~russian Получить процессорное время.
		\~english Get processor time.
		\~
	*/	
	template<class T = fseconds> 
	float proc() {
		return p / (double)CLOCKS_PER_SEC * (double)(typename T::period().den) / (double)(typename T::period().num);
	} 

	/** \brief
		\~russian Получить системное время.
		\~english Get system time.
		\~
	*/	
	template<class T = fseconds> 
	float wall() {
		using dur = std::chrono::duration<float, typename T::period>;
		return std::chrono::duration_cast<dur>(w).count();
	} 

	/** \brief
		\~russian Прибавить отрезки времени.
		\~english Add time interval.
		\~
		\param[in] t
			\~russian отрезок времени.
			\~english time interval.
			\~
	*/	
	elapsed_time_t& operator += (const elapsed_time_t& t) {
		p += t.p;
		w += t.w;
		return *this;
	}
	
	/** \brief
		\~russian Сложить отрезки времени.
		\~english Sum time interval.
		\~
		\param[in] t
			\~russian отрезок времени.
			\~english time interval.
			\~
	*/	
	elapsed_time_t operator + (const elapsed_time_t& t) const {
		elapsed_time_t res;
		res.p = p + t.p;
		res.w = w + t.w;
		return res;
	}

private:
	proc_t_elapsed p;
		/**<\~russian Прошедшее процессорное время.
			\~english Elapsed processor time.
			\~
		*/  
	wall_t_elapsed w;
		/**<\~russian Прошедшее системное время.
			\~english Elapsed system time.
			\~
		*/  
};

/** \brief
	\~russian Получить время между двумя временн'ыми точками.
	\~english Measure time between two time points.
	\~
	\param[in] t1
		\~russian временн'ая точка.
		\~english a time point.
		\~
	\param[in] t2
		\~russian временн'ая точка.
		\~english a time point.
		\~
*/
static elapsed_time_t operator - (const time_point_t& t1, const time_point_t& t2) {
	return elapsed_time_t(t1.p - t2.p, t1.w - t2.w);
}

} // mdcraft::tools