#pragma once
#include <mpi.h>
#include <mdcraft/data/atom.h>
#include <mdcraft/tools/stopwatch.h>
#include <mdcraft/tools/threads.h>
#include <mdcraft/tools/network.h>
#include <mdcraft/decomp/router.h>

namespace mdcraft::decomp {

using data::Atom;
using data::Atoms;
using tools::Stopwatch;
using tools::Threads;
using tools::Network;

/**
\brief
~english Base class for data decomposition across MPI processes.
~
*/
class Decomp {
public:
	/// @brief
	///   @english Class constructor.
	/// @param comm
	///   @english MPI-Communicator for group of processes
	/// @param atoms
	///   @english Storage with elements for which decomposition is performed.
	/// @param pool
	/// @english Multithreading implementation
	Decomp(
	Network  comm,
	Atoms&   atoms,
	        Threads& pool
	);
	/// @brief Default virtual destructor
	virtual ~Decomp() = default;
	
	/// @brief MPI-rank of the process
	int rank() const { return m_comm.rank(); }
	
	/// @brief
	///   @english Number of MPI-processes
	int size() const { return m_comm.size(); }
	
	/// @brief
	///   @english Reference to elements that are stored by the process.
	Atoms& locals() { return m_locals; }
	
	/// @brief
	///   @english Reference to elements that are not stored by the process,
	///   but are used by it.
	Atoms& aliens() { return m_aliens; }
	
	/// @brief
	///   @english Reference to elements that are sent to other processes during
	///   exchanges, elements can be duplicated.
	Atoms& border() { return m_border; }
	
	/// @brief Set measurer for load balancing. Available types:
	///   "time" -- Useful computational time
	///   "size" -- Number of local atoms
	void set_measurer(const std::string& type);
	
	/// @brief Update decomposition according to input Storage changes.
	void update(bool verbose);
	
	/// @brief Exchange elements near edges of decomposition subdomains.
	void exchange();
	
	/* \brief
	~english Start asynchronous exchange of elements near boundaries
	of decomposition subdomains.
	~
	\see exchange.
	*/
	void exchange_start();
	
	/* \brief
	~english Finish asynchronous exchange of elements near boundaries
	of decomposition subdomains.
	~
	\see exchange.
	*/
	void exchange_end();
	
	/// @brief
	///   @english Perform load balancing. It does not perform actual
	///   redistribution of the elements between MPI-processes.
	virtual void balancing(const std::vector<double>& w) = 0;
	
	/// @brief Perform several load balancing iterations based on the number
	/// of elements in the process. Then redistribute elements.
	void prebalancing(int n_iters = 15, bool verbose = false);
	
	/* \brief
	~english Function to call elements permutation function.
	~
	*/
	void redistribute();

	/// @brief Redistribute particles according to the ranks
	/// that decomposition sets
	void exchange_migrants();
	
	/// @brief Collect global information from local particles
	virtual void collect_locals_info() = 0;
	    
    /// @brief Collect global information from remote particles
    virtual void collect_aliens_info() = 0;
    
    /// @brief Important function. Determine the rank of the process to which
    /// an arbitrary point v belongs. In practice, point v is usually the position
    /// of a particle or the center of a calculation cell.
    /// The function is used for redistribution of elements between processes.
    /// @param v An arbitrary point in the entire domain
    virtual int rank(const data::vector& v) const = 0;
    
    /// @brief Important function. Determine if a certain point v
    /// of this (!) process is in the neighborhood of a process with rank
    /// neib_rank. The function is used for forming aliens lists.
    /// The maximum search radius is usually chosen as the search radius
    /// for elements from this process and the neighboring one (neib_rank)
    /// @param v Point from the subdomain of this process
    /// @param neib_rank Rank of the process in whose neighborhood the point is searched
    /// @return true if the point belongs to the neighborhood
    virtual bool is_near(const data::vector& v, int neib_rank) const = 0;
	
protected:
	    
    /// @brief Collect indices of particles to send,
    void prepare_border();

	Atoms&   m_locals;  ///< internal particles
	Atoms    m_aliens;  ///< particles from the external exchange layer
	Network  m_comm;    ///< MPI-communicator
	Threads& m_pool;    ///< Multithreading implementation
    // Storage for particles to send. Cells that are sent to one
    // process are arranged in a continuous block. A cell can be included
    // in the array twice if sent to multiple processes.
    Atoms m_border;
    // Indices of cells that make up the m_border storage
    std::vector<size_t> m_border_indices;
    Router   m_router;    ///< Routes for sending from m_tourists
    Requests m_send_req;  ///< Requests from isend
	Requests m_recv_req;  ///< Requests from irecv
	
	/// @brief Use computational time for load balancing,
	/// use number of locals atoms if 'false'.
	bool m_time_measurer = false;
	Stopwatch useful;
	
	/**<~english Stopwatch that measures time spent on
	"useful" calculations.
	~
	*/
	Stopwatch elapsed;
	/**<~english Stopwatch that measures all time.
	After calling the \p start() function, it does not stop.
	~
	*/
};
} // namespace mdcraft::decomp