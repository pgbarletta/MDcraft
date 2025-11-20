#include <mdcraft/decomp/decomp.h>

namespace mdcraft::decomp {

inline std::tuple<double, double> max_avg(const std::vector<double>& loads) {
    double max_load = *std::max_element(loads.begin(), loads.end());
    double avg_load = std::accumulate(loads.begin(), loads.end(), 0.0) / loads.size();

    return {max_load, avg_load};
}

inline void print_loads(const std::vector<double>& loads, bool time_measurer) {
    auto[max_load, avg_load] = max_avg(loads);

    std::cout << std::setprecision(1) << std::fixed;
    if (time_measurer) {
        std::cout << "  Max load: " << 1.0e-3 * max_load << " sec;";
        std::cout << "  Avg load: " << 1.0e-3 * avg_load << " sec;";
    }
    else {
        std::cout << "  Max load: " << 1.0e-3 * max_load << "k;";
        std::cout << "  Avg load: " << 1.0e-3 * avg_load << "k;";
    }
    std::cout << "  Imbalance: " << 100.0 * (max_load / avg_load - 1.0) << "%\n";
    std::cout << "  Loads: [ ";
    for (auto load: loads) {
        std::cout << 100.0 * load / max_load << "% ";
    }
    std::cout << "]\n";
}

Decomp::Decomp(
	Network  comm,
	Atoms&   atoms,
	Threads& pool
)
    : m_locals(atoms),
      m_comm(comm),
      m_pool(pool),
      m_router(comm)
{

}

void Decomp::set_measurer(const std::string& type) {
    if (type != "time" && type != "size") {
        throw std::invalid_argument("Decomposition available measurers: {'time', 'size'}, get unknown '" + type + "'");
    }
    m_time_measurer = type[0] == 't';
}

void Decomp::update(bool verbose) {
    double load = m_locals.size();
    if (m_time_measurer) {
        load = useful.milliseconds();

        useful.stop();
        elapsed.stop();
    }

    auto loads = m_comm.all_gather(load);
    if (verbose && m_comm.master()) {
        print_loads(loads, m_time_measurer);
    }
    balancing(loads);
    redistribute();

    if (m_time_measurer) {
        useful.start();
        elapsed.start();
    }
}

void Decomp::exchange() {
    exchange_start();
    exchange_end();
}

void Decomp::exchange_start() {
    // Prepare array for border
    m_border.resize(m_router.send_buffer_size());

    // Copy particles to border layer
    m_pool.parallel_for(
        size_t{0}, m_border_indices.size(),
        [this](size_t i) {
            m_border[i] = m_locals[m_border_indices[i]];
        });

    m_aliens.resize(m_router.recv_buffer_size());

    m_send_req = m_router.isend(m_border);
    m_recv_req = m_router.irecv(m_aliens);
}

void Decomp::exchange_end() {
    useful.stop(); // Stop waiting

    // Wait for the end of send/recv
    m_send_req.wait();
    m_recv_req.wait();

    useful.resume(); // Launch again
}

void Decomp::exchange_migrants() {
    // Step 1. Get new ranks, count for
    // number of particles to send

    // New ranks of particles
    std::vector<int> ranks(m_locals.size());

    // The number to send
    std::vector<size_t> send_count(m_comm.size(), 0);

    // TODO: multithread version
    for (size_t i = 0; i < m_locals.size(); ++i) {
        int new_rank = rank(m_locals[i].r);

        ranks[i] = new_rank;
        ++send_count[new_rank];
    }

    // Step 2. Get data on message passes
    Router router(m_comm);
    router.set_send_count(send_count);

    // Fill in
    router.fill_partial();

    //router.print();

    // Step 3. Get the permutations indices

    // Array of permutations
    std::vector<size_t> indices(m_locals.size());

    // Shifts of indices
    std::vector<size_t> offsets = router.send_offset();

    // TODO: multithread version
    for (size_t i = 0; i < m_locals.size(); ++i) {
        indices[i] = offsets[ranks[i]]++;
    }

    constexpr bool debug = false;
    if constexpr (debug) {
        // Check the indices run all range from 0 to m_locals.size()
        auto indices2 = indices;
        std::sort(indices2.begin(), indices2.end());
        for (size_t i = 0; i < m_locals.size(); ++i) {
            if (indices2[i] != i) {
                throw std::runtime_error("NO #1205");
            }
        }
    }

    // Step 4. Copy particles to an aux buffer
    Atoms migrants(m_locals.size());
    m_pool.parallel_for(
        size_t{0}, m_locals.size(),
        [this, &migrants, &indices](size_t i) {
            migrants[indices[i]] = m_locals[i];
        });

    // Step 5. Enlarge the buffer to receive data
    m_locals.resize(router.recv_buffer_size());

    // Do message passes
    auto send_req = router.isend(migrants);
    auto recv_req = router.irecv(m_locals);

    send_req.wait();
    recv_req.wait();
}

void Decomp::prepare_border() {
    // Step 1. Evaluate the number to send
    // TODO: multithread version
    std::vector<size_t> send_count(m_comm.size(), 0);
    for (int neib_rank = 0; neib_rank < m_comm.size(); ++neib_rank) {
        if (neib_rank == m_comm.rank()) {
            continue;
        }
        for (auto & atom: m_locals) {
            if (is_near(atom.r, neib_rank)) {
                send_count[neib_rank] += 1;
            }
        }
    }

    // Step 2. Fill in router.send_count

    // Set the number to exchange
    m_router.set_send_count(send_count);

    // Fill in the message passes table
    m_router.fill_partial();

    //m_router.print();

    // Step 3. Fill in border-indices
    m_border_indices.clear();
    m_border_indices.reserve(m_router.send_buffer_size());
    for (int neib_rank = 0; neib_rank < m_comm.size(); ++neib_rank) {
        if (neib_rank == m_comm.rank()) {
            continue;
        }
        for (size_t i = 0; i < m_locals.size(); ++i) {
            if (is_near(m_locals[i].r, neib_rank)) {
                m_border_indices.push_back(i);
            }
        }
    }
}

void Decomp::redistribute() {
    useful.stop();  // Stop timers
    elapsed.stop(); // do not account for redistribution

    exchange_migrants();   // A simple particles redistributor
    collect_locals_info(); // Gather the subdomains info
    prepare_border();      // Indices to send, set the Router

    exchange(); // exchange_start, exchange_end
                // exchange_start: pack the border, isend, irecv
                // exchange_end:   wait for isend and irecv.

    collect_aliens_info(); // Gather the subdomains info

    // Run timers from zero
    useful.start();
    elapsed.start();
}

inline double imbalance(const std::vector<double>& loads) {
    double max_load = *std::max_element(loads.begin(), loads.end());
    double avg_load = std::accumulate(loads.begin(), loads.end(), 0.0) / loads.size();
    return max_load / avg_load - 1.0;
}

void Decomp::prebalancing(int n_iters, bool verbose) {
    std::vector<size_t> count(m_comm.size());
    std::vector<double> loads(m_comm.size());

    for (int k = 0; k < n_iters; ++k) {
        std::fill(count.begin(), count.end(), 0);

        // TODO: multithread version
        for (auto& atom: m_locals) {
            ++count[rank(atom.r)];
        }
        count = m_comm.sum(count);

        for (int r = 0; r < m_comm.size(); ++r) {
            loads[r] = static_cast<double>(count[r]);
        }
        if (verbose && m_comm.master()) {
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "  Prebalancing " << k << " / " << n_iters
                      << ". Imbalance: " << 100 * imbalance(loads) << "%\n";
        }

        balancing(loads);
        redistribute();
    }
}

} // namespace mdcraft::decomp