#include <memory>
#include <thread>
#include <cmath>
#include <iostream>

#include <mdcraft/tools/threads.h>

namespace mdcraft::tools {

inline int get_maximum() {
    return static_cast<int>(std::thread::hardware_concurrency());
}

#ifdef mdcraft_ENABLE_MPI
// get_maximum taking into account the number of MPI-processes
inline int get_maximum(MPI_Comm comm) {
    Network net(comm);
    int n_tasks = net.n_tasks();

    // Проблема не решена, округляю в большую сторону,
    // как вариант разделить между процессами, чтобы все
    // точно складывалось. Но как тогда балансировать?
    // Problem not solved. Rounding up. 
    // Possibly could divide among processes, 
    // so it sums up precisely. But how to balance then? 
    int res = std::ceil(get_maximum() / double(n_tasks));

    return res;
}
#endif

Threads::Threads(int count) {
    m_count = 1;
#ifndef mdcraft_ENABLE_MPI
    m_max_count = get_maximum();
#else
    m_max_count = get_maximum(MPI_COMM_WORLD);
#endif
    resize(count);
}

#ifdef mdcraft_ENABLE_MPI
Threads::Threads(int count, MPI_Comm comm) {
    m_count = 1;
    m_max_count = get_maximum(comm);
    resize(count);
}
#endif

void Threads::on() {
    resize(m_max_count);
}

void Threads::off() {
    resize(1);
}

void Threads::resize(int count) {
    m_count = std::max(1, std::min(count, m_max_count));

#ifdef mdcraft_ENABLE_TBB
    if (m_count < 2) {
        m_control = nullptr;
    } else {
        m_control = std::make_shared<tbb::global_control>(tbb::global_control::max_allowed_parallelism, m_count);
    }
#else
    if (m_count < 2) {
        m_pool = nullptr;
    } else {
        m_pool = std::make_shared<ThreadPool>(m_count);
    }
#endif
}

Threads dummy_pool = Threads(1);

} // namespace mdcraft::tools