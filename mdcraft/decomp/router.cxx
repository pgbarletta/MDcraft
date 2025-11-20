#include <mdcraft/decomp/router.h>

namespace mdcraft::decomp {

using data::Atoms;

inline std::vector<size_t> accumulate(const std::vector<size_t> &arr) {
    std::vector<size_t> res(arr.size());
    res[0] = 0;
    for (size_t i = 1; i < arr.size(); ++i) {
        res[i] = res[i - 1] + arr[i - 1];
    }
    return res;
}

inline size_t sum(const std::vector<size_t> &arr) {
    return std::accumulate(arr.begin(), arr.end(), size_t{0});
}

inline std::ostream &operator<<(std::ostream &os, const std::vector<size_t> &arr) {
    os << "[";
    for (size_t i = 0; i < arr.size() - 1; ++i) {
        os << arr[i] << ", ";
    }
    if (!arr.empty()) {
        os << arr.back() << "]";
    }
    return os;
}

void Requests::wait() {
    for (auto& r: m_requests) {
        if (r != MPI_REQUEST_NULL) {
            MPI_Wait(&r, MPI_STATUS_IGNORE);
            r = MPI_REQUEST_NULL;
        }
    }
}

Router::Router(MPI_Comm comm) : m_comm(comm) {
    MPI_Type_contiguous(sizeof(data::Atom), MPI_BYTE, &MPI_Atom);
    MPI_Type_commit(&MPI_Atom);

    m_send_count  = std::vector<size_t>(m_comm.size(), 0);
    m_send_offset = std::vector<size_t>(m_comm.size(), 0);
    m_recv_count  = std::vector<size_t>(m_comm.size(), 0);
    m_recv_offset = std::vector<size_t>(m_comm.size(), 0);
}

Router::~Router() {
    MPI_Type_free(&MPI_Atom);
}

void Router::set_send_count(const std::vector<size_t> &send_count) {
    assert(m_comm.size() == send_count.size());

    m_send_count  = send_count;
    m_send_offset = accumulate(m_send_count);
}

void Router::set_send_count(const std::atomic_size_t* send_count, int size) {
    assert(m_comm.size() == size);

    m_send_count.resize(m_comm.size());
    for (size_t i = 0; i < m_comm.size(); ++i) {
        m_send_count[i] = send_count[i].load();
    }
    m_send_offset = accumulate(m_send_count);
}

void Router::set_recv_count(const std::vector<size_t> &recv_count) {
    assert(m_comm.size() == recv_count.size());

    m_recv_count  = recv_count;
    m_recv_offset = accumulate(recv_count);

    m_send_recv.clear();
}

void Router::fill_partial() {
    // На получение
    m_comm.all_to_all(m_send_count, m_recv_count);

    // Посчитать смещения
    m_recv_offset = accumulate(m_recv_count);
}

void Router::fill_complete() {
    // Полные обмены числами (все со всеми)
    m_send_recv.resize(m_comm.size() * m_comm.size());

    // Полные обмены числами (все со всеми)
    MPI_Allgather(m_send_count.data(), m_comm.size(), tools::mpi_type<size_t>(),
                  m_send_recv.data(),  m_comm.size(), tools::mpi_type<size_t>(),
                  m_comm);

    // Соберем массив recv_count
    m_recv_count.resize(m_comm.size());
    for (int r = 0; r < m_comm.size(); ++r) {
        m_recv_count[r] = get(r, m_comm.rank());
    }

    // Посчитать смещения
    m_recv_offset = accumulate(m_recv_count);
}

size_t Router::get(int i, int j) const {
    return m_send_recv[m_comm.size() * i + j];
}

size_t Router::send_buffer_size() const {
    return m_send_offset.back() + m_send_count.back();
}

size_t Router::recv_buffer_size() const {
    return m_recv_offset.back() + m_recv_count.back();
}

void Router::print() const {
    if (complete()) {
        print_complete();
    } else {
        print_partial();
    }
}

void Router::print_partial() const {
    std::cout << "Rank " << m_comm.rank() << ". send count: " << m_send_count << ". send offset: " << m_send_offset << "\n";
    std::cout << "        recv count: " << m_recv_count << ". recv offset: " << m_recv_offset << "\n";
}

void Router::print_complete() const {
    int n = 7;
    std::cout << "from \\ to |";
    for (int i = 0; i < m_comm.size(); ++i) {
        std::cout << std::setw(n) << i << " |";
    }
    std::cout << "\n";

    for (int i = 0; i < m_comm.size(); ++i) {
        std::cout << "   " << i << "      |";
        for (int j = 0; j < m_comm.size(); ++j) {
            std::cout << std::setw(n) << get(i, j) << " |";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

Requests Router::isend(const Atoms& src) {
    Requests send_req(m_comm.size());
    for (int r = 0; r < m_comm.size(); ++r) {
        if (m_send_count[r] > 0) {
            MPI_Isend(src.data() + m_send_offset[r], m_send_count[r],
                      MPI_Atom, r, 1371, m_comm, &send_req[r]);
        }
    }
    return send_req;
}

Requests Router::irecv(Atoms& dst) {
    Requests recv_req(m_comm.size());
    for (int r = 0; r < m_comm.size(); ++r) {
        if (m_recv_count[r] > 0) {
            MPI_Irecv(dst.data() + m_recv_offset[r], m_recv_count[r],
                      MPI_Atom, r, 1371, m_comm, &recv_req[r]);
        }
    }
    return recv_req;
}

} // namespace mdcraft::decomp