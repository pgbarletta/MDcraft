#include <mdcraft/tools/network.h>

namespace mdcraft::tools {

#ifdef mdcraft_ENABLE_MPI

MPI_Datatype MPI_VEC3 = MPI_DATATYPE_NULL;

inline void init_types() {
    if (MPI_VEC3 == MPI_DATATYPE_NULL) {
        MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_VEC3);
        MPI_Type_commit(&MPI_VEC3);
    }
}

inline void free_types() {
    if (MPI_VEC3 != MPI_DATATYPE_NULL) {
        MPI_Type_free(&MPI_VEC3);
    }
}

Network::Network(int& argc, char**& argv) {
    int flag;
    MPI_Initialized(&flag);
    if (!flag) {
        MPI_Init(&argc, &argv);
    }

    m_comm = MPI_COMM_WORLD;

    MPI_Comm_size(m_comm, &m_size);
    MPI_Comm_rank(m_comm, &m_rank);

    init_types();
}

Network::Network(MPI_Comm comm) {
    int flag;
    MPI_Initialized(&flag);
    if (!flag) {
        MPI_Init(nullptr, nullptr);
    }

    m_comm = comm;

    MPI_Comm_size(comm, &m_size);
    MPI_Comm_rank(comm, &m_rank);

    init_types();
}

void Network::finalize() {
    free_types();
    MPI_Finalize();
}

std::vector<std::string> Network::proc_names() {
    const int max_length = MPI_MAX_PROCESSOR_NAME;

    int length = 0;
    char name[max_length];

    MPI_Get_processor_name(name, &length);

    std::vector<int> lengths = all_gather(length);
    std::vector<int> offsets(lengths.size());
    offsets[0] = 0;
    for (int i = 1; i < m_size; ++i) {
        offsets[i] = offsets[i - 1] + lengths[i - 1];
    }

    int buff_size = offsets.back() + lengths.back();
    char *all_names = new char[buff_size + 1];

    MPI_Allgatherv(name, length, MPI_CHAR, all_names,
                   lengths.data(), offsets.data(), MPI_CHAR, m_comm);

    all_names[buff_size] = '\0';

    std::vector<std::string> names(m_size);
    for (int r = 0; r < m_size; ++r) {
        names[r] = std::string(all_names + offsets[r], lengths[r]);
    }

    //std::cout << "Process names:\n";
    //for (int r = 0; r < m_size; ++r) {
    //    std::cout << "  '" << names[r] << "'\n";
    //}

    return names;
}

int Network::n_tasks() {
    auto names = proc_names();
    std::string my_name = names[m_rank];

    int n_tasks = 0;
    for (int r = 0; r < m_size; ++r) {
        if (names[r] == my_name) {
            n_tasks += 1;
        }
    }

    return n_tasks;
}
#endif

} // namespace mdcraft::tools
