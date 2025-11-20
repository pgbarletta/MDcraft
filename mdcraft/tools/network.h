#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include <mdcraft/configuration.h>

#ifdef mdcraft_ENABLE_MPI
#include <mpi.h>
#include <mdcraft/data/vector.h>

namespace mdcraft::tools {

/// @brief Simplified interface for basic MPI
class Network {
private:
    int      m_rank;  ///< Rank of MPI process
    int      m_size;  ///< Number of processes
    MPI_Comm m_comm;  ///< Network communicator
    
public:
    /// @{
    /// @name Basic functions

    /// @brief Classic MPI-initialization
    Network(int& argc, char**& argv);
    
    /// @brief Create from MPI_Comm
    Network(MPI_Comm comm = MPI_COMM_WORLD);

    /// @brief Explicit conversion to MPI_Comm
    operator MPI_Comm() { return m_comm; }

    /// @brief Free types and finalize MPI
    void finalize();

    // @brief Rank of MPI process
    int rank() const { return m_rank; }

    /// @brief Number of MPI processes
    int size() const { return m_size; }

    /// @brief Is it master process?
    bool master() const { return m_rank == 0; }

    /// @brief Single process execution?
    bool single() const { return m_size == 1; }

    /// @brief MPI communicator
    MPI_Comm comm() const { return m_comm; };

    /// @brief 
    void barrier() const { MPI_Barrier(m_comm); }

    /** \brief
        \~russian Выполнить функцию последовательно (!) на каждом процессе
        \~english Execute function consecutively (!) on each process. 
        \~
    */  
    template <class Func>
    void for_each(Func&&) const;

    /// @brief Collect unique names of nodes
    std::vector<std::string> proc_names();

    /// @brief Number of MPI tasks per node
    int n_tasks();

    /// @}

    /// @{
    /// @name Collective reduce operations

    /** \brief
        \~russian Коллективная операция. Минимальное значение по всем процессам.
        \~english Collective operation. Minimum value over all processes. 
        \~
    */  
    template <class T>
    T min(const T& value);

    /** \brief
        \~russian Коллективная операция. Покомпонентный минимум для каждого 
        элемента вектора среди всех процессов сети.
        \~english Collective operation. Componentwise minimum 
        for each vector element among all processes.
        \~
    */  
    template <class T>
    std::vector<T> min(const std::vector<T>& values);

    /** \brief
        \~russian Коллективная операция. Максимальное значение по всем процессам.
        \~english Collective operation. Maximum value over all processes. 
        \~
    */  
    template <class T>
    T max(const T& value);

    /** \brief
        \~russian Коллективная операция. Покомпонентный максимум для каждого 
        элемента вектора среди всех процессов сети.
        \~english Collective operation. Componentwise maximum 
        for each vector element among all processes.
        \~
    */  
    template <class T>
    std::vector<T> max(const std::vector<T>& values);

    /** \brief
        \~russian Коллективная операция. Сумма величин со всех процессов сети.
        \~english Collective operation. Sum of values from all network processes.
        \~
    */  
    template <class T>
    T sum(const T& value);


    /** \brief
        \~russian Коллективная операция. Покомпонентная сумма элементов
        векторов со всех процессов сети.
        \~english Collective operation. Componentwise sum of 
        vector elements from all network processes.
        \~
    */  
    template <class T>
    std::vector<T> sum(const std::vector<T>& values);

    /** \brief
        \~russian Номер процесса и значение.
        \~english Process number and value.
        \~
    */  
    template <class T>
    struct value_rank { T value; int rank; };

    /** \brief
        \~russian Минимальное значение и номер процесса, на котором достигается.
        \~english Minimum value and number of process, where it is reached.
        \~
    */  
    template <class T>
    value_rank<T> argmin(const T& value);

    /** \brief
        \~russian Максимальное значение и номер процесса, на котором достигается.
        \~english Maximum value and number of process, where it is reached.
        \~
    */ 
    template <class T>
    value_rank<T> argmax(const T& value);

    /// @}

    /// @{
    /// @name Collective communication operations

    /** 
        \~russian 
        @brief Отправляет сообщение с "корневого" процесса всем процессам.
        @param root Ранг "корневого" процесса.
        @param value Величина шаблонного типа по ссылке, переменная также
        принимает пересланное значение.
        \~english 
        @brief Sends a message from "root" process to all processes.
        @param root "Root" process rank.
        @param value Size of template type by reference. 
        The variable also accepts the sent value.
        \~
    */  
    template <class T>
    void bcast(int root, T& value);

    /** \brief
        \~russian Коллективная операция. Собирает величину value со всех
        процессов, записывает в массив и раздает этот массив всем процессам.
        \~english Collective opration. Collects "value" from all processes,
        writes to an array, sends this array to all processes.
        \~
    */ 
    template <class T>
    std::vector<T> all_gather(const T& value);

    /** \brief
        \~russian Коллективная операция. Собирает величину value со всех процессов,
        записывает в массив и раздает этот массив всем процессам сети.
        \~english Collective opration. Collects "value" from all processes,
        writes to an array, sends this array to all network processes.
        \~
    */ 
    template <class T>
    void all_gather(const T& value, std::vector<T>& values);

    /** \brief
        \~russian Коллективная операция. Размеры буферов send и recv совпадают
        и равны числу процессов. Каждая величина из send отправляется своему
        процессу, каждая величина в recv получается с определенного процесса.
        \~english Collective opration. Sizes of buffers send and recv coincide 
        and are equal in terms of processes number. Every value from send is sent to 
        its own process, each value in recv is received from a certain process.
        \~
    */ 
    template <class T>
    void all_to_all(const std::vector<T>& send, std::vector<T>& recv);

    /// @}
};

/** \brief
    \~russian Примитивные MPI-типы по шаблону
    По умолчанию используются типы contiguous, кратные 4 байтам
    \~english Primitive MPI-types by template.
    By default contiguous types, divisible by 4 bytes, are used.
    \~
*/ 
template<class T>
MPI_Datatype mpi_type();

extern MPI_Datatype MPI_VEC3;

template <> inline MPI_Datatype mpi_type<int>()    { return MPI_INT;    }
template <> inline MPI_Datatype mpi_type<double>() { return MPI_DOUBLE; }
template <> inline MPI_Datatype mpi_type<float>()  { return MPI_FLOAT;  }
template <> inline MPI_Datatype mpi_type<long>()   { return MPI_LONG;   }
template <> inline MPI_Datatype mpi_type<short>()  { return MPI_SHORT;  }
template <> inline MPI_Datatype mpi_type<char>()   { return MPI_BYTE;   }
template <> inline MPI_Datatype mpi_type<size_t>() { return MPI_UNSIGNED_LONG; }
template <> inline MPI_Datatype mpi_type<data::vector> () { return MPI_VEC3; }


template <class Func>
void Network::for_each(Func&& f) const {
    if (single()) {
        f(); std::cout.flush();
        return;
    }
    for (int r = 0; r < m_size; ++r) {
        if (r == m_rank) {
            f(); std::cout.flush();
        }
        barrier();
    }
}

template <class T>
T Network::min(const T& value) {
    if (single()) {
        return value;
    }
    T g_value;
    MPI_Allreduce(&value, &g_value, 1, mpi_type<T>(), MPI_MIN, m_comm);
    return g_value;
}

template <class T>
std::vector<T> Network::min(const std::vector<T>& values) {
    if (single()) {
        return values;
    }
    std::vector<T> g_values(values.size());
    MPI_Allreduce(values.data(), g_values.data(), int(values.size()),
                  mpi_type<T>(), MPI_MIN, m_comm);
    return g_values;
}

template <class T>
T Network::max(const T& value) {
    if (single()) {
        return value;
    }
    T g_value;
    MPI_Allreduce(&value, &g_value, 1, mpi_type<T>(), MPI_MAX, m_comm);
    return g_value;
}

template <class T>
std::vector<T> Network::max(const std::vector<T>& values) {
    if (single()) {
        return values;
    }
    std::vector<T> g_values(values.size());
    MPI_Allreduce(values.data(), g_values.data(), int(values.size()),
                  mpi_type<T>(), MPI_MAX, m_comm);
    return g_values;
}

template <class T>
T Network::sum(const T& value) {
    if (single()) {
        return value;
    }
    T g_value;
    MPI_Allreduce(&value, &g_value, 1, mpi_type<T>(), MPI_SUM, m_comm);
    return g_value;
}

template <class T>
std::vector<T> Network::sum(const std::vector<T>& values) {
    if (single()) {
        return values;
    }
    std::vector<T> g_values(values.size());
    MPI_Allreduce(values.data(), g_values.data(), int(values.size()),
                  mpi_type<T>(), MPI_SUM, m_comm);
    return g_values;
}

template <class T>
Network::value_rank<T> Network::argmin(const T& value) {
    static_assert(std::is_same<T, double>::value, "Only doubles");
    value_rank<T> in {.value=value, .rank=m_rank};
    value_rank<T> out{.value=value, .rank=m_rank};
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    return out;
}

template <class T>
Network::value_rank<T> Network::argmax(const T& value) {
    static_assert(std::is_same<T, double>::value, "Only doubles");
    value_rank<T> in {.value=value, .rank=m_rank};
    value_rank<T> out{.value=value, .rank=m_rank};
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    return out;
}

template <class T>
void Network::bcast(int root, T& value) {
    MPI_Bcast((void*)&value, sizeof(value), MPI_CHAR, root, m_comm);
}

template <class T>
std::vector<T> Network::all_gather(const T& value) {
    std::vector<T> values(m_size);
    MPI_Allgather(&value, 1, mpi_type<T>(), values.data(), 1, mpi_type<T>(), m_comm);
    return values;
}

template <class T>
void Network::all_gather(const T& value, std::vector<T>& values) {
    values.resize(m_size);
    MPI_Allgather(&value, 1, mpi_type<T>(), values.data(), 1, mpi_type<T>(), m_comm);
}

template <class T>
void Network::all_to_all(const std::vector<T>& send, std::vector<T>& recv) {
    recv.resize(m_size);
    MPI_Alltoall(send.data(), 1, mpi_type<T>(), recv.data(), 1, mpi_type<T>(), m_comm);
}

} // namespace mdcraft::tools

#else

// Однопроцессорная заглушка при компиляции без mpi
// Single-process stub for compilation without MPI

namespace mdcraft::tools {

/// @brief Stub for MPI interface
class Network {
public:
    /// @{
    /// @name Basic functions

    Network() = default;

    /// @brief Classic MPI-initialization
    Network(int& argc, char**& argv) { }

    /// @brief Free types and finalize MPI
    void finalize() const { }

    // @brief Rank of MPI process
    constexpr int rank() const { return 0; }

    /// @brief Number of MPI processes
    constexpr int size() const { return 1; }

    /// @brief Is it master process?
    constexpr bool master() const { return true; }

    /// @brief Single process execution?
    constexpr bool single() const { return true; }

    /// @brief
    void barrier() const { }

    /** \brief
    \~russian Выполнить функцию последоватьно (!) на каждом процессе
    \~english Execute function consecutively (!) on each process.
    \~
    */ 
    template <class Func>
    void for_each(Func&& func) const { func(); }

    /// @brief Number of MPI tasks per node
    constexpr int n_tasks() const { return 1; }

    /// @}

    /// @{
    /// @name Collective reduce operations

    /** \brief
        \~russian Коллективная операция. Минимальное значение по всем процессам.
        \~english Collective operation. Minimum value over all processes. 
        \~
    */      
    template <class T>
    T min(const T& value) const { return value; }

    /** \brief
        \~russian Коллективная операция. Покомпонентный минимум для каждого 
        элемента вектора среди всех процессов сети.
        \~english Collective operation. Componentwise minimum 
        for each vector element among all processes.
        \~
    */  
    template <class T>
    std::vector<T> min(const std::vector<T>& values) const {
        return values;
    }

    /** \brief
        \~russian Коллективная операция. Максимальное значение по всем процессам.
        \~english Collective operation. Maximum value over all processes. 
        \~
    */     
    template <class T>
    T max(const T& value) const { return value; }

    /** \brief
        \~russian Коллективная операция. Покомпонентный максимум для каждого 
        элемента вектора среди всех процессов сети.
        \~english Collective operation. Componentwise maximum 
        for each vector element among all processes.
        \~
    */  
    template <class T>
    T max(const std::vector<T>& values) const { return values; }

    /** \brief
        \~russian Коллективная операция. Сумма величин со всех процессов сети.
        \~english Collective operation. Sum of values from all network processes.
        \~
    */    
    template <class T>
    T sum(const T& value) const { return value; }

    /** \brief
        \~russian Коллективная операция. Покомпонентная сумма элементов
        векторов со всех процессов сети.
        \~english Collective operation. Componentwise sum of 
        vector elements from all network processes.
        \~
    */  
    template <class T>
    std::vector<T> sum(const std::vector<T>& values) const { return values; }

    /** \brief
        \~russian Номер процесса и значение.
        \~english Process number and value.
        \~
    */  
    template <class T>
    struct value_rank { T value; int rank; };

    /** \brief
        \~russian Минимальное значение и номер процесса, на котором достигается.
        \~english Minimum value and number of process, where it is reached.
        \~
    */      
    value_rank<double> argmin(double value) const {
        return {.value=value, .rank=0};
    }

    /** \brief
        \~russian Максимальное значение и номер процесса, на котором достигается.
        \~english Maximum value and number of process, where it is reached.
        \~
    */     
    value_rank<double> argmax(double value) const {
        return {.value=value, .rank=0};
    }

    /// @}

    /// @{
    /// @name Collective communication operations

    /** 
        \~russian 
        @brief Отправляет сообщение с "корневого" процесса всем процессам.
        @param root Ранг "корневого" процесса.
        @param value Величина шаблонного типа по ссылке, переменная также
        принимает пересланное значение.
        \~english 
        @brief Sends a message from "root" process to all processes.
        @param root "Root" process rank.
        @param value Size of template type by reference. 
        The variable also accepts the sent value.
        \~
    */  
    template <class T>
    void bcast(int root, T& value) const { }

    /** \brief
        \~russian Коллективная операция. Собирает величину value со всех
        процессов, записывает в массив и раздает этот массив всем процессам.
        \~english Collective opration. Collects "value" from all processes,
        writes to an array, sends this array to all processes.
        \~
    */ 
    template <class T>
    std::vector<T> all_gather(const T& value) const { return {value}; }

    /** \brief
        \~russian Коллективная операция. Собирает величину value со всех процессов,
        записывает в массив и раздает этот массив всем процессам сети.
        \~english Collective opration. Collects "value" from all processes,
        writes to an array, sends this array to all network processes.
        \~
    */ 
    template <class T>
    void all_gather(const T& value, std::vector<T>& values) const {
        values[0] = value;
    }

    /** \brief
        \~russian Коллективная операция. Размеры буферов send и recv совпадают
        и равны числу процессов. Каждая величина из send отправляется своему
        процессу, каждая величина в recv получается с определенного процесса.
        \~english Collective opration. Sizes of buffers send and recv coincide 
        and are equal in terms of processes number. Every value from send is sent to 
        its own process, each value in recv is received from a certain process.
        \~
    */     
    template <class T>
    void all_to_all(const std::vector<T>& send, std::vector<T>& recv) const {
        recv = send;
    }

    /// @}
};

} // namespace mdcraft::tools

#endif