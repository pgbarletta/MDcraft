#pragma once
#include <mdcraft/lattice/domain.h>
#include <mdcraft/decomp/decomp.h>
namespace mdcraft::decomp {
using lattice::Domain;
using data::vector;
/** ~russian Класс декомпозиции данных по MPI-процессам
на основе подвижной диаграммы Вороного.
~english Data decomposition among MPI processes based
on moving Voronoi diagram.
~
*/
class VD3 : public Decomp {
public:
    /// @brief Default values for balancing parameters
    struct params {
         /// @brief Decomposition dimension {1, 2, 3}
        int dimension = 1;
        /// @brief Generator displacement speed (0, 1)
        double mobility = 0.2;
        /// @brief Parameter determining the fraction of displacement toward the geometric center
        /// in the generator displacement expression.
        double centroidal = 0.25;
        /// @brief Generator weight change rate (0, 1)
        double growth_rate = 0.02;
    };
    VD3(
     Network             comm,
     Atoms&              atoms,
     Domain&             domain,
     Threads&            pool,
     params              opts,
     std::vector<vector> centers = {}
    );
    vector center(int i) const;
    std::vector<vector> centers() const;
    /* \brief
    ~russian Выполнить балансировку нагрузки.
    ~english Perform load balancing.
    ~
    */
    void balancing(const std::vector<double>& w) override;
    /* ~russian Сместить центр подобласти диаграммы.
    ~english Move the subdomain center.
    ~
    \param[in] w
    ~russian вектор нагрузки всех подобластей.
    ~english load vector for all subdomains.
    ~
    */
    void correct_cell_center(const std::vector<double>& w);
    /// @brief Move the i-th generator by vector v. The new value
    /// is placed inside the domain.
    void move_center(int i, const vector& v);
    /// @brief Weight of the i-th generator
    double weight(int i) const;
    /// @brief Vector in the direction from the i-th generator to the j-th generator
    /// taking periodicity into account.
    vector direction(int i, int j) const;
    /// @brief Standard distance between two generators taking
    /// periodicity into account.
    double distance(int i, int j) const;
    /// @brief Weighted "distance" from the i-th generator to point v
    /// taking periodicity into account
    double wdistance(int i, const vector& v) const;
    /// @brief Synchronize the diagram across all processes
    void synchronize();
    /// @brief Determine the rank of the process to which an arbitrary
    /// point v belongs. Used for redistributing particles among processes.
    /// @param v An arbitrary point in the entire domain
    int rank(const data::vector& v) const final;
    /// @brief Determine if a certain point v of this (!) process
    /// lies in the neighborhood of the process with rank neib_rank. The function is used for
    /// forming aliens lists. The search radius is selected as
    /// the maximum search radius for particles from this process and the neighboring one.
    /// @param v_in Point from the subdomain of this process
    /// @param neib_rank Rank of the process in whose neighborhood the point is searched
    /// @return true if the point belongs to the neighborhood
    bool is_near(const data::vector& v_in, int neib_rank) const final;
    /// @brief Collect general information from local particles
    ///     vector<double> m_search_radius; -- neighbor search radii
    ///     vector<double> m_max_radius;    -- circumscribed circle radii
    ///     vector<vector> m_centroids;     -- Voronoi cell centroids
    void collect_locals_info() final;
    /// @brief Collect global information from remote particles
    void collect_aliens_info() final;
private:
    Domain&             m_domain;   ///< Computational domain
    std::vector<vector> m_centers;  ///< Diagram generators
    std::vector<double> m_weights;  ///< Diagram weights
    /// @brief Decomposition dimension {1, 2, 3}
    int m_dim;
    /// @brief Generator displacement speed (0, 1)
    double m_mobility;
    /// @brief Parameter determining the fraction of displacement toward the geometric center
    /// in the generator displacement expression.
    double m_centroidal;
    /// @brief Weight change rate (0, 1)
    double m_growth_rate;
    /// Subdomain information
    /// @brief List of neighbors of the current process, i.e., ranks of processes,
    /// with which boundary element data exchange occurs
    std::vector<int> m_neibs;
    /// @brief Maximum neighbor search radius
    std::vector<double> m_search_radius;
    /// @brief Inscribed circle radius
    std::vector<double> m_min_radius;
    /// @brief Circumscribed circle radius
    std::vector<double> m_max_radius;
    /// @brief Subdomain barycenter
    std::vector<vector> m_centroids;
};
} // namespace mdcraft::decomp