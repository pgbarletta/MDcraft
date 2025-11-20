#include <random>
#include <mdcraft/decomp/VD3.h>

namespace mdcraft::decomp {

VD3::VD3(
	Network             comm,
    Atoms&              atoms,
	Domain&             domain,
	Threads&            pool,
	params              opts,
	std::vector<vector> centers
) : Decomp(comm, atoms, pool),
	m_domain(domain),
	m_dim(std::min(opts.dimension, domain.dim())),
	m_mobility(opts.mobility),
	m_centroidal(opts.centroidal),
	m_growth_rate(opts.growth_rate) {
    if (centers.size() == m_comm.size()) { // use user initial positions
        m_centers = centers;
    } else { // fill random within domain
        m_centers.resize(m_comm.size());
        m_weights.resize(m_comm.size());

        std::mt19937 gen(57391);
        std::uniform_real_distribution x(m_domain.xmin(), m_domain.xmax());
        std::uniform_real_distribution y(m_domain.ymin(), m_domain.ymax());
        std::uniform_real_distribution z(m_domain.zmin(), m_domain.zmax());

        for (int i = 0; i < m_comm.size(); ++i) {
            vector probe;
            do {
                probe.x() = x(gen);
                probe.y() = m_dim > 1 ? y(gen) : 0.0;
                probe.z() = m_dim > 2 ? z(gen) : 0.0;
            } while (!m_domain.belong(probe));

            m_centers[i] = probe;
            m_weights[i] = 0.0;
        }
    }

    // init neibs
    m_neibs.resize(m_comm.size());
    // add all ranks as neibs
    for (int i = 0; i < m_comm.size(); ++i) {
        m_neibs[i] = i;
    }

    synchronize();
    redistribute();
}

vector VD3::center(int i) const {
	return m_centers[i];
}

std::vector<vector> VD3::centers() const {
	return m_centers;
}

void VD3::balancing(const std::vector<double>& w) {
    correct_cell_center(w);
    synchronize();
}

inline double imb_func(double Imax) {
    // The value of disbalance which provides 1/2
    // an acceptable parameter value
    const double I0 = 1.0e-2;
    return Imax / (Imax + I0);
}

void VD3::correct_cell_center(const std::vector<double>& loads) {
    const double xi_x = m_mobility;
    const double sigma = m_centroidal;

    vector dr = {0.0, 0.0, 0.0};

    double min_imb = 1.0e20;
    double max_imb = 0.0;

    int iGen = m_comm.rank();
    for (auto jGen: m_neibs) {
        if (iGen == jGen) {
            continue;
        }
        double imb = (loads[jGen] - loads[iGen]) / (loads[jGen] + loads[iGen]);

        // Check for NaN, emerges when loads[iGen] = loads[jGen] = 0.0
        if (std::isnan(imb)) { imb = 0.0; }

        vector dir = direction(iGen, jGen).normalized();

        dr += imb * dir;

        min_imb = std::min(min_imb, std::abs(imb));
        max_imb = std::max(max_imb, std::abs(imb));
    }

    // Renew the coordinates

    // Minimum distance to the cell boundary
    double DR = m_min_radius[iGen];

    double theta_x = xi_x * imb_func(max_imb);

    vector new_coords = center(iGen) + theta_x * DR * dr;

    // Shift to the center of mass
    new_coords = sigma * m_centroids[iGen] + (1.0 - sigma) * new_coords;

    // The generator into the domain
    m_domain.fit_in_period(new_coords);

    // Set the new coordinate
    m_centers[iGen] = vector::Zero();
    for (int i = 0; i < m_dim; ++i) {
        m_centers[iGen][i] = new_coords[i];
    }
}

void VD3::move_center(int i, const vector& v) {

}

void VD3::synchronize() {
    double w = m_weights[m_comm.rank()];
    vector c = m_centers[m_comm.rank()];

    m_comm.all_gather(w, m_weights);
    m_comm.all_gather(c, m_centers);
}

double VD3::weight(int i) const {
    return m_weights[i];
}

vector VD3::direction(int i, int j) const {
    return m_domain.shortest(center(j) - center(i));
}

double VD3::distance(int i, int j) const {
    return m_domain.shortest(center(i) - center(j)).norm();
}

double VD3::wdistance(int i, const vector& v) const {
    return (m_domain.shortest(center(i) - v)).norm();
}

int VD3::rank(const vector& v) const {
    int rank = 0;
    double dist = std::numeric_limits<double>::max();
    for (int i = 0; i < m_comm.size(); ++i) {
        double d = wdistance(i, v);
        if (d < dist) {
            dist = d;
            rank = i;
        }
    }
    return rank;
}

bool VD3::is_near(const vector& v_in, int neib_rank) const {
    // Neibs search radius
    double search_rad = std::max(
            m_search_radius[m_comm.rank()],
            m_search_radius[neib_rank]
    );

    /*
    // Domains radii
    double R1 = m_max_radius[m_comm.rank()];
    double R2 = m_max_radius[neib_rank];

    // Far-away procs
    double dist = distance(m_comm.rank(), neib_rank);
    if (dist > R1 + R2 + 2 * search_rad) { return false; }
    */

    // Need to transform the point and neighboring generator into the coordinate system
    // associated with the own generator considering periodicity. In this system, the own
    // generator will have coordinates (0, 0, 0).

    vector g = m_centers[neib_rank] - m_centers[m_comm.rank()];
    vector v = v_in - m_centers[m_comm.rank()];

    // Учтем периодичность. Областью периодичности считаем |x| < Lx/2, |y| < Ly/2, |z| < Lz/2.
    // Сам генератор будет располагаться в центре. Таким образом, даже с учетом периодов ячейка
    // Вороного для нашего процесса будет выпуклой.
    for (int i = 0; i < m_dim; ++i) {
        if (m_domain.periodic(i)) {
            if (v[i] < -0.5 * m_domain.size(i)) {
                v[i] += m_domain.size(i);
            }
            else if (v[i] > 0.5 * m_domain.size(i)) {
                v[i] -= m_domain.size(i);
            }

            if (g[i] < -0.5 * m_domain.size(i)) {
                g[i] += m_domain.size(i);
            }
            else if (g[i] > 0.5 * m_domain.size(i)) {
                g[i] -= m_domain.size(i);
            }
        }
    }

    // Without periodicity, a single check is sufficient:
    // return g.normalized().dot(0.5 * g - v) <= search_rad;
    // This line checks the proximity of the point to the neighboring generator g.
    // For periodic boundaries, it is necessary to check not only the proximity of the point to generator g
    // but also to its periodic images. However, for the check, not all possible images need to be considered,
    // only those that "surround" our own generator.
    // For example, in the 2D case with periodicity, it is sufficient to check 4 images of generator g
    // that form a rectangle around the coordinate origin (our own generator).
    // In the 3D case, it is sufficient to check 8 generators that form a parallelepiped
    // around our own generator.

    // These are the coordinates of the proposed images of generator g
    // { {x1, x2}, {y1, y2}, {z1, z2} }
    // The value NAN indicates that the check with such a generator image should not be performed.
    std::array<std::array<double, 2>, 3> coords {
            std::array<double, 2>{g.x(), NAN},
            std::array<double, 2>{g.y(), NAN},
            std::array<double, 2>{g.z(), NAN}
    };

    for (int i = 0; i < m_dim; ++i) {
        // If |g[i]| is small, skip the image check. This would occur if two generators
        // shared nearly the same coordinate, which is theoretically possible but highly unlikely.
        // As a precaution, you could add a check: && std::abs(g[i]) > search_rad

        // Add the image coordinate g if the boundary is periodic.
        if (m_domain.periodic(i)) {
            coords[i][1] = g[i] > 0.0 ? g[i] - m_domain.size(i) : g[i] + m_domain.size(i);
        }
    }

    // {x1, x2}
    for (double x: coords[0]) {
        if (std::isnan(x)) { continue; }

        // {y1, y2}
        for (double y: coords[1]) {
            if (std::isnan(y)) { continue; }

            // {z1, z2}
            for (double z: coords[2]) {
                if (std::isnan(z)) { continue; }

                // Изображение точки g
                vector G = {x, y, z};
                if (G.normalized().dot(0.5 * G - v) <= search_rad) {
                    return true;
                }
            }
        }
    }
    return false;
}

void VD3::collect_locals_info() {
    // NativeInfo
    // Contains information about a single cell
    // max_rad  -- distance from the element to the center of the subdomain, after reduce - radius
    // of the circle circumscribed around the subdomain.
    // search_r -- neighbor search radius, after reduce - maximum neighbor
    // search radius.
    // centroid -- cell coordinates, after reduce -- sum of coordinates of all
    // elements in the subdomain
    struct NInfo {
        double max_rad  = 0.0;
        double search_r = 0.0;
        vector centroid = {0.0, 0.0, 0.0};

        /// reduce function
        void operator&=(const NInfo& res) {
            max_rad = std::max(max_rad, res.max_rad);
            search_r = std::max(search_r, res.search_r);
            centroid += res.centroid;
        }
    };

    vector center = m_centers[m_comm.rank()];
    auto func = [this, &center](Atom& atom) -> NInfo {
        vector dr = m_domain.shortest(atom.r - center);
        return {.max_rad  = dr.norm(),
                .search_r = atom.rns,
                .centroid = dr};
    };

    NInfo res = m_pool.reduce(m_locals.begin(), m_locals.end(), NInfo{}, func);

    res.centroid /= m_locals.size();
    res.centroid += center;

    // gather from all the procs
    m_comm.all_gather(res.max_rad,  m_max_radius);
    m_comm.all_gather(res.search_r, m_search_radius);
    m_comm.all_gather(res.centroid, m_centroids);
}

void VD3::collect_aliens_info() {
    // list of the neighbour procs
    m_neibs.clear();
    for (int r = 0; r < m_comm.size(); ++r) {
        if (m_router.send_count(r) > 0) {
            m_neibs.push_back(r);
        }
    }

    // find the inscribed circle in a Voronoi cell
    vector center = m_centers[m_comm.rank()];
    auto func = [this, &center](Atom& atom) -> double {
        return (m_domain.shortest(atom.r - center)).norm();
    };

    double minR = m_pool.min<2>(m_aliens.begin(), m_aliens.end(), func);

    m_comm.all_gather(minR, m_min_radius);
}

} // namespace mdcraft::decomp
