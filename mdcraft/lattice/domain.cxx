#include <mdcraft/lattice/domain.h>

#include <iostream>

namespace mdcraft::lattice {

using data::vector;

Domain::Domain(
	double xmin, double xmax,
	double ymin, double ymax,
	double zmin, double zmax
) : m_bmin{xmin, ymin, zmin},
    m_bmax{xmax, ymax, zmax},
    m_size{xmax - xmin, ymax - ymin, zmax - zmin},
    m_periodic{false, false, false} {

	m_dim = std::fabs(zmin - zmax) > 1e-15 ? 3 :
			std::fabs(ymin - ymax) > 1e-15 ? 2 : 1;
}

void Domain::reshape(
	double xmin, double xmax,
	double ymin, double ymax,
	double zmin, double zmax
) {
	m_bmin = {xmin, ymin, zmin};
	m_bmax = {xmax, ymax, zmax};
	m_size = {xmax - xmin, ymax - ymin, zmax - zmin};
	m_dim = zmin == zmax ? (ymin == ymax ? 1 : 2) : 3;
}

bool Domain::belong(vector v) const {
	bool inside = true;
	for (int i = 0; i < 3; i++) {
		inside = inside && (m_bmin[i] <= v[i] && v[i] <= m_bmax[i]);
	}
	return inside;
}

void Domain::fit_in_period(vector& v) const {
	for (int i = 0; i < 3; i++) {
		if (m_periodic[i])
			v[i] = v[i] < m_bmin[i] ? v[i] + m_size[i] :
				       v[i] > m_bmax[i] ? v[i] - m_size[i] : v[i];
	}
}

double Domain::nearest_distance(vector v, int axis) const {
	double dist = std::numeric_limits<double>::max();
	if (m_periodic[axis]) {
		double dist_min = std::fabs(v[axis] - m_bmin[axis]);
		double dist_max = std::fabs(v[axis] - m_bmax[axis]);
		dist = std::min(dist_min, dist_max);
	}
	return dist;
}

vector Domain::shortest(vector v) const {
	for (int i = 0; i < 3; ++i) {
		if (!m_periodic[i]) continue;
		v[i] = std::fabs(v[i] + m_size[i]) < std::fabs(v[i]) ? v[i] + m_size[i] : v[i];
		v[i] = std::fabs(v[i] - m_size[i]) < std::fabs(v[i]) ? v[i] - m_size[i] : v[i];
	}
	return v;
}

} // namespace mdcraft::lattice