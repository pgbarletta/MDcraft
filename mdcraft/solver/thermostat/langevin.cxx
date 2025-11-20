#include <mdcraft/solver/thermostat/langevin.h>

namespace mdcraft {
namespace solver {
namespace thermostat {

Langevin::Langevin(
	double beta,
	double temperature,
	double time_step,
	int    heat_x,
	int    heat_y,
	int    heat_z,
	double Ux,
	double Uy,
	double Uz,
	double xmin,
	double xmax
) : Base(),
    beta(beta), T(temperature), dt(time_step), bxmin(xmin), bxmax(xmax),
    heat_axis(heat_x, heat_y, heat_z),
    U(Ux, Uy, Uz)
{
	sigma = std::sqrt(2.0*beta*Kb*T/dt);
	generator.seed(random_dev());
	distribution = std::normal_distribution<double>(0.0, sigma);
}

void Langevin::apply_one(Atoms::iterator atomit) {
	auto& atom = *atomit;
	if ((atom.r.x() < bxmin) || (atom.r.x() > bxmax)) return;
	vector random_force(
		distribution(generator), 
		distribution(generator), 
		distribution(generator)
	);
	const double sqrt_mass = std::sqrt(atom.m);
	vector therm_force = random_force*sqrt_mass;
	therm_force -= (atom.v - U)*beta*atom.m;
	therm_force = therm_force.array() * heat_axis.array();
	atom.f += therm_force;
}

void Langevin::set_average_velocity(vector _U) { U = _U; }

void Langevin::set_friction(double _beta) {
	beta = _beta;
	sigma = std::sqrt(2.0*beta*Kb*T/dt);
	generator.seed(random_dev());
	distribution = std::normal_distribution<double>(0.0, sigma);
}

void Langevin::set_temperature(double _T) {
	T = _T;
	sigma = std::sqrt(2.0*beta*Kb*T/dt);
	generator.seed(random_dev());
	distribution = std::normal_distribution<double>(0.0, sigma);	
}

void Langevin::set_heating_axes(
	int heat_x,
	int heat_y,
	int heat_z
) {
	heat_axis = vector(heat_x, heat_y, heat_z);
}

void Langevin::set_dimensions(double _xmin, double _xmax) {
	bxmin = _xmin;
	bxmax = _xmax;
}

double Langevin::friction() { return beta;}
double Langevin::temperature() { return T;}
vector Langevin::average_velocity() { return U;}
double Langevin::xmin() { return bxmin; }
double Langevin::xmax() { return bxmax; }

} // thermostat
} // solver
} // mdcraft