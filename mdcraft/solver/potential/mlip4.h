#pragma once

#include <string>

#include <mdcraft/configuration.h>

#ifdef mdcraft_ENABLE_MLIP4

#include <mdcraft/solver/potential/base.h>

//#include <external/json-io/vector3.h>
#include <vector3.h>


#include <json.h>
//#include <pots/combined_pot/combined_pot.h>
#include <basic_pots.h>


namespace mdcraft {
namespace solver {
namespace potential {

class MLIP4Potential : public BaseManyBody {
public:
	MLIP4Potential(const std::string& filename, Domain& domain);
	~MLIP4Potential();

	//void prepare( Atoms::iterator atom, Atoms& neibs, const NeibsList& nlist) override;
	void force(   Atoms& atoms, Atoms& neibs, NeibsList& nlist) override;
	// void virial(  Atoms& atoms, Atoms& neibs, NeibsList& nlist) override;

protected:
	struct NBhood {
	public:
		int size() const 
		{
		    return size_;
		}


		int my_ind() const
		{
		    return my_ind_;
		}


		int ind(int j) const 
		{
		    return inds[j];
		}


		int self_species() const 
		{
		    return my_kind_;
		}

		int nbhr_species(int j) const
		{
			return kinds[j];
		}


		Vec3 rel_pos(int j) const
		{
			auto dr = rel_pos_[j];
			return Vec3({dr(0), dr(1), dr(2)});
		}

	public:
		int my_ind_;
		int my_kind_;
		int size_;

		std::vector<int> inds;
		std::vector<int> kinds;

		::mdcraft::data::vector    my_pos;
		std::vector<::mdcraft::data::vector> rel_pos_;
	};

private:
	double rcut_{0.};
    
    std::shared_ptr<AnyLocalPot> mtp_pot{nullptr};
    // std::unordered_map<std::thread::id, std::shared_ptr<AnyLocalPot> p_mlip;

    std::vector<NBhood> frame_atoms;

	void create_nbhoods(Atoms& atoms, Atoms& neibs, NeibsList& nlist);
};

} // potential
} // solver
} // mdcraft

#endif