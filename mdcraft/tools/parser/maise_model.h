#pragma once

#include <memory>
#include <vector>

#include <mdcraft/tools/parser/base.h>

namespace mdcraft {
namespace parser {

class MaiseModelHeadContent;
class MaiseModelBodyContent;

class MaiseModelParser: public IParser {
public:
	MaiseModelParser(std::string_view filename);

	virtual ~MaiseModelParser();

	void readHead() override;
    void readBody() override;

    bool headRead() const override;
    bool bodyRead() const override;

    const std::unique_ptr<MaiseModelHeadContent>&  getHead() const;
    const std::unique_ptr<MaiseModelBodyContent>&  getBody() const;

    const std::vector<double>& getWeights() const;

private:
	std::unique_ptr<MaiseModelHeadContent>  Head;
	std::unique_ptr<MaiseModelBodyContent>  Body;
};

struct MaiseModelHeadContent: public IFileContent {
	unsigned 						    id;		// model unique ID
	// atomic info
	int 				   	          nspc;	    // number of species
	std::vector<int> 	  	 		 types;     // species types
	std::vector<std::string> 	 	 names;		// species names
	// NN info
	int 				   		        nl;	    // number of layers
	std::vector<int> 	    			NN;	    // architecture
	std::size_t 					    nw;		// number of weights
};

struct MaiseModelBodyContent: public IFileContent {
	using TableMaps = std::vector<std::vector<std::vector<int>>>;
	using TableVals = std::vector<std::vector<std::vector<double>>>;

	TableMaps  				  	  		table_maps;
	TableVals					  		table_vals;
	std::vector<double> 			     	   b2a;		// bohr to angstrom
	std::vector<double>				 	   scaling;		// scaling
	std::vector<double>              	   weights;     // neural network weights
	std::vector<std::vector<double>>         vRmin;     // input component range: minimum
	std::vector<std::vector<double>>         vRmax;     // input component range: maximum
	std::vector<std::vector<double>>         range;     // input component range: range
};

} // parser
} // mdcraft
