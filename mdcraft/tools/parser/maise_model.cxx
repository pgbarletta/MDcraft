#include <cassert>

#include <mdcraft/tools/split.h>
#include <mdcraft/tools/parser/maise_model.h>

namespace mdcraft {
namespace parser {

using ::mdcraft::tools::custom_split;

const std::string ICR = "|                          input component range                            |";

auto convert_vector = [] (const std::vector<std::string>& vs) {
	std::vector<int> vi;
	for (auto&& e : vs)
		vi.push_back(std::stoi(e));
	return vi;
};

MaiseModelParser::MaiseModelParser(std::string_view filename)
 : IParser(filename)
 , Head(std::make_unique<MaiseModelHeadContent>())
 , Body(std::make_unique<MaiseModelBodyContent>())
{}

bool MaiseModelParser::headRead() const { return IParser::headRead(); }
bool MaiseModelParser::bodyRead() const { return IParser::bodyRead(); }

void MaiseModelParser::readHead() {
	auto extract_text_value = [this] (std::string& line, std::size_t shift = 10) {
		std::getline(inputFile_, line);
		return line.substr(line.find("|", 1) + 1, line.find("|", 1) + shift);
	};

    std::string line;
    for (std::size_t i = 0; i < 3; ++i)
    	std::getline(inputFile_, line);
   	
	Head->id    = std::stoul(extract_text_value(line), nullptr, 16);
	Head->nspc  = std::stoi(extract_text_value(line));
	Head->types = convert_vector(custom_split(extract_text_value(line)));
	Head->names = custom_split(extract_text_value(line));
	Head->nl    = std::stoi(extract_text_value(line));
	Head->NN    = convert_vector(custom_split(extract_text_value(line)));
	Head->nw    = std::stoi(extract_text_value(line));

	for (; line != ICR; std::getline(inputFile_, line)) {}
	is_head_read_ = true;
}

void MaiseModelParser::readBody() {
	assert(is_head_read_);
	std::string line;
	std::getline(inputFile_, line);
	std::getline(inputFile_, line);

	auto extract_double_value = [this] (std::string& line) {
		std::getline(inputFile_, line);
		return std::stod(line);
	};

	auto convert_vectord = [] (const std::vector<std::string>& vs) {
		std::vector<double> vd;
		for (auto&& e : vs)
			vd.push_back(std::stod(e));
		return vd;
	};

	auto inDim = Head->NN[0];
	Body->vRmin.resize(Head->nspc, std::vector<double>(inDim, 0.0));
	Body->vRmax.resize(Head->nspc, std::vector<double>(inDim, 0.0));
	Body->range.resize(Head->nspc, std::vector<double>(inDim, 0.0));
	for (std::size_t ispc = 0; ispc < Head->nspc; ++ispc)
		for (std::size_t i = 0; i < inDim; ++i) {
			std::getline(inputFile_, line);
			auto dline = convert_vectord(custom_split(line));
			Body->vRmin[ispc][i] = dline[0]; 
			Body->vRmax[ispc][i] = dline[1]; 
			Body->range[ispc][i] = dline[2]; 
		}
	for (std::size_t i = 0; i < 3; ++i) std::getline(inputFile_, line);

	Body->weights.resize(Head->nw);
	for (std::size_t i = 0; i < Head->nw; ++i)
		Body->weights[i] = extract_double_value(line);

	Body->b2a.resize(       Head->nspc);
	Body->scaling.resize(   Head->nspc);
	Body->table_maps.resize(Head->nspc);
	Body->table_vals.resize(Head->nspc);

	for (std::size_t ispc = 0; ispc < Head->nspc; ++ispc) {
		//if we read until the end of file, fill in tables for other species with the same values
		if (inputFile_.peek() == EOF) {
		 	Body->b2a[ispc] = Body->b2a[0];
		 	Body->scaling[ispc] = Body->scaling[0];
		 	Body->table_maps[ispc] = Body->table_maps[0];
		 	Body->table_vals[ispc] = Body->table_vals[0];
		} else {
			for (std::size_t i = 0; i < 3; ++i) std::getline(inputFile_, line);
			auto b2a_sc   = convert_vectord(custom_split(line));
			Body->b2a[ispc]     = b2a_sc[0];
			Body->scaling[ispc] = b2a_sc[1];
			for (std::size_t i = 0; i < 3; ++i) std::getline(inputFile_, line);

			std::vector<std::vector<int>> tmp;
			std::getline(inputFile_, line);
			while (!line.empty()) {
				auto hline = convert_vector(custom_split(line));
				tmp.push_back(hline);
				std::getline(inputFile_, line);
			}
			Body->table_maps[ispc] = tmp;
			
			assert(line.empty());
			std::vector<std::vector<double>> tmpd;
			for (std::size_t i = 0; i < 7; ++i) {
				std::getline(inputFile_, line);
				auto hline = custom_split(line);
				assert(std::stoi(hline[1]) == hline.size() - 2);
				auto hlined = convert_vectord(std::vector(hline.begin() + 2, hline.end()));
				tmpd.push_back(hlined);
			}
			Body->table_vals[ispc] = tmpd;
		}
	}

	is_body_read_ = true;
}

const std::unique_ptr<MaiseModelHeadContent>& MaiseModelParser::getHead() const { return Head; }

const std::unique_ptr<MaiseModelBodyContent>& MaiseModelParser::getBody() const { return Body; }

const std::vector<double>& MaiseModelParser::getWeights() const { return Body->weights; }

MaiseModelParser::~MaiseModelParser() {}

} // parser
} // mdcraft
