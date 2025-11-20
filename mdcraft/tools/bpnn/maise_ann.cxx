#include <mdcraft/tools/bpnn/maise_ann.h>

#include <iostream>

namespace mdcraft {
namespace bpnn {

MaiseBPNN::MaiseBPNN(std::string_view filename)
 : parser(std::make_unique<MaiseModelParser>(filename))
{
	parser->readHead();
	parser->readBody();
	nKinds_ = parser->getHead()->nspc;

	architec_    = parser->getHead()->NN;
	auto nLayers = architec_.size();

	auto& Body = parser->getBody();

	anns_.resize(nKinds_);

	std::size_t i = 0;
	for (std::size_t kind = 0; kind < nKinds_; ++kind) {

		std::vector<Tens2<double>> biases;
		biases.resize(nLayers - 1);
		for (std::size_t l = 0; l < nLayers - 1; ++l) {
			Tens2<double> biases_l(neurons(l + 1), 1);
			for (std::size_t n = 0; n < architec_[l + 1]; ++n)
				biases_l(n, 0) = Body->weights[i++];
			biases[l] = biases_l;
		}

		std::vector<Tens2<double>> weights;
		weights.resize(nLayers - 1);
		for (std::size_t l = 0; l < nLayers - 1; ++l) {
			Tens2<double> weights_l(neurons(l + 1), neurons(l));
			for (std::size_t nnext = 0; nnext < architec_[l + 1]; ++nnext)
				for (std::size_t ncur = 0; ncur < architec_[l]; ++ncur)
					weights_l(nnext, ncur) = Body->weights[i++];
			weights[l] = weights_l;
		}

		/* see maise nutl.c lines 621-630 */
		for (std::size_t in = 0; in < neurons(0); ++in)
			for (std::size_t k = 0; k < neurons(1); ++k)
				weights[0](k, in) *= Body->range[kind][in];
		/* actually vRmin always zeros in all potentials in maise examples */
		for (std::size_t in = 0; in < neurons(0); ++in)
			Body->vRmin[kind][in] /= Body->range[kind][in];

		anns_[kind] = {biases, weights};
	}

}

const std::vector<std::vector<double>>&  MaiseBPNN::vRmin() const {
	return parser->getBody()->vRmin;
}

const std::vector<std::vector<double>>&  MaiseBPNN::vRmax() const {
	return parser->getBody()->vRmax;
}

const std::vector<std::vector<double>>&  MaiseBPNN::range() const {
	return parser->getBody()->range;
}

} // bpnn
} // mdcraft