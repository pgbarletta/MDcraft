#pragma once

#include <Eigen/Eigen>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/CXX11/Tensor>

namespace mdcraft {
namespace bpnn {
namespace autodiff {


template <typename T>
using Tens2 = Eigen::Tensor<T, 2>;

template <typename T>
using TM2 = Eigen::TensorMap<Tens2<T>>;

using TM2d = TM2<double>;

template <typename T>
using ADiff = Eigen::AutoDiffScalar<Eigen::Matrix<T, Eigen::Dynamic, 1>>;

using ADiffd = ADiff<double>;

template <typename T>
T sigmoid(T t) {
	return 1 / (exp(-t) + 1);
}

template <typename T>
Tens2<T> sigmoid_af(Tens2<T>& P) {
	Tens2<T> res = P.unaryExpr(std::ref(sigmoid<T>));
	return res;
}

template <typename T>
T tanhyp(T t) {
	return tanh(t);
}

template <typename T>
Tens2<T> tanh_af(Tens2<T>& P) {
	Tens2<T> res = P.unaryExpr(std::ref(tanhyp<T>));
	return res;
}


template <typename T>
Tens2<T> linear_af(Tens2<T>& P) {
	return P;
}

/* not necessary to make it a template function
   it can be made as a template class instead
   'In' and 'Bias' are assumed to be given as Tens2 of shapes (1, N)
	E.g.,
	      In = { {G0, G1, ..., G144} },
	      Bias = { {B0, ..., B9} },
		  W = { {w0.0, w0.1, ..., w0.9},
		        {w1.0, w1.1, ..., w1.9},
		                ...
		     {w144.0, w144.1, ... w144.9} }
*/
template <typename T>
Tens2<T> Layer(
	Tens2<T>& In,
	Tens2<T>& Bias,
	Tens2<T>& W,
	std::function<Tens2<T>(Tens2<T>&)> activation
) {
	Eigen::array<Eigen::IndexPair<Eigen::Index>, 1> dims = {
		Eigen::IndexPair<Eigen::Index>(1, 0)
	};
	Tens2<T> Z = (W.contract(In, dims) + Bias);
	Tens2<T> res = activation(Z);
	return res;
}

template <typename T>
Tens2<ADiff<T>> prepare(Eigen::TensorMap<Tens2<T>> tensor, int offset = 0, int size = 0) {
	const int rows = tensor.dimension(0);
	const int cols = tensor.dimension(1);

	Tens2<ADiff<T>> res(rows, cols);

	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j) {
			int index = i * cols + j;
			res(i,j).value() = tensor(i,j);
			if (size)
				res(i,j).derivatives() = Eigen::VectorXd::Unit(size, offset + index);
		}

	return res;
};

template <typename T>
Tens2<ADiff<T>> prepare(const Tens2<T>& tensor, int offset = 0, int size = 0) {
	const int rows = tensor.dimension(0);
	const int cols = tensor.dimension(1);

	Tens2<ADiff<T>> res(rows, cols);

	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j) {
			int index = i * cols + j;
			res(i,j).value() = tensor(i,j);
			if (size)
				res(i,j).derivatives() = Eigen::VectorXd::Unit(size, offset + index);
		}

	return res;
};

template <typename T>
Tens2<T> grads(const ADiff<T>& Out, const Tens2<ADiff<T>>& In) {
	auto derivatives = Out.derivatives();
	int index = 0;

	Tens2<T> res(In.dimension(0), In.dimension(1));
	for (int i = 0; i < In.dimension(0); ++i)
		for (int j = 0; j < In.dimension(1); ++j) {
			T val = derivatives[index];
			res(i, j) = val;
			index++;
		}

	return res;
};


} // autodiff
} // bpnn
} // mdcraft