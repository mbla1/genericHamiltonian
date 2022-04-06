#include <vector>
#include <array>
#include <complex>

/* The smallest subdivision of the Hamiltonian:
 * a vector coupling electronic states with values on the grid.
 * The potentials, NACs and dipoles will be based on those Elements.
 */

struct Element{
	int es_init;
	int es_fin;
	std::vector<double> values;

	Element opposite() const;
	Element symmetric() const;
	double operator[](int index) const;
};

std::ostream& operator<<(std::ostream& ost, Element el);

/* The potential energy vector. It is simply diagonal, so it behaves
 * much like a wavefunction. I simply define access operators
 * It is not directly based on elements (though it could be?).
 */

struct Potential{
	int grid_size;
	int elec_size;
	std::vector<double> coeff;

	double operator[](int index) const;
	double operator()(int grid_index, int elec_index) const;
};

/* A Nac is a vector of Elements. As there are more couplings
 * than electronic states, the access is a bit different:
 * you seach inside the vector the Element that represents your 
 * coupling. The Nac are antisymmetric, thus 
 * search(es1, es2) = -search(es2, es1) in terms of values
 */

struct Nac{
	int n_couplings;
	int grid_size;
	std::vector<Element> couplings;

	Element search(int es1, int es2) const;
};

/* A structure for the Dipoles, which is similar to Nac.
 * The main difference is that it is supposed to contain
 * the permanent dipoles as well (and have both indices
 * (1,1) and (2,1)), and that the dipoles are symmetric.
 */

struct Dipoles{
	int n_couplings;
	int grid_size;
	std::vector<Element> couplings;

	Element search(int es1, int es2) const;
};

/* A structure to describe the wavefunction and
 * access to its elements.
 */

struct Wavefunction{ // Add more features: norm, dot product, ...
	int grid_size;
	int elec_size;
	std::vector<std::complex<double>> coeff;

	std::complex<double> operator[](int index) const;
	std::complex<double> operator()(int grid_index, int elec_index) const;
};

Wavefunction operator+(const Wavefunction& wf1, const Wavefunction& wf2);
Wavefunction operator*(const Wavefunction& wf, double d);
Wavefunction operator*(double d, const Wavefunction& wf);
Wavefunction operator/(const Wavefunction& wf, double d);
Wavefunction operator/(const Wavefunction& wf, std::complex<double> cd);
