#include <iostream>
#include <vector>
#include <array>
#include <complex>
#include "Hamiltonian.h"

using namespace std;

// Methods for elements

ostream& operator<<(ostream& ost, Element el){
	return ost << "(" << el.es_init << "," << el.es_fin << ")";
}

Element Element::opposite() const{ // For the NAC, where the order matters
	vector<double> temp;
	for(double d : values){
		temp.push_back(-d);
	}

	return Element{es_fin, es_init, temp};
}

Element Element::symmetric() const{ // For the dipole, where the order does not matter
	return Element{es_init, es_fin, values};
}

double Element::operator[](int index) const{
	return values[index];	
}

/* The potential energy vector. It is simply diagonal, so it behaves
 * much like a wavefunction. I simply define access operators
 * It is not directly based on elements (though it could be?).
 */

double Potential::operator[](int index) const{
	return coeff[index]; 
}

double Potential::operator()(int grid_index, int elec_index) const{
	return coeff[elec_index * grid_size + grid_index];
}

/* A Nac is a vector of Elements. As there are more couplings
 * than electronic states, the access is a bit different:
 * you seach inside the vector the Element that represents your 
 * coupling. The Nac are antisymmetric, thus 
 * search(es1, es2) = -search(es2, es1) in terms of values
 */

Element Nac::search(int es1, int es2) const{
	for(Element el : couplings){
		if(el.es_init == es1 && el.es_fin == es2)
			return el;
		else if(el.es_init == es2 && el.es_fin == es1)
			return el.opposite();
	}
	cerr << "Element not found\n";
}

/* A structure for the Dipoles, which is similar to Nac.
 * The main difference is that it is supposed to contain
 * the permanent dipoles as well (and have both indices
 * (1,1) and (2,1)), and that the dipoles are symmetric.
 */

Element Dipoles::search(int es1, int es2) const{
	for(Element el : couplings){
		if(es1 == es2 && el.es_init == es1) return el; // Might not impact much the performance
		else if(el.es_init == es1 && el.es_fin == es2) return el;
		else if(el.es_init == es2 && el.es_fin == es1)
			return el.symmetric();
	}
	cerr << "Element not found\n";
}

/* A structure to describe the wavefunction and
 * access to its elements.
 */

complex<double> Wavefunction::operator[](int index) const{
	return coeff[index]; 
}

complex<double> Wavefunction::operator()(int grid_index, int elec_index) const{
	return coeff[elec_index * grid_size + grid_index];
}

Wavefunction operator+(const Wavefunction& wf1, const Wavefunction& wf2){
	Wavefunction somme{wf1};
	int grid_size = wf1.grid_size;
	int elec_size = wf1.elec_size;

	for(int i = 0; i < grid_size * elec_size; ++i)
		somme.coeff[i] = wf1[i] + wf2[i];
	
	return somme;
}

Wavefunction operator*(const Wavefunction& wf, double d){
	Wavefunction mult{wf};
	for(int i = 0; i < mult.grid_size * mult.elec_size; ++i)
		mult.coeff[i] *= d;

	return mult;
}

Wavefunction operator*(double d, const Wavefunction& wf){
	Wavefunction mult{wf};
	for(int i = 0; i < mult.grid_size * mult.elec_size; ++i)
		mult.coeff[i] *= d;

	return mult;
}

Wavefunction operator/(const Wavefunction& wf, complex<double> cd){
	Wavefunction div{wf};
	for(int i = 0; i < div.grid_size * div.elec_size; ++i)
		div.coeff[i] /= cd;

	return div;
}

Wavefunction operator/(const Wavefunction& wf, double d){
	Wavefunction div{wf};
	for(int i = 0; i < div.grid_size * div.elec_size; ++i)
		div.coeff[i] /= d;

	return div;
}
