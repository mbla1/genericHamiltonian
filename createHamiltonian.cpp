#include <iostream>
#include <fstream>
#include <vector>
#include "lectureFichiers.h"
#include "Hamiltonian.h"
#include "methodeMatrice.h"

using namespace std;

Potential make_potential(const vector<double>&vec, int grid_size, int elec_size){
	vector<double> pot_vec;
	for(int i = 0; i < grid_size * elec_size; ++i)
		pot_vec.push_back(vec[i]);
	return Potential{grid_size, elec_size, pot_vec};
}

Nac make_nac(const vector<double>& vec, int grid_size, int elec_size){
	vector<Element> nac_vec;
	int cpt = 0; // Number of couplings
	vector<double> temp;
	for(int es1 = 1; es1 < elec_size; ++es1){
		for(int es2 = 0; es2 < es1; ++es2){
			for(int i = 0; i < grid_size; ++i){
				temp.push_back(vec[cpt * grid_size + i]);
			}
			nac_vec.push_back(Element{es1, es2, temp}); // Each coupling pair
			temp.clear(); // Important
			++cpt;
		}
	}
	return Nac{cpt, grid_size, nac_vec};
}

Dipoles make_dipoles(const vector<double>& perm, const vector<double>& trans, int grid_size, int elec_size){
	vector<Element> dip_vec;
	int cpt_perm = 0;
	int cpt_trans = 0;
	vector<double> temp;
	for(int es1 = 0; es1 < elec_size; ++es1){
		for(int es2 = 0; es2 <= es1; ++es2){
			if(es1 == es2){
				for(int i = 0; i < grid_size; ++i)
					temp.push_back(perm[cpt_perm * grid_size + i]);
				++cpt_perm;
			} else {
				for(int i = 0; i < grid_size; ++i)
					temp.push_back(trans[cpt_trans * grid_size + i]);
				++cpt_trans;
			}
			dip_vec.push_back(Element{es1, es2, temp}); // Each coupling pair
			temp.clear(); // Important
		}
	}
	return Dipoles{cpt_perm + cpt_trans, grid_size, dip_vec};
}

int main() {
	constexpr long coord1 = 512; //The order of the coordinates is crucial!!!
    constexpr long coord2 = 1;
    constexpr long elec_size = 4;
    constexpr long grid_size = coord1 * coord2;
    constexpr long size = coord1 * coord2 * elec_size;
    constexpr long hamiltonian_size = size * size;
	constexpr double ang2au = 1.88973;
    constexpr double inv_step_coord1 = 1 / (0.010176125245); //Inverse of NP with the correction from Angströms to Bohr
    constexpr double inv_step_coord2 = 1; //NP
	constexpr double saw2au = 1822.889950851334;
	constexpr double nitrogen_mass = saw2au * 14.003;
	const string output_path = "output/hamiltonianN2.txt";

	cout.precision(12);

	// Potential
	vector<double> potential = prepare_file("source/adiab_sigma_pot512", grid_size, elec_size + 1); // first column is the radius
	Potential pot = make_potential(potential, grid_size, elec_size); // For N2, do not put the potentials to 0

	// Tau
	vector<double> tau = prepare_file("source/adiab_sigma_nac512", grid_size, 4);
	vector<double> temp;
	for(int i = 0; i < grid_size; ++i) // Manual filling because of the data
		temp.push_back(tau[0 * grid_size + i]);
	Element nac12{1, 2, temp};
	temp.clear();

	for(int i = 0; i < grid_size; ++i) 
		temp.push_back(tau[1 * grid_size + i]);
	Element nac13{1, 3, temp};
	temp.clear();

	for(int i = 0; i < grid_size; ++i) 
		temp.push_back(tau[2 * grid_size + i]);
	Element nac23{2, 3, temp};
	vector<Element> vecElem{nac12, nac13, nac23};
	temp.clear();

	Nac nac_vec{3, grid_size, vecElem};

	vector<double> hamiltonien(hamiltonian_size);

	vector<double> gradient1(grid_size * grid_size);

	// Computation of the reduced masses
	double mu1inv = (nitrogen_mass + nitrogen_mass)/(nitrogen_mass * nitrogen_mass); // Reduced mass 
	//double mu1inv = 1./1617.95; //Reduced mass 
	double mu2inv = 0;
	double mu12inv = 0.;

	cout << "mu1inv : " << mu1inv << endl;
	cout << "mu2inv : " << mu2inv << endl;
	cout << "mu12inv : " << mu12inv << endl;

	cout << "Building the kinetic elements" << endl;
	makeKineticOrder4(hamiltonien, elec_size, coord1, coord2, inv_step_coord1, inv_step_coord2, mu1inv, mu2inv, mu12inv);

	cout << "Constructing the gradients" << endl;
	makeGradient6Order1(gradient1, coord1, coord2, inv_step_coord1, mu1inv);

	cout << "Checking out the anti-symmetry" << endl;
	if(isAntiSymmetric(gradient1, grid_size)) {
		cout << "The first gradient is antisymmetric!" <<  endl;
	}

	cout << "Building the scalar products between Tau vectors and gradients" << endl;
	makeNac(hamiltonien, gradient1, nac_vec, grid_size, elec_size);

	cout << "Adding the potential vector" << endl;

	addPotential(hamiltonien, pot, size);	

	cout << "Symmetry test of the hamiltonian" << endl;

	bool testSym = isSymmetric(hamiltonien, size);

	if (!testSym) {
		cout << "This matrix is not symmetrical!" << endl;
		return -1;
	} else {
		cout << "This matrix is symmetrical! Isn't it nice?" << endl;
	}

	cout << "Exporting the hamiltonian in DOK format" << endl;

	bool testWriting = writeSparseForm(hamiltonien, size, output_path);

	if (!testWriting) {
		cout << "An error occured while trying to write the file!" << endl;
		return -1;
	} else {
		cout << "Everything worked correctly! It's done!" << endl;
	}

	return 0;
}
