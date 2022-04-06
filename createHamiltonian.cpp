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
    constexpr double inv_step_coord1 = 1 / (0.040625 * ang2au); //Inverse of NP with the correction from Angströms to Bohr
    constexpr double inv_step_coord2 = 1; //NP
	double saw2au = 1822.889950851334;
	double nitrogen_mass = saw2au * 6.941;

	cout.precision(12);

	// Potential
	vector<double> potential = prepare_file("source/adiab_sigma_pot512", grid_size, elec_size + 1); // first column is the radius
	Potential pot = make_potential(potential, grid_size, elec_size);

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


	// Transition dipoles
	vector<double> t_dip = prepare_file("source/adiab_sigma_mom512", grid_size, 4); // Strange set of data
	vector<double> p_dip(size); // Zero-filled vector
	Dipoles dip_vec = make_dipoles(p_dip, t_dip, grid_size, elec_size);
	
	// Allocating the memory to make the code a bit faster
	vector<double> hamiltonien;
	hamiltonien.reserve(hamiltonian_size);

	vector<double> gradient1;
	gradient1.reserve(grid_size * grid_size);

	return 0;
/*
	//Computation of the reduced masses
	double mu1inv = (nitrogen + nitrogen)/(nitrogen * nitrogen); //Reduced mass 
	//double mu1inv = 1./1617.95; //Reduced mass 
	double mu2inv = 0;
	double mu12inv = 0.;

	cout << "mu1inv : " << mu1inv << endl;
	cout << "mu2inv : " << mu2inv << endl;
	cout << "mu12inv : " << mu12inv << endl;

	//Kinetic energy
	cout << "Building the kinetic elements" << endl;
	makeKineticOrder4(hamiltonien, elec_size, coord1, coord2, inv_step_coord1, inv_step_coord2, mu1inv, mu2inv, mu12inv);

	//Gradients
	cout << "Constructing the gradients" << endl;
	makeGradient6Order1(gradient1, coord1, coord2, inv_step_coord1, mu1inv);
	//makeGradient1(gradient1, coord1, coord2, inv_step_coord1, mu1inv);

	cout << "Checking out the anti-symmetry" << endl;
	if(isAntiSymmetric(gradient1, grid_size)) {
		cout << "The first gradient is antisymmetric!" <<  endl;
	}

	return 0;

	//Tau vectors
	cout << "Building Tau vectors" << endl;

	//Scalar product between the Tau vectors and gradients
	cout << "Building the scalar products between Tau vectors and gradients" << endl;
	//To Loïc: those should alternate between the coordinates (q1 and q2) and the electronic couplings: (0,1),(0,2),(1,2)

	makeNac(nac21, gradient1, vec21, grid_size);
	makeNac(nac31, gradient1, vec31, grid_size);
	makeNac(nac32, gradient1, vec32, grid_size);
	makeNac(nac41, gradient1, vec41, grid_size);
	makeNac(nac42, gradient1, vec42, grid_size);
	makeNac(nac43, gradient1, vec43, grid_size);
	makeNac(nac51, gradient1, vec51, grid_size);
	makeNac(nac52, gradient1, vec52, grid_size);
	makeNac(nac53, gradient1, vec53, grid_size);
	makeNac(nac54, gradient1, vec54, grid_size);
	makeNac(nac61, gradient1, vec61, grid_size);
	makeNac(nac62, gradient1, vec62, grid_size);
	makeNac(nac63, gradient1, vec63, grid_size);
	makeNac(nac64, gradient1, vec64, grid_size);
	makeNac(nac65, gradient1, vec65, grid_size);
	makeNac(nac71, gradient1, vec71, grid_size);
	makeNac(nac72, gradient1, vec72, grid_size);
	makeNac(nac73, gradient1, vec73, grid_size);
	makeNac(nac74, gradient1, vec74, grid_size);
	makeNac(nac75, gradient1, vec75, grid_size);
	makeNac(nac76, gradient1, vec76, grid_size);

	cout << "Adding the NAC" << endl;

	addNac(hamiltonien, nac21, 1, 0, size, grid_size);
	addNac(hamiltonien, nac31, 2, 0, size, grid_size);
	addNac(hamiltonien, nac32, 2, 1, size, grid_size);
	addNac(hamiltonien, nac41, 3, 0, size, grid_size);
	addNac(hamiltonien, nac42, 3, 1, size, grid_size);
	addNac(hamiltonien, nac43, 3, 2, size, grid_size);
	addNac(hamiltonien, nac51, 4, 0, size, grid_size);
	addNac(hamiltonien, nac52, 4, 1, size, grid_size);
	addNac(hamiltonien, nac53, 4, 2, size, grid_size);
	addNac(hamiltonien, nac54, 4, 3, size, grid_size);
	addNac(hamiltonien, nac61, 5, 0, size, grid_size);
	addNac(hamiltonien, nac62, 5, 1, size, grid_size);
	addNac(hamiltonien, nac63, 5, 2, size, grid_size);
	addNac(hamiltonien, nac64, 5, 3, size, grid_size);
	addNac(hamiltonien, nac65, 5, 4, size, grid_size);
	addNac(hamiltonien, nac71, 6, 0, size, grid_size);
	addNac(hamiltonien, nac72, 6, 1, size, grid_size);
	addNac(hamiltonien, nac73, 6, 2, size, grid_size);
	addNac(hamiltonien, nac74, 6, 3, size, grid_size);
	addNac(hamiltonien, nac75, 6, 4, size, grid_size);
	addNac(hamiltonien, nac76, 6, 5, size, grid_size);

	cout << "Loading the potential energy files" << endl;

	double potential[size];

	string nomsFichiersPot[elec_size] = {"1", "2", "3", "4", "5", "6", "7"}; //Autant de potentiels que d'etats electroniques

	bool testPotential = loadPotential(potential, grid_size, "source/pot", nomsFichiersPot, ".txt", elec_size);

	if (!testPotential) {
		cout << "An error occured while reading the potential files!" << endl;
		return -1;
	}

	cout << "Adding the potential vector" << endl;

	addPotential(hamiltonien, potential, size);	

	cout << "Symmetry test of the hamiltonian" << endl;

	bool testSym = isSymmetric(hamiltonien, size);

	if (!testSym) {
		cout << "This matrix is not symmetrical!" << endl;
		return -1;
	} else {
		cout << "This matrix is symmetrical! Isn't it nice?" << endl;
	}

	cout << "Exporting the hamiltonian in DOK format" << endl;

	bool testWriting = writeSparseForm(hamiltonien, size, "hamiltonienLiHAlt4O.txt");

	if (!testWriting) {
		cout << "An error occured while trying to write the file!" << endl;
		return -1;
	} else {
		cout << "Everything worked correctly! It's done!" << endl;
		return 0;
	}*/
}
