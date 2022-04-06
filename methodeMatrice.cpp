#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "methodeMatrice.h"

using namespace std;

void setZero(double vecteur[], long dimension) {
	for (long i = 0; i < dimension; i++) {
		vecteur[i] = 0.;
	}	
}

double calculMoments(double coordinates1[], double coordinates2[], double masses[])
{
	double somme = 0;

    for(int iter = 0; iter < 15; iter++) {
		somme += coordinates1[iter] * coordinates2[iter] / masses[iter / 3];	
    }
    return somme;
}

void makeKineticOrder4(vector<double>& matrice, long elecStates, long coordinate1,
			   	long coordinate2, double step1, double step2, double invReducedMass, double mom2, double mom12) {
	long taille = elecStates * coordinate1 * coordinate2;
	//First coordinate
	for (int i = 0; i < elecStates; i++) {
        for (int j = 0; j < coordinate1; j++) {
            for (int k = 0; k < coordinate2; k++) {
                matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille 
						+ (i * coordinate1 * coordinate2 + j * coordinate2 + k)] +=
                        (-1 / 2.) * invReducedMass * step1 * step1 * (-5.) / 2.; //j=0
                if (j > 0) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + (j - 1) * coordinate2 + k)] +=
                            (-1 / 2.) * invReducedMass * step1 * step1 * (4.) / 3.; // j-1
                }
                if (j > 1) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + (j - 2) * coordinate2 + k)] +=
                            (-1 / 2.) * invReducedMass * step1 * step1 / (-12.); //j-2
                }
                if (j < coordinate1 - 1) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + (j + 1) * coordinate2 + k)] +=
                            (-1 / 2.) * invReducedMass * step1 * step1 * (4.) / 3.; //j+1
                }
                if (j < coordinate1 - 2) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + (j + 2) * coordinate2 + k)] +=
                            (-1 / 2.) * invReducedMass * step1 * step1 / (-12.); //j+2
                }
            }
        }
    }

	//Second derivative in q2 : c(a,b+1) = M[a*coordinate2+b,a*coordinate2+b+1]
    for (int i = 0; i < elecStates; i++) {
        for (int j = 0; j < coordinate1; j++) {
            for (int k = 0; k < coordinate2; k++) {
                matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille 
						+ (i * coordinate1 * coordinate2 + j * coordinate2 + k)] +=
                        (-1 / 2.) * mom2 * step2 * step2 * (-5.) / 2.;
                if (k > 0) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + j * coordinate2 + k - 1)] +=
                            (-1 / 2.) * mom2 * step2 * step2 * (4.) / 3.;
                }
                if (k > 1) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + j * coordinate2 + k - 2)] +=
                            (-1 / 2.) * mom2 * (-step2 * step2) / 12.;
                }
                if (k < coordinate2 - 1) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + j * coordinate2 + k + 1)] +=
                            (-1 / 2.) * mom2 * step2 * step2 * (4.) / 3.;
                }
                if (k < coordinate2 - 2) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + j * coordinate2 + k + 2)] +=
                            (-1 / 2.) * mom2 * (-step2 * step2) / 12.;
                }
            }
        }
    }

	//Cross-derivative : c(a+1,b+1) = M[a*coordinate2+b,(a+1)*coordinate2+b+1]
	//This contribution is present twice within the hamiltonian, hence the difference
	// with the formula in the article
    for (int i = 0; i < elecStates; i++) {
        for (int j = 0; j < coordinate1; j++) {
            for (int k = 0; k < coordinate2; k++) {
                matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille 
						+ (i * coordinate1 * coordinate2 + j * coordinate2 + k)] +=
                        (-1 / 2.) * mom12 * step1 * step2 * 2.;
                if (j > 0 && k > 0) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + (j - 1) * coordinate2 + k - 1)] +=
                            (-1 / 2.) * mom12 * step1 * step2;
                }
                if (j < coordinate1 - 1 && k < coordinate2 - 1) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + (j + 1) * coordinate2 + k + 1)] +=
                            (-1 / 2.) * mom12 * step1 * step2;
                }
                if (j > 0) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + (j - 1) * coordinate2 + k)] +=
                            (-1 / 2.) * mom12 * (-step1 * step2);
                }
                if (k > 0) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + j * coordinate2 + k - 1)] +=
                            (-1 / 2.) * mom12 * (-step1 * step2);
                }
                if (j < coordinate1 - 1) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + (j + 1) * coordinate2 + k)] +=
                            (-1 / 2.) * mom12 * (-step1 * step2);
                }
                if (k < coordinate2 - 1) {
                    matrice[(i * coordinate1 * coordinate2 + j * coordinate2 + k) * taille +
                                (i * coordinate1 * coordinate2 + j * coordinate2 + k + 1)] +=
                            (-1 / 2.) * mom12 * (-step1 * step2);
                }
            }
        }
    }
}

void makeGradient1(vector<double>& matrice, long coordinate1, long coordinate2, double step1, double invReducedMass) {
	long tailleGrille = coordinate1 * coordinate2;
	for (int j = 0; j < coordinate1; j++) {
        for (int k = 0; k < coordinate2; k++) {
            if (j > 0) {
                matrice[(j * coordinate2 + k) * tailleGrille + ((j - 1) * coordinate2 + k)] += -invReducedMass * step1 * (-2.) / 3.;  //j-1 pas de 0.5 !!!
            }
            if (j > 1) {
                matrice[(j * coordinate2 + k) * tailleGrille + ((j - 2) * coordinate2 + k)] += -invReducedMass * (step1) * (1.) / 12.; //j-2
            }
            if (j < coordinate1 - 1) {
                matrice[(j * coordinate2 + k) * tailleGrille + ((j + 1) * coordinate2 + k)] += -invReducedMass * step1 * (2.) / 3.; //j+1
            }
            if (j < coordinate1 - 2) {
                matrice[(j * coordinate2 + k) * tailleGrille + ((j + 2) * coordinate2 + k)] += -invReducedMass * (step1) * (-1.) / 12.; //j+2
            }
        }
    }
}

void makeGradient2(vector<double>& matrice, long coordinate1, long coordinate2, double step2) {
	long tailleGrille = coordinate1 * coordinate2;
	for (int j = 0; j < coordinate1; j++) {
        for (int k = 0; k < coordinate2; k++) {
            if (k > 0) {
                matrice[(j * coordinate2 + k) * tailleGrille + (j * coordinate2 + k - 1)] += step2 * (-2.) / 3.;
            }                                                                             
            if (k > 1) {                                                                  
                matrice[(j * coordinate2 + k) * tailleGrille + (j * coordinate2 + k - 2)] += (step2) * (1.) / 12.;
            }                                                                             
            if (k < coordinate2 - 1) {                                                         
                matrice[(j * coordinate2 + k) * tailleGrille + (j * coordinate2 + k + 1)] += step2 * (2.) / 3.;
            }                                                                             
            if (k < coordinate2 - 2) {                                                         
                matrice[(j * coordinate2 + k) * tailleGrille + (j * coordinate2 + k + 2)] += (step2) * (-1.) / 12.;
            }
        }
    }
}

void makeGradient6Order1(vector<double>& matrice, long coordinate1, long coordinate2, double step1, double reducedMass) {
	long tailleGrille = coordinate1 * coordinate2;
	for (int j = 0; j < coordinate1; j++) {
        for (int k = 0; k < coordinate2; k++) {
            if (j > 0) {
                matrice[(j * coordinate2 + k) * tailleGrille + ((j - 1) * coordinate2 + k)] += -reducedMass * step1 * (-3.) / 4.; //j-1
            }
            if (j > 1) {
                matrice[(j * coordinate2 + k) * tailleGrille + ((j - 2) * coordinate2 + k)] += -reducedMass * step1 * (3.) / 20.; //j-2
            }
            if (j > 2) {
                matrice[(j * coordinate2 + k) * tailleGrille + ((j - 3) * coordinate2 + k)] += -reducedMass * step1 * (-1.) / 60.; //j-3
            }
            if (j < coordinate1 - 1) {
                matrice[(j * coordinate2 + k) * tailleGrille + ((j + 1) * coordinate2 + k)] += -reducedMass * step1 * (3.) / 4.; //j+1
            }
            if (j < coordinate1 - 2) {
                matrice[(j * coordinate2 + k) * tailleGrille + ((j + 2) * coordinate2 + k)] += -reducedMass * step1 * (-3.) / 20.; //j+2
            }
            if (j < coordinate1 - 3) {
                matrice[(j * coordinate2 + k) * tailleGrille + ((j + 3) * coordinate2 + k)] += -reducedMass * step1 * (1.) / 60.; //j+3
            }
        }
    }
}

void makeGradient6Order2(vector<double>& matrice, long coordinate1, long coordinate2, double step2) {
	long tailleGrille = coordinate1 * coordinate2;
	for (int j = 0; j < coordinate1; j++) {
        for (int k = 0; k < coordinate2; k++) {
            if (k > 0) {
                matrice[(j * coordinate2 + k) * tailleGrille + (j * coordinate2 + k - 1)] += step2 * (3.) / 4.;
            }                                                                             
            if (k > 1) {                                                                  
                matrice[(j * coordinate2 + k) * tailleGrille + (j * coordinate2 + k - 2)] += (step2) * (-3.) / 20.;
            }                                                                             
            if (k > 2) {                                                                  
                matrice[(j * coordinate2 + k) * tailleGrille + (j * coordinate2 + k - 3)] += (step2) * (1.) / 60.;
            }                                                                             
            if (k < coordinate2 - 1) {                                                         
                matrice[(j * coordinate2 + k) * tailleGrille + (j * coordinate2 + k + 1)] += step2 * (-3.) / 4.;
            }                                                                             
            if (k < coordinate2 - 2) {                                                         
                matrice[(j * coordinate2 + k) * tailleGrille + (j * coordinate2 + k + 2)] += (step2) * (3.) / 20.;
            }
            if (k < coordinate2 - 3) {                                                         
                matrice[(j * coordinate2 + k) * tailleGrille + (j * coordinate2 + k + 3)] += (step2) * (-1.) / 60.;
            }
        }
    }
}

bool isAntiSymmetric(const vector<double>& matrice, long dimension) {
	int compteur = 0;

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            if (matrice[i * dimension + j] == -matrice[j * dimension + i]) {
                compteur += 1;
            }
        }
    }

	if (compteur == dimension*dimension) {
        return true;
    } else {
		return false;
	}
}

bool loadTau(double vecteur[], long dimension, string directory) {
	fstream entreeNAC;
    string nombre;
    entreeNAC.open(directory);
	if(entreeNAC){
    	for (int jter = 0; jter < dimension; jter++) {
        	getline(entreeNAC, nombre);
			//cout << "Nombre : " << nombre << endl;
        	vecteur[jter] = stod(nombre);
    	}
	} else {	
		return false;
		cerr << "Could not open: " << directory << endl;
	}
    entreeNAC.close();
}

void makeNac(double matrice[], double gradient[], double vecteur[], long dimension) {
	for (long i = 0; i < dimension; i++) {
		for (long j = 0; j < dimension; j++) {
			matrice[j * dimension + i] += gradient[j * dimension + i] * vecteur[j];
		}
	}
}

void makeNacAlt(double matrice[], double gradient[], double vecteur[], long dimension) {
	for (long i = 0; i < dimension; i++) {
		for (long j = 0; j < dimension; j++) {
			if(i > j)
				matrice[j * dimension + i] += gradient[j * dimension + i] * vecteur[j];
			else if(i < j && j > 0)
				matrice[j * dimension + i] += gradient[j * dimension + i] * vecteur[j-1];
		}
	}
}

void printMat(double matrice[], long dimension) {
	for (long i = 0; i < dimension; i++) {
		for (long j = 0; j < dimension; j++) {
			cout << matrice[i * dimension + j] << " ";
		}
		cout << endl;
	}	
}

/*I need to pass through the elements (elecState1 * gridDimension + i, elecState2 * gridDimension + j)
 *This leads to the iteration on the vector element
 *[(elecState1 * gridDimension + i) * fullDimension + (elecState2 * gridDimension + j)]
 */
void addNac(double matrice[], double matriceNac[], long elecState1, long elecState2, long fullDimension, long gridDimension) {
	for (long i = 0; i < gridDimension; i++) {
		for (long j = 0; j < gridDimension; j++) {
			matrice[(elecState1 * gridDimension + i) * fullDimension + 
					(elecState2 * gridDimension + j)] +=
					matriceNac[i * gridDimension + j];
			matrice[(elecState2 * gridDimension + j) * fullDimension + 
					(elecState1 * gridDimension + i)] +=
					matriceNac[i * gridDimension + j]; // transpose
		}
	}
}

void addNacAlt(double matrice[], double matriceNac[], long elecState1, long elecState2, long fullDimension, long gridDimension) {
	for (long i = 0; i < gridDimension; i++) {
		for (long j = 0; j < gridDimension; j++) {
			matrice[(elecState1 * gridDimension + i) * fullDimension + 
					(elecState2 * gridDimension + j)] +=
					matriceNac[i * gridDimension + j] - matriceNac[j * gridDimension + i];
			matrice[(elecState2 * gridDimension + j) * fullDimension + 
					(elecState1 * gridDimension + i)] +=
					matriceNac[i * gridDimension + j] - matriceNac[j * gridDimension + i]; // transpose
		}
	}
}

bool isSymmetric(double matrice[], long dimension) {
	long compteur = 0;
	for (long i = 0; i < dimension; i++) {
		for (long j = 0; j < dimension; j++) {
			if (matrice[i * dimension + j] == matrice[j * dimension + i]) {
				compteur++;
			}
		}
	}

	if (compteur == dimension * dimension) {
		return true;
	} else {
		return false;
	}
}

bool loadPotential(double vecteur[], long gridDimension, string directory, 
				string filenames[], string extension, int nFiles) {
	ifstream entree;
	string valeur;
	for (int i = 0; i < nFiles; i++) {
		entree.open(directory + filenames[i] + extension);
		if (entree) {
			for (long j = 0; j < gridDimension; j++) {
				getline(entree, valeur);
				vecteur[i * gridDimension + j] = stod(valeur);	
			}
		} else {
			cerr << "Could not open: " << directory + filenames[i] + extension << endl;
			return false;	
		}
		entree.close();
	}

	return true;
}

void addPotential(double matrice[], double vecteur[], long fullDimension) {
	for (long i = 0; i < fullDimension; i++) {
		matrice[i * fullDimension + i] += vecteur[i];
	}
}

bool writeSparseForm(double matrice[], long dimension, string filename) {
	ofstream sortie;
	sortie << setprecision(12);
	sortie.open(filename);
	if (sortie) {
		for (long i = 0; i < dimension; i++) {
			for (long j = 0; j < dimension; j++) {
				if (matrice[i * dimension + j] != 0.) {
					sortie << i << " " << j << " " << matrice[i * dimension + j] << endl;
				}
			}
		}
	} else {
		cerr << "Could not create the file: " << filename << endl;
		return false;
	}
	sortie.close();

	return true;
}

void printSparseMat(double matrice[], long dimension) {
	for (long i = 0; i < dimension; i++) {
		for (long j = 0; j < dimension; j++) {
			if (matrice[i * dimension + j] != 0.) {
				cout << matrice[i * dimension + j] << endl;	
			}
		}
	}
}
