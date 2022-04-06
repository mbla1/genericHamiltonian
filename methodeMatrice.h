#ifndef METHODEMATRICE_H_INCLUDED
#define METHODEMATRICE_H_INCLUDED
#include <string>
#include <vector>

void setZero(double vecteur[], long dimension);

double calculMoments(double coordinates1[], double coordinates2[], double masses[]);

void makeKineticOrder4(std::vector<double>& matrice, long elecStates, long coordinate1,
			   	long coordinate2, double step1, double step2, double mom1, double mom2, double mom12);

void makeGradient1(std::vector<double>& matrice, long coordinate1, long coordinate2, double step1, double reducedMass);

void makeGradient2(std::vector<double>& matrice, long coordinate1, long coordinate2, double step2);

void makeGradient6Order1(std::vector<double>& matrice, long coordinate1, long coordinate2, double step1, double reducedMass);

void makeGradient6Order2(std::vector<double>& matrice, long coordinate1, long coordinate2, double step2);

bool isAntiSymmetric(const std::vector<double>& matrice, long dimension);

bool loadTau(double vecteur[], long dimension, std::string directory);

void makeNac(std::vector<double>& matrice, const std::vector<double>& gradient, const Nac& vec_nac, int grid_size, int elec_size);

void makeNacAlt(double matrice[], double gradient[], double vecteur[], long dimension);

void printMat(double matrice[], long dimension);

void addNac(double matrice[], double matriceNac[], long elecState1, long elecState2, long fullDimension, long gridDimension);

void addNacAlt(double matrice[], double matriceNac[], long elecState1, long elecState2, long fullDimension, long gridDimension);

bool isSymmetric(const std::vector<double>& matrice, long dimension);

void addPotential(std::vector<double>& matrice, const Potential& pot, long fullDimension);

bool writeSparseForm(const std::vector<double>& matrice, long dimension, std::string filename);

void printSparseMat(double matrice[], long dimension);

#endif
