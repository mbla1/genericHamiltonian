#ifndef METHODEMATRICE_H_INCLUDED
#define METHODEMATRICE_H_INCLUDED
#include <string>

void setZero(double vecteur[], long dimension);

double calculMoments(double coordinates1[], double coordinates2[], double masses[]);

void makeKineticOrder2(double matrice[], long elecStates, long coordinate1,
			   	long coordinate2, double step1, double step2, double mom1, double mom2, double mom12);

void makeKineticOrder4(double matrice[], long elecStates, long coordinate1,
			   	long coordinate2, double step1, double step2, double mom1, double mom2, double mom12);

void makeGradient1(double matrice[], long coordinate1, long coordinate2, double step1, double reducedMass);

void makeGradient2(double matrice[], long coordinate1, long coordinate2, double step2);

void makeGradient6Order1(double matrice[], long coordinate1, long coordinate2, double step1, double reducedMass);

void makeGradient6Order2(double matrice[], long coordinate1, long coordinate2, double step2);

bool isAntiSymmetric(double matrice[], long dimension);

bool loadTau(double vecteur[], long dimension, std::string directory);

void makeNac(double matrice[], double gradient[], double vecteur[], long dimension);

void makeNacAlt(double matrice[], double gradient[], double vecteur[], long dimension);

void printMat(double matrice[], long dimension);

void addNac(double matrice[], double matriceNac[], long elecState1, long elecState2, long fullDimension, long gridDimension);

void addNacAlt(double matrice[], double matriceNac[], long elecState1, long elecState2, long fullDimension, long gridDimension);

bool isSymmetric(double matrice[], long dimension);

bool loadPotential(double vecteur[], long gridDimension, std::string directory, 
				std::string filenames[], std::string extension, int nFiles);

void addPotential(double matrice[], double vecteur[], long fullDimension);

bool writeSparseForm(double matrice[], long dimension, std::string filename);

void printSparseMat(double matrice[], long dimension);

#endif
