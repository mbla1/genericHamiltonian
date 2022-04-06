#include <iostream>
#include <vector>
#include <fstream>
#include "lectureFichiers.h"

using namespace std;

vector<double> load_file(string path) {
	fstream entree;
	vector<double> dipoles;
	string temp;

	entree.open(path);
	if(entree) {
		getline(entree, temp); // Skip first line
		for(double d; entree >> d; ){
			dipoles.push_back(d);
		}
		return dipoles;
	} else {
		cerr << "Could not open the file: " << path << endl;
	}
}

vector<double> transpose(vector<double> &vec, int row_init, int col_init){
	vector<double> temp;

	for(int j = 0; j < col_init; j++)
		for(int i = 0; i < row_init; i++)
			temp.push_back(vec[i * col_init + j]);

	return temp;
}

vector<double> remove_radius(vector<double> &vec, int row_init, int col_init){
	vector<double> temp;

	for(int j = 1; j < col_init; j++) // Remove the first column
		for(int i = 0; i < row_init; i++)
			temp.push_back(vec[j * row_init + i]); // Order

	return temp;
}

vector<double> prepare_file(string path, int row_init, int col_init){
	vector<double> data;

	data = load_file(path);
	data = transpose(data, row_init, col_init);
	data = remove_radius(data, row_init, col_init);

	return data;
}

