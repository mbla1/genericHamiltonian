#include <vector>

std::vector<double> load_file(std::string path);

std::vector<double> transpose(std::vector<double> &vec, int row_init, int col_init);

std::vector<double> remove_radius(std::vector<double> &vec, int row_init, int col_init);

std::vector<double> prepare_file(std::string path, int row_init, int col_init);
