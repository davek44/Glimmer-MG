#include <vector>
#include <string>
using namespace::std;

vector<string> split(string s, char c);
vector<string> split(string s);

string join(vector<string> v, char c);

string upper(string s);

void kernel_smooth(vector<float> & counts, float sigma);
void kernel_smooth(vector<double> & counts, float sigma, unsigned int max_count=0);

double log_add(const double l1, const double l2);
double coeff_log_add(const double l1, const double l2, const double coeff);

void gamma_ml(double & k, double & theta, vector<double> & dist);
void geom_ml(double & p, vector<double> & dist);

void normalize(vector<double> & dist, unsigned int min_l);
void log_normalize(vector<double> & dist, unsigned int min_l);
