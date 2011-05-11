#include "kelley.hh"
#include "delcher.hh"
#include <limits>

////////////////////////////////////////////////////////////////////////////////
// split
//
// Split on the character c, trying to match Python's split method
////////////////////////////////////////////////////////////////////////////////
vector<string> split(string s, char c)
{
  vector<string> splits;
  splits.push_back("");
  int split_num = 0;

  for(unsigned int i = 0; i < s.size(); i++) {
    if(s[i] == c) {
      split_num++;
      splits.push_back("");
    } else {
      splits[split_num] += s[i];
    }
  }
  
  return splits;
}


////////////////////////////////////////////////////////////////////////////////
// split
//
// Split on whitespace
////////////////////////////////////////////////////////////////////////////////
vector<string> split(string s) {
  vector<string> splits;
  unsigned int split_num = 0;
  bool last_space = true;

  for(unsigned int i = 0; i < s.size(); i++) {
    if(s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r') {
      if(!last_space)
	split_num++;
      last_space = true;
    } else {
      if(split_num == splits.size())
	splits.push_back("");
      splits[split_num] += s[i];
      last_space = false;
    }
  }

  return splits;
}


////////////////////////////////////////////////////////////////////////////////
// join
//
// Join a vector of strings separating them with the character c, trying to
// match Python's join method
////////////////////////////////////////////////////////////////////////////////
string join(vector<string> v, char c) {
  string s = v[0];
  for(unsigned int i = 1; i < v.size(); i++) {
    s += c;
    s += v[i];
  }
  return s;
}


////////////////////////////////////////////////////////////////////////////////
// upper
//
// Return an uppercase version of the string s.
////////////////////////////////////////////////////////////////////////////////
string upper(string s) {
     string u(s);
     for(unsigned int i = 0; i < s.size(); i++)
	  u[i] = toupper(s[i]);
     return u;
}

////////////////////////////////////////////////////////////////////////////////
// kernel_smooth
//
// Using Gaussian kernel
////////////////////////////////////////////////////////////////////////////////
void kernel_smooth(vector<float> & counts, float sigma)
{
     vector<double> dcounts(counts.size());
     for(unsigned int i = 0; i < counts.size(); i++)
	  dcounts[i] = (double)counts[i];

     kernel_smooth(dcounts, sigma);

     for(unsigned int i = 0; i < counts.size(); i++)
	  counts[i] = (float)dcounts[i];
}


////////////////////////////////////////////////////////////////////////////////
// kernel_smooth
//
// Using Gaussian kernel
////////////////////////////////////////////////////////////////////////////////
void kernel_smooth(vector<double> & counts, float sigma, unsigned int max_count)
{
     float sigma2 = pow(sigma, 2);
     int band = (int)(4*sigma);
     if(max_count == 0)
	  max_count = counts.size();

     vector<double> smooth_counts(counts);
     double num, den;

     // precompute Gaussians
     vector<double> gauss(band+1);
     for(unsigned int i = 0; i < gauss.size(); i++)
	  gauss[i] = exp(-pow((double)i, 2)/(2*sigma2));
     
     for(int l = 0; l < (int)max_count; l++) {
	  num = 0;
	  den = 0;

	  int lk_start = Max(0, l - band);
	  int lk_end = Min((int)max_count, l + band);
	  for(int lk = lk_start; lk < lk_end; lk++) {
	       num += counts[lk]*gauss[abs(lk - l)];
	       den += gauss[abs(lk - l)];
	  }
	  smooth_counts[l] = num/den;
     }

     for(unsigned int l = 0; l < max_count; l++)
	  counts[l] = smooth_counts[l];
}


////////////////////////////////////////////////////////////////////////////////
// log_add
//
// Add two numbers whose logarithms are given and return the logaritm.
////////////////////////////////////////////////////////////////////////////////
double log_add(const double l1, const double l2)
{
     if (l1 < -numeric_limits<double>::max() && l2 < -numeric_limits<double>::max())
	  return l1;
     else if(l1 > l2)
	  return l1 + log(1.0 + exp(l2 - l1));
     else
	  return l2 + log(1.0 + exp(l1 - l2));
}


////////////////////////////////////////////////////////////////////////////////
// coeff_log_add
//
// Add coeff times one number and (1-coeff) times the second when the numbers'
// logarithms are given and return the logaritm.
////////////////////////////////////////////////////////////////////////////////
double coeff_log_add(const double l1, const double l2, const double coeff)
{
     if (l1 < numeric_limits<double>::min() && l2 < numeric_limits<double>::min())
	  return l1;
     else if(l1 > l2)
	  return l1 + log(coeff + (1.0 - coeff)*exp(l2 - l1));
     else
	  return l2 + log(1.0 - coeff + coeff*exp(l1 - l2));
}


////////////////////////////////////////////////////////////////////////////////
// gamma_ml
//
// Compute maximum likelihood Gamma parameters
////////////////////////////////////////////////////////////////////////////////
void gamma_ml(double & k, double & theta, vector<double> & dist)
{
     unsigned int l;
     double N = 0;
     double sum_x = 0;
     double sum_lnx = 0;
     for(l = 1; l < dist.size(); l++) {
	  N += dist[l];
	  sum_x += l*dist[l];
	  sum_lnx += log((double)l)*dist[l];
     }
     double s = log(sum_x/N) - sum_lnx/N;
     k = (3.0 - s + sqrt((s - 3)*(s - 3) + 24*s)) / (12*s);
     theta = sum_x / (N*k);
}

////////////////////////////////////////////////////////////////////////////////
// geom_ml
//
// Compute maximum likelihood Geometric parameters, ignoring the first samples
// equal to 1, which tend to be unreliable here.
////////////////////////////////////////////////////////////////////////////////
void geom_ml(double & p, vector<double> & dist)
{
     double N = 0.0;
     double sum_x = 0.0;
     for(unsigned int l = 2; l < dist.size(); l++) {
	  N += dist[l];
	  sum_x += l*dist[l];
     }
     p = N / (sum_x + N);
}


////////////////////////////////////////////////////////////////////////////////
// normalize
//
// Normalizes probabilites
////////////////////////////////////////////////////////////////////////////////
void normalize(vector<double> & dist, unsigned int min_l)
{
     unsigned int l;
     double length_sum = 0;
     for(l = min_l; l < dist.size(); l++)
	  length_sum += dist[l];
     for(l = min_l; l < dist.size(); l++)
	  dist[l] /= length_sum;
}


////////////////////////////////////////////////////////////////////////////////
// log_normalize
//
// Normalizes probabilites in log-space
////////////////////////////////////////////////////////////////////////////////
void log_normalize(vector<double> & dist, unsigned int min_l)
{
     unsigned int l;
     double log_length_sum;
     double length_sum = 0;
     for(l = min_l; l < dist.size(); l++)
	  length_sum += exp(dist[l]);
     log_length_sum = log(length_sum);
     for(l = min_l; l < dist.size(); l++)
	  dist[l] -= log_length_sum;
}
