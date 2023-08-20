#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>

const int N_conf = 10000;
const int N_dump = 5000; 
const int N_skip = 10;
const int M = (N_conf-N_dump)/N_skip;  
const int T = 50;
const int N_t = 500;   
const double a = T/N_t;
const double delta = 0.05;

double jaknf(const std::vector<double>& u) {
  double sum = 0;
  for (auto x : u) sum += x;
  return (sum - u.size()*u[0]) / (u.size() - 1);
}

double modelFunc(double t, const std::vector<double>& p) {
  return p[0] * exp(-p[1] * t);  
}

// void nonlinearFit(const std::vector<double>& x, const std::vector<double>& y,
//                   const std::vector<double>& yerr, std::function<double(double,const std::vector<double>&)> modelfunc) {

//   // prior guess 
//   std::vector<double> prior = {1, 1};
  
//   // fit model
//   // get params and errors 
  
//   // print results
// }

double weight(const std::vector<double>& u, int t) {
  double u_b_t = u[(t-1+N_t)%N_t];
  double u_f_t = u[(t+1)%N_t];
  double u_t = u[t];
   
  double S = (u_f_t-u_t)*(u_f_t-u_t)/(2.0*a) + 
            a*(u_f_t*u_f_t + u_t*u_t)/4.0 +
            (u_t-u_b_t)*(u_t-u_b_t)/(2.0*a) +  
            a*(u_t*u_t + u_b_t*u_b_t)/4.0;
           
  return exp(-S);
}


int main() {

  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  std::vector<std::vector<double>> Un(N_conf, std::vector<double>(N_t, 0));
  std::vector<double> old_u(N_t, 0);
  std::vector<double> new_u(N_t, 0);

  for (int n = 0; n < N_conf; n++) {
    std::vector<int> t_iter(N_t);
    std::iota(t_iter.begin(), t_iter.end(), 0);
    std::random_shuffle(t_iter.begin(), t_iter.end());

    for (int t : t_iter) {
      for (int i = 0; i < N_t; i++) {
        new_u[i] = old_u[i] + delta * dist(generator);
      }

      if (std::uniform_real_distribution<double>(0,1)(generator) < 
          weight(new_u,t) / weight(old_u,t)) {
        old_u[t] = new_u[t];
      }
      else {
        new_u[t] = old_u[t];
      }
    }

    Un[n] = old_u;
  }

  std::vector<std::vector<double>> Um(M, std::vector<double>(N_t, 0));
  for (int i = 0; i < M; i++) {
    Um[i] = Un[N_dump + i*N_skip];
  }

  // 计算Um后100个数组的平均值
  const int n = 100; 
  double tmp=0;
  std::vector<double> mean_last(n, 0);
for (int j = M-n; j < M; j++) {
  tmp=0
  for (int i = 0; i < N_t; i++) {
      tmp += Um[j][i];
    }
    mean_last[j] = tmp/N_t; 
  }

  // 打印输出
  for (auto x : mean_last) {
    std::cout << x << " ";
  }
  std::cout << std::endl;
  
  return 0;
}