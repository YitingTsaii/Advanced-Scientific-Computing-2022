#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// f function
double f(int d, std::vector<double> x_sample){
    double total = 0;
    double ans = 0;
    for (int i=0; i<d; i++) {
        total += pow(x_sample[i], 2);
    }
    ans = pow(2, -d) * total;
    return ans;
}

// g function
double g(int d, std::vector<double> x_sample){
    for (int i=0; i<d; i++){
        if(x_sample[i] < -1 || x_sample[i] > 1){
            return 0;
        }
    }
    return f(d, x_sample);
}

// q_star function
double q_star(int d, std::vector<double> x_sample){
    for (int i=0; i<d; i++){
        if(x_sample[i] < -1 || x_sample[i] > 1){
            return 0;
        }
    }
    return 0.1 + f(d, x_sample);
}


int main()
{
    // variables that use across questions
    std::vector<int> M_vector{1, 10, 100, 1'000, 10'000, 100'000, 1'000'000};
    std::vector<int> d_vector{1, 2, 4, 8, 16};
    int M = 0;
    int d = 0;
    double I_MC_b = 0;
    double I_IS_c = 0;
    double I_IS_d = 0;
    std::mt19937 gen;
    std::random_device dev;
    gen.seed(dev());
    //std::default_random_engine gen;
    std::string name_part;
    std::string file_name;
    double result[100][7]; 
    
    // Q4(b)
    std::cout << "======== b =========" << "\n";
    for(int d_idx = 0; d_idx < d_vector.size(); d_idx ++){
        for(int M_idx = 0; M_idx < M_vector.size(); M_idx ++){
            M = M_vector[M_idx];
            d = d_vector[d_idx];
            std::cout << "M = " << M << ", ";
            std::cout << "d = " << d << "\n";

            std::vector<double> x_sample(d); // a vector storing the d-dimensional random samples
            std::uniform_real_distribution<double> dist(-1.0,1.0);
            
            for(int sample_idx = 0; sample_idx < 100; sample_idx ++){
                double current_f = 0;
                double sum_f = 0;

                for(int j = 0; j < M; j++){
                    // sample from uniform(0,1) to form a d-dimensional vector x_sample
                    for(int i = 0; i < x_sample.size(); i++){
                        x_sample[i] = dist(gen);
                        //std::cout << x_sample[i];
                    }
                    current_f = f(d, x_sample);
                    sum_f += current_f;
                }
                I_MC_b = sum_f * pow(2, d) / M;
                result[sample_idx][M_idx] = I_MC_b;
            }
            //std::cout << "I_MC_b = " << I_MC_b << "\n";
        }
        // save the result to a csv file
        name_part = "b_d";
        file_name = name_part + std::to_string(d) + ".csv";
        std::ofstream out(file_name);
        out << "M0,M1,M2,M3,M4,M5,M6," << std::endl;
        for (auto& row : result) {
            for (auto col : row){
                out << col <<',';
            }  
            out << '\n';
        }
        // print out the result matrix
        for (int i = 0; i < 100; i++){
            for (int j = 0; j < 7; j++){
                std::cout << result[i][j] << ' ';
            }
            std::cout << std::endl;
        }
    }
    

    // Q4(c)
    std::cout << "======== c =========" << "\n";
    //for(int d_idx = 4; d_idx < d_vector.size(); d_idx ++){  // only for d=16
    for(int d_idx = 0; d_idx < d_vector.size(); d_idx ++){
        for(int M_idx = 0; M_idx < M_vector.size(); M_idx ++){
            M = M_vector[M_idx];
            d = d_vector[d_idx];
            std::cout << "M = " << M << ", ";
            std::cout << "d = " << d << "\n";
            double mean = 0.0;
            double sigma = 0.5; 
            std::vector<double> x_sample(d); // a vector storing the d-dimensional random samples
            std::normal_distribution<double> distribution(mean, sigma);
            
            for(int sample_idx = 0; sample_idx < 100; sample_idx ++){
                double current_g_N = 0;
                double sum_g_N = 0;
                
                for(int j = 0; j < M; j++){
                    double normal_density = 1.0;
                    // sample from normal(0,0.5) to form a d-dimensional vector x_sample
                    for(int i = 0; i < x_sample.size(); i++){
                        x_sample[i] = distribution(gen);
                        //std::cout << x_sample[i] << " ";
                        normal_density *= ( 1/( pow(2*M_PI, 0.5) * sigma) ) *  exp(- pow(x_sample[i] - mean, 2)/(2*pow(sigma,2)));
                    }
                    current_g_N = g(d, x_sample)/normal_density;
                    sum_g_N += current_g_N;
                }
                I_IS_c = sum_g_N  / M;
                result[sample_idx][M_idx] = I_IS_c;
                //std::cout << "I_IS_c = " << I_IS_c << "\n";
            }
        }
        // save the result to a csv file
        name_part = "c_d";
        file_name = name_part + std::to_string(d) + ".csv";
        std::ofstream out(file_name);
        out << "M0,M1,M2,M3,M4,M5,M6," << std::endl;
        for (auto& row : result) {
            for (auto col : row){
                out << col <<',';
            }  
            out << '\n';
        }
        // print out the result matrix
        for (int i = 0; i < 100; i++){
            for (int j = 0; j < 7; j++){
                std::cout << result[i][j] << ' ';
            }
            std::cout << std::endl;
        }
    }
    
    
    // Q4(d)
    std::cout << "======== d =========" << "\n";
    for(int d_idx = 0; d_idx < d_vector.size(); d_idx ++){
        for(int M_idx = 0; M_idx < M_vector.size(); M_idx ++){
            M = M_vector[M_idx];
            d = d_vector[d_idx];
            std::cout << "M = " << M << ", ";
            std::cout << "d = " << d << "\n";
            // rejection sampling
            double gamma = pow(2, -d) * d;
            std::vector<double> x_sample(d); // a vector storing the d-dimensional random samples
            std::uniform_real_distribution<double> unif_dist(-1.0,1.0);
            std::uniform_real_distribution<double> unif_dist_u(0,gamma);
            
            for(int sample_idx = 0; sample_idx < 100; sample_idx ++){
            
                double u = 0;
                double current_w = 0;
                double sum_w = 0;
                double sum_w_g = 0;

                int j = 0;
                while (j < M){
                    // sample from uniform(0,1) to form a d-dimensional vector x_sample
                    for(int i = 0; i < x_sample.size(); i++){
                        x_sample[i] = unif_dist(gen);
                        //std::cout << x_sample[i];
                    }
                    u = unif_dist_u(gen); // u ~ uniform (0, gamma)
                    if ( u > q_star(d, x_sample) ){
                        // reject this x_sample
                        continue;
                    } else {
                        // accept this x_sample
                        current_w = 1/q_star(d, x_sample);
                        sum_w_g += current_w * g(d, x_sample);
                        sum_w += current_w;
                        j++;
                    }
                }
                I_IS_d = pow(2, d) * (sum_w_g  / sum_w); 
                result[sample_idx][M_idx] = I_IS_d;
                //std::cout << "I_IS_d = " << I_IS_d << "\n";
            }
        }
        // save the result to a csv file
        name_part = "d_d";
        file_name = name_part + std::to_string(d) + ".csv";
        std::ofstream out(file_name);
        out << "M0,M1,M2,M3,M4,M5,M6," << std::endl;
        for (auto& row : result) {
            for (auto col : row){
                out << col <<',';
            }  
            out << '\n';
        }
        // print out the result matrix
        for (int i = 0; i < 100; i++){
            for (int j = 0; j < 7; j++){
                std::cout << result[i][j] << ' ';
            }
            std::cout << std::endl;
        }
    }
    
    
    return 0;
}




