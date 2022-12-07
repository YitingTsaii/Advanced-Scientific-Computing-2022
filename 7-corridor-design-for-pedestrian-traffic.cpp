#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>


double rij(std::vector<double> Xi, std::vector<double> Xj){ 
    // Xi = (Xi[0], Xi[1]) = (xi, yi)
    // Xj = (Xj[0], Xj[1]) = (xj, yj)
    double x = Xi[0] - Xj[0];
    double y = Xi[1] - Xj[1];
    double dist = sqrt( x*x + y*y);
    return dist;
}

// Fij function (1-dim, for x or y)
double Fij(std::vector<double> Xi, std::vector<double> Xj, double A, double B, int x_or_y){ 
    // x_or_y: 0 for x, 1 for y
    double zi = Xi[x_or_y];
    double zj = Xj[x_or_y];
    double r_ij = rij(Xi, Xj);
    double F_ij = A * exp(-r_ij/B) * (zi-zj)/r_ij;
    return F_ij;
}

// Fik function (1-dim, for x or y)
double Fik(std::vector<double> Xi, std::vector<double> Xk, double C, double D, int x_or_y){ 
    // x_or_y: 0 for x, 1 for y
    double zi = Xi[x_or_y];
    double zk = Xk[x_or_y];
    double r_ik = rij(Xi, Xk);
    double F_ik = C * exp(-r_ik/D) * (zi-zk)/r_ik;
    return F_ik;
}

// Fi function (1-dim, for x or y)
double Fi(int i, std::vector<double> v0, std::vector<double> v, double ped_loc[10][2], double obs_loc[5][2], double tau, double A, double B, double C, double D, int x_or_y){
    // v0 and v are of size 1*2, only pass in the specific v0, v for that particular i
    // i in 0,1,2,..,9
    double F_i = 0;
    F_i += (1/tau) * (v0[x_or_y] - v[x_or_y]);
    for(int j=0; j<10; j++){
        if(i == j){
            continue;
        }else{
            std::vector<double> Xi{ped_loc[i][0], ped_loc[i][1]};
            std::vector<double> Xj{ped_loc[j][0], ped_loc[j][1]};
            F_i += Fij(Xi, Xj, A, B, x_or_y);
        }
    }
    for(int k=0; k<5; k++){
        std::vector<double> Xi{ped_loc[i][0], ped_loc[i][1]};
        std::vector<double> Xk{obs_loc[k][0], obs_loc[k][1]};
        F_i += Fik(Xi, Xk, C, D, x_or_y);
    }
    return F_i;
}

// average displacement in the x direction
double get_x_displace(double ped_loc[10][2], double ped_loc0[10][2]){
    double x_displace = 0;
    for(int i=0; i<10; i++){
        x_displace += abs(ped_loc[i][0] - ped_loc0[i][0]);
    }
    x_displace = x_displace/10;
    return x_displace;
}

// cost function, input: the obstacle location (10 values), output: the average x displacement
double cost_function(double obs_loc[5][2], double D){
    // initialization
    double tau = 0.2;
    double A = 20;
    double B = 0.5;
    double C = 10;
    double dt = 0.05;
    double T = 50; // time duration

    double ped_loc0[10][2]; // pedestrian location at time 0
    double ped_loc[10][2]; // pedestrian location
    double ped_v[10][2];   // pedestrian velocity
    double ped_v0[10][2];  // pedestrian desired velocity
    double ped_F[10][2];   // force
    // double obs_loc[5][2];  // obstacle loction. This is now an input
    
    // check obstacles in dashed area
    for(int i = 0; i < 5; i++){
        if(obs_loc[i][0] < -5+D || obs_loc[i][0] > 5-D){ //x
            //std::cout << "x out of bound" << "\n";
            return 10000;
        }
        if(obs_loc[i][1] < -3+D || obs_loc[i][1] > 3-D){ //y
            //std::cout << "y out of bound" << "\n";
            return 10000;
        }
    }
    // check the obstacles do not overlap
    for(int i = 0; i < 4; i++){
        for(int j = i+1; j < 5; j++){
            //use rij function
            std::vector<double> Xi{obs_loc[i][0], obs_loc[i][1]};
            std::vector<double> Xj{obs_loc[j][0], obs_loc[j][1]};
            if(rij(Xi, Xj)<2*D){
                //std::cout << "i=" << i << ", j=" << j << ": overlap" << "\n";
                return 10000;
            }
        }
    }

    //  Q2(a)
    // the pedestrain initial location
    std::vector<double> x_initial{-28.5, -27.0, -25.5, -24.0, -22.5, 22.5, 24.0, 25.5, 27.0, 28.5};
    std::vector<double> y_initial{0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    for(int i = 0; i < 10; i++){
        ped_loc[i][0] = x_initial[i];
        ped_loc[i][1] = y_initial[i];
        ped_loc0[i][0] = x_initial[i]; // for calculating x displacement
        ped_loc0[i][1] = y_initial[i];
    }
    // the pedestrain initial velocity and the desired velocity
    std::vector<double> x_initial_v{1, 1, 1, 1, 1, -1, -1, -1, -1, -1};
    std::vector<double> y_initial_v{0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    for(int i = 0; i < 10; i++){
        ped_v[i][0] = x_initial_v[i];
        ped_v[i][1] = y_initial_v[i];
        ped_v0[i][0] = x_initial_v[i];
        ped_v0[i][1] = y_initial_v[i];
    }

    // initialize F (using x, v)
        for(int i = 0; i < 10; i++){
            std::vector<double> v0{ped_v0[i][0], ped_v0[i][1]};
            std::vector<double> v{ped_v[i][0], ped_v[i][1]};
            for(int x_or_y = 0; x_or_y < 2; x_or_y ++){
                ped_F[i][x_or_y] = Fi(i, v0, v, ped_loc, obs_loc, tau, A, B, C, D, x_or_y);
                //std::cout << ped_F[i][x_or_y] << " ";
            }
            //std::cout << "\n";
        }

    // the obstacle location. These are unknown for (b).
    // std::vector<double> x_obs{-4, -3, 0, 2, 5};
    // std::vector<double> y_obs{0, 3, -1, -2, 0};
    // for(int i = 0; i < 5; i++){
    //    obs_loc[i][0] = x_obs[i];
    //    obs_loc[i][1] = y_obs[i];
    //}

    // update F, v, x
    double current_t = 0;
    while(true){
        current_t += dt;
        //std::cout << "current_t = " << current_t << " ";
        if(current_t > T){
            break;
        }
        // update x
        for(int i = 0; i < 10; i++){
            for(int x_or_y = 0; x_or_y < 2; x_or_y ++){
                ped_loc[i][x_or_y] = ped_loc[i][x_or_y] + ped_v[i][x_or_y] * dt;
                //std::cout << ped_F[i][x_or_y] << " ";

                // put back to domain if pedestrians penetrates the wall
                if(x_or_y == 0){  //x
                    if(ped_loc[i][x_or_y] > 30){ ped_loc[i][x_or_y] = 30; }
                    if(ped_loc[i][x_or_y] < -30){ ped_loc[i][x_or_y] = -30; }
                }
                if(x_or_y == 1){  //y
                    if(ped_loc[i][x_or_y] > 3){ ped_loc[i][x_or_y] = 3; }
                    if(ped_loc[i][x_or_y] < -3){ ped_loc[i][x_or_y] = -3; }
                }
            }
            //std::cout << "\n";
        }
        // update v
        for(int i = 0; i < 10; i++){
            for(int x_or_y = 0; x_or_y < 2; x_or_y ++){
                ped_v[i][x_or_y] = ped_v[i][x_or_y] + ped_F[i][x_or_y] * dt;
                //std::cout << ped_F[i][x_or_y] << " ";
            }
            //std::cout << "\n";
        }
        // update F
        for(int i = 0; i < 10; i++){
            std::vector<double> v0{ped_v0[i][0], ped_v0[i][1]};
            std::vector<double> v{ped_v[i][0], ped_v[i][1]};
            for(int x_or_y = 0; x_or_y < 2; x_or_y ++){
                ped_F[i][x_or_y] = Fi(i, v0, v, ped_loc, obs_loc, tau, A, B, C, D, x_or_y);
                //std::cout << ped_F[i][x_or_y] << " ";
            }
            //std::cout << "\n";
        }
    }

    // print the average displacement in the x direction
    double x_displace = 0;
    x_displace = get_x_displace(ped_loc, ped_loc0);
    //std::cout << "x displace = " << x_displace << "\n";

    return -x_displace;
}

int main()
{
    double D = 0.6; // radius of the obstacle
    double obs_loc[5][2]; // x1 in python
    double x1[7][10]; //shape = NP*dimension (x1, y1, x2, y2, x3, y3, ...,x5, y5)
    double x2[7][10];

    int NP = 7;
    double F = 0.5;
    double CR = 0.1;
    int dimension = 10; // D in python
    int max_gen = 1000; // change here
    double trial[10]; // dimension=10
    double cost[7]; // since NP=7
    int count = 0;
    int a;
    int b;
    int c;
    int k;
    double score;

    // Initialize the particles using uniform distribution
    std::mt19937 gen(345);
    std::uniform_real_distribution<double> dist_x(-5.0+D, 5.0-D);
    std::uniform_real_distribution<double> dist_y(-3.0+D, 3.0-D);
    std::uniform_real_distribution<double> dist_uni(0, 1);
    for(int np=0; np<NP; np++){
        for(int i = 0; i < 5; i++){
            x1[np][2*i] = dist_x(gen);
            x1[np][2*i+1] = dist_y(gen);
            if(i>0){
                while(true){
                    int legal = 0;  // make sure the initial positions of the obstacles are legal
                    for(int j=i-1;j>=0;j--){
                        std::vector<double> Xi{x1[np][2*i], x1[np][2*i+1]};
                        std::vector<double> Xj{x1[np][2*j], x1[np][2*j+1]};
                        if(rij(Xi, Xj)>2*D){
                            legal += 1;
                        }
                    }
                    if(legal == i){
                        break;
                    }
                    x1[np][2*i] = dist_x(gen);
                    x1[np][2*i+1] = dist_y(gen);
                }
            }
        //std::cout << x1[np][2*i] << " " << x1[np][2*i+1] << " ";
        }
        //std::cout << "\n";
    }

    // initialize the cost
    for(int i=0; i<NP; i++){
        for(int m=0; m<5; m++){
                obs_loc[m][0] = x1[i][2*m];
                obs_loc[m][1] = x1[i][2*m+1];
        }
        cost[i] = cost_function(obs_loc, D);
    }
    std::cout << "iteration = " << count << ": ";
    std::cout << "cost = ";
    for(int i=0;i<NP;i++){
        std::cout << cost[i] << " ";
    }
    std::cout << "\n";

    

    while(true){
        count += 1;
        if(count==max_gen){
            break;
        }
        for(int i=0; i<NP; i++){
            while(true){
                a = rand() % NP;  // generate random integer from 0 to NP-1
                while(true){    // make sure b!= a
                    b = rand() % NP;
                    if(b != a){
                        break;
                    }
                }
                while(true){    // make sure b!= a
                    c = rand() % NP;
                    if(c != a && c!= b){
                        break;
                    }
                }
                if(i != a && i!= b && i!= c){
                    break;
                }
            }
            //std::cout << "a,b,c = " << a << " " << b << " " << c << "\n";
            k = rand() % dimension;
            
            for(int j=0; j<dimension; j++){
                if(j==k || dist_uni(gen)<CR){
                    trial[j] = x1[c][j] + F*(x1[a][j] - x1[b][j]);
                }else{
                    trial[j] = x1[i][j];
                }
            }
            // reshape the 10*1 matrix back to 5*2 matrix
            for(int m=0; m<5; m++){
                obs_loc[m][0] = trial[2*m];
                obs_loc[m][1] = trial[2*m+1];
            }
            score = cost_function(obs_loc, D);
            //std::cout << "score = " << score << " ";
            if(score <= cost[i]){
                for(int m=0; m<dimension; m++){
                    x2[i][m] = trial[m];
                }
                cost[i] = score;
            }else{
                for(int m=0; m<dimension; m++){
                    x2[i][m] = x1[i][m];
                }
            }
        }
        //assign x2 to x1
        for(int i=0; i<NP; i++){
            for(int j=0; j<dimension; j++){
                x1[i][j] = x2[i][j];
            }
        }
        std::cout << "iteration = " << count << ": ";
        std::cout << "cost = ";
        for(int i=0;i<NP;i++){
            std::cout << cost[i] << " ";
        }
        std::cout << "\n";
    }

    // choose the particle with the smallest cost function value
    double best_obs_loc[5][2];
    double best_cost = 0;
    for(int np=0; np<NP; np++){
        for(int i=0; i<5; i++){
            obs_loc[i][0] = x1[np][2*i];
            obs_loc[i][1] = x1[np][2*i+1];
        }
        if(np == 0){
            best_cost = cost_function(obs_loc, D);
            for(int i=0; i< 5; i++){
                best_obs_loc[i][0] = obs_loc[i][0];
                best_obs_loc[i][1] = obs_loc[i][1];
            }
        }else{
            if(cost_function(obs_loc, D) < best_cost){
                best_cost = cost_function(obs_loc, D);
                for(int i=0; i< 5; i++){
                    best_obs_loc[i][0] = obs_loc[i][0];
                    best_obs_loc[i][1] = obs_loc[i][1];
                }
            }
        }
    }
    
    std::cout << "best cost = " << best_cost << "\n";
    std::cout << "best obstacle positions: " << "\n";
    for(int i=0; i<5; i++){
        std::cout << "(" << best_obs_loc[i][0] << ", " << best_obs_loc[i][1] << ")\n";
    }

    return 0;
}
