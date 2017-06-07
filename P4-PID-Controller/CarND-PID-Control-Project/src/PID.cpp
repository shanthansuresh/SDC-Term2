#include "PID.h"
#include <ctime>
#include <iostream>
#include <math.h>

using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {}

PID::~PID() {}

void PID::Init(double kp, double ki, double kd) {
    Kp = kp;
    Ki = ki;
    Kd = kd;

    num_steps = 0;

    total_error = 0;
    best_error = 1e8;

    dp[0] = 0.01;
    dp[1] = 0.0001;
    dp[2] = 0.01;

    go_up[0] = go_up[1] = go_up[2] = 0;

    idx = 0;
}

double cte_prev = 0;  //TODO: Ok to set this to 0?
double cte_sum = 0; 
 
void PID::UpdateError(double cte) {
    /* 
     * P error.
     */
    p_error = cte;

    /*
     * D error
     */
    d_error = (cte - cte_prev); 
    cte_prev = cte;
   
    /*
     * I error
     */
    cte_sum += cte; 
    i_error = cte_sum; 

    /* Update number of steps. */
    num_steps += 1;

    /*
     * Update total error
     */
    total_error += pow(cte, 2);
}

double PID::TotalError() {
    return total_error/(num_steps);
}

int PID::get_num_steps() {
    return num_steps;
}

double PID::getOutput() {
    double pid_out;
 
    pid_out = (-1*Kp*p_error) + (-1*Kd*d_error) + (-1*Ki*i_error);
    if (pid_out > 1) {
        pid_out = 1;
    } else if (pid_out < -1) {
        pid_out = -1;
    }
    
    return pid_out; 
}

/*
 * This function is called to tweak the gain parameters.
 * In the first step, update the 
 * gain parameters. Run the twiddle algorithm
 * after 300 iterations.
 * At the end of the twiddle algorithm the num_steps
 * is reset to 0 and the process repeats.
 */
void PID::twiddle() {
    double p[3] = {Kp, Ki, Kd};

    if ((best_error > 1e7) && (num_steps == 300)) {
        best_error = TotalError();
        return;
    }

    if (go_up[idx] == 0) {
        p[idx] += dp[idx];
        go_up[idx] = -1; 
    }

    if (num_steps < 300) {
        return;
    }

    double err = TotalError();

    std::cout << "idx: " << idx << " err " << err << " best err " << best_error << std::endl;

    if (err < best_error) {
        best_error = err;
        dp[idx] *= 1.1;
        
        idx = (idx+1)%3;
        go_up[idx] = 0;

        //std::cout << __LINE__ << std::endl;
    } else {
        if (go_up[idx] == -1) {
            p[idx] -= 2*dp[idx];
            go_up[idx] = +1; 
        
            //std::cout << __LINE__ << std::endl;
        } else if (go_up[idx] == 1) {
            p[idx] += dp[idx];
            dp[idx] *= 0.9; 
        
            //total_error = 0;
            idx = (idx+1)%3;
            go_up[idx] = 0;

            //std::cout << __LINE__ << std::endl;
        } else {
            std::cout << "Error, shouldn't land here " << std::cout;
        }
    }
    
    total_error = 0;
    num_steps = 0;
   
    Kp = p[0];
    Ki = p[1];
    Kd = p[2];

    std::cout << "Gain: " << Kp << " " << Ki << " " << Kd << std::endl;
}
