#ifndef PID_H
#define PID_H

class PID {
public:
  /*
  * Errors
  */
  double p_error;
  double i_error;
  double d_error;

  /*
  * Coefficients
  */ 
  double Kp;
  double Ki;
  double Kd;

  /*
  * Total error
  */
  double total_error; 

  /*
   * Best error so far
   */
  double best_error;

  /*
  * Number of steps seen by the pid after twiddle
  */
  int num_steps;

  /*
   * delta values for twiddle 
   */
  double dp[3];

  /*
   * Index for tweaking p, i, d gains.
   */
  int idx;
  int go_up[3];

  /*
  * Constructor
  */
  PID();

  /*
  * Destructor.
  */
  virtual ~PID();

  /*
  * Initialize PID.
  */
  void Init(double Kp, double Ki, double Kd);

  /*
  * Update the PID error variables given cross track error.
  */
  void UpdateError(double cte);

  /*
  * Calculate the total PID error.
  */
  double TotalError();

  /*
   * Get pid output value.
   */
  double getOutput();

  /*
   * Get the step count.
   */
  int get_num_steps();

  /*
   * Set the step count to a value.
   */
  void set_num_steps(int value);

  /*
   * twiddle 
   */
  void twiddle();
};

#endif /* PID_H */
