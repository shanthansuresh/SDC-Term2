#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"]; 
          double v = j[1]["speed"];

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          double steer_value;
          double throttle_value;

          /*
           * TODO: Incorporate latency.
           * There is a latency associated with the actuator commands being
           * fired and actually getting exectued. In this project we assume 
           * this latency to be 100ms. This latency can be taken into account
           * by considering the vehicle model(kinematic/dynamic). Note that
           * it is not possible to model this latency in the PID controller.
           * We use the kinematic model to predict the vehicle state at latency
           * time and feed that as the initial state to the model solver(lpopt)
           *
           * velocity: The model needs the velocity converted from mph to meters per second
           */
          /* 
           * Note that the steering command angle(delta) and vehicl model angle
           * conversions are reverse. Hence the negative sign here.
           */
          double latency = 0.1; // Set the expected latency here.
          double delta = j[1]["steering_angle"];
          delta *= -1;
          double a = j[1]["throttle"];
          double Lf = 2.67;
          v *= 0.44704;
          px += v*cos(psi)*latency;
          py += v*sin(psi)*latency;
          psi += v*delta*latency/Lf;
          v += a*latency;
 
          /* Cordinate Transformation.
           * For simplifying the problem, we convert the waypoints from
           * map cordinates to vehicle cordinates. Note that this is the
           * the reverse of what was done for particle filter project.
           * 
           * The transformation equations are :
           * x_waypoint_vehicle = (x_waypoint_map - x_vehicle_map) * cos(psi) + (y_waypoint_map - y_vehicle_map) * sin(psi)
           * y_waypoint_vehicle = -(x_waypoint_map - x_vehicle_map) * sin(psi) + (y_waypoint_map - y_vehicle_map) * cos(psi)
           *
           * where, (x_waypoint_map, y_waypoint_map) are waypoint cordinates in map framework,
           *        (x_vehicle_map, y_vehicle_map) are vehicle cordinates in vehicle framework,
           *        psi is the vehicle orientation w.r.t map framework,
           *        (x_waypoint_vehicle, y_waypoint_vehicle) are the waypoint cordinates in vehicle framework.
           *
           * Reference:
           * (a) https://discussions.udacity.com/t/do-you-need-to-transform-coordiantes/256483/11
           * (b) https://cdn-enterprise.discourse.org/udacity/uploads/default/original/4X/3/0/f/30f3d149c4365d9c395ed6103ecf993038b3d318.png
           */
          vector<double> ptsx_vehicle = vector<double>(ptsx.size());
          vector<double> ptsy_vehicle = vector<double>(ptsy.size());
          for (int j=0; j<ptsx.size(); j++) {
              double xdiff = ptsx[j] - px;
              double ydiff = ptsy[j] - py;
              ptsx_vehicle[j] = (xdiff* cos(psi)) + (ydiff*sin(psi));
              ptsy_vehicle[j] = -1*(xdiff* sin(psi)) + (ydiff*cos(psi));
          }

          /*
           * Fit a polynomial.
           * We fit a first order polynomial as described in the lesson.
           * The STL vectors need to be converted to Eigen::VectorXd format
           * as polyfit() expects in Eigen format.
           */
          Eigen::VectorXd ptsx_vehicle_ = Eigen::VectorXd(ptsx_vehicle.size());
          Eigen::VectorXd ptsy_vehicle_ = Eigen::VectorXd(ptsy_vehicle.size());

          for (int j=0; j<ptsx_vehicle.size(); j++) {
              ptsx_vehicle_[j] = ptsx_vehicle[j];
              ptsy_vehicle_[j] = ptsy_vehicle[j];
          }

          auto coeffs = polyfit(ptsx_vehicle_, ptsy_vehicle_, 1);

          /*
           * The cross track error is calculated by evaluating at polynomial at x, f(x)
           * and subtracting y.
           * Note that since we are working in vehicle cordinates(centered about the vehicle),
           * we evaluate the polynomial at x=0 and y=0)
           */
          double cte = polyeval(coeffs, 0) - 0;
   
          /*
           * Due to the sign starting at 0, the orientation error is -f'(x).
           * derivative of coeffs[0] + coeffs[1] * x -> coeffs[1]
           * Also the psi=0 as the x axis of vehicle framework is along the vehicle orientation.
           */
          double epsi = 0 - atan(coeffs[1]);

          Eigen::VectorXd state(6);
          state << 0, 0, 0, v, cte, epsi;

          vector<double> vars = mpc.Solve(state, coeffs);
          steer_value = -1 * vars[0]; //-1 to account for the right turn being positive angle.
          throttle_value = vars[1];
          steer_value /= deg2rad(25); //Normalize steer_value between [-1, 1]

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line
          for (int j=1; j<(vars.size()/2); j++) {
              mpc_x_vals.push_back(vars[2*j]);
              mpc_y_vals.push_back(vars[2*j+1]);
          }

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line
          int num_ref_points = 15;
          double ref_point_dist = 1.5;
          
          for (int j=0; j<num_ref_points; j++) {
              double x_ = j*ref_point_dist;
              next_x_vals.push_back(x_);
              next_y_vals.push_back(polyeval(coeffs, x_));
          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100)); //TODO: Revert this.
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
