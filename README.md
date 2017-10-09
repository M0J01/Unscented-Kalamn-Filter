
# Unscented Kalman Filter 

## Description

This project uses an Unscented Kalman Filter to track and predict the path of a bicycle around a moving car.
The Filter takes in noisy Radar and Lidar data, and uses this to guess position, velocity, and turn rate of the target bike.
In this project, we perform the update step at the same time as the prediction step, however this could be changed to
a constant prediction at certain time itnervals, and perform an update step once triggered by an incoming measurement.

The project relies on a simulator to provide input data, and returns information to the simulator.
The simulator can be found [here](https://github.com/udacity/self-driving-car-sim/releases)

This project was completed as part of the UDACITY Self-Driving Car Engineer Nanodegree Program. For more information, or to enroll today, please check [here.](https://www.udacity.com/course/self-driving-car-engineer-nanodegree--nd013) 


## Process Flow

The filter performs the following major steps. These steps are called from the ukf.cpp file. It is important to note that
the Main.cpp file performs administrative tasks such as communicating with the simulator, however, it does not perform
calculations.

### Initialization
Initialization is called when the first measurement from Lidar or Radar is recieved. This function, along with internal state variables, provide the state vector x_ and the Covariance matrix P_ which represents our tracking targets current or future position. 

Initialization is performed in the Process Measurement function, howeve rit is performed first, and only performed once, so it has been addresses separately in this document.

### Process Measurement
This function performs basic oversight on our process. It checks initializations status, calculates the time between the 
last measurement, calls our prediction step, and determines which type of sensor update to perform.


### Prediction
The Prediction function calculates Sigma Points in a 7 dimensional vector. The sigma points are a generated through the 
use of our 5 state dimension values, 2 noise values, a noise modified Process covariance matrix, and a "spreading constant". 
These Sigma points allow us to then perform time step predictions on where our target will be after a given amount of 
time has passed. 

### Update (Lidar or Radar)
This step combines our predicted sigma points from the prior step, and our actual measurment from a sensor (Lidar or Radar).

The update functions convert our predicted sigma points into the state space of the measuring sensor (2D for Lidar, 3D for Radar),
and use these in combination with our current state, covariance, and sensor noise variables, to update the belief about the
 trajectory of our target.
 
 ### Calculate RMSE
 
 The Root Mean Square Error of the Kalman filter belief vs the ground truth is calculated in the tools.cpp file. It is 
 called by Main.cpp after the prdecition and update steps are called.
 
 
 ## Increasing Performance
 
 Aside from performing all math calculations correctly, and implemneting a proper model, the initialization of the 
 Covariance Matrix P_, state matrix x_, and noise variables std_a_ and std_yawdd_ have the largest affect on your 
 filter accuracy. Playing with these variables will yield significantly different results when running.
 
 **Note :** Careful about setting the values of std_a_. std_yawdd_ or P_ too high. It is possible that the code will hang
 as the covarance matrix and state vector can explode rather quickly. A possible soultion would be implementing a hard
 reset on state and covariance if this occurs.





## Dependencies / Installation
This repository includes two files that can be used to set up and intall [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows you can use either Docker, VMware, or even [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) to install uWebSocketIO. Please see [this concept in the classroom](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/16cf4a78-4fc7-49e1-8621-3450ca938b77) for the required version and installation scripts.

Once the install for uWebSocketIO is complete, the main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./UnscentedKF

Tips for setting up your environment can be found [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/23d376c7-0195-4276-bdf0-e02f1f3c665d)

Note that the programs that need to be written to accomplish the project are src/ukf.cpp, src/ukf.h, tools.cpp, and tools.h

The program main.cpp has already been filled out, but feel free to modify it.

Here is the main protcol that main.cpp uses for uWebSocketIO in communicating with the simulator.


INPUT: values provided by the simulator to the c++ program

["sensor_measurement"] => the measurment that the simulator observed (either lidar or radar)


OUTPUT: values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x
["estimate_y"] <= kalman filter estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]

---


## Other Important Dependencies
* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF` Previous versions use i/o from text files.  The current state uses i/o
from the simulator.



## Generating Additional Data

This is optional!

If you'd like to generate your own radar and lidar data, see the
[utilities repo](https://github.com/udacity/CarND-Mercedes-SF-Utilities) for
Matlab scripts that can generate additional data.
