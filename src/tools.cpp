#include <iostream>
#include "tools.h"
#include <time.h>


using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

	//time_t time1;
	//time(&time1);

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	if(estimations.size() != ground_truth.size()
		 || estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
		//std::cout << " Estimation : \n" << estimations[i] << "\n\nGround Truth : \n" << ground_truth[i] << std::endl << std::endl;

	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//time_t time2;
	//time(&time2);
	//std::cout <<"Calculating RMSE took  " << (time2 - time1) << " Seconds \n";
	//return the result
	return rmse;


}