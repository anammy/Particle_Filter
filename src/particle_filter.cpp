/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 10;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	for (int i=0; i < num_particles; i++){
		Particle pt;
		pt.id = i+1;
		pt.x = dist_x(gen);
		pt.y = dist_y(gen);
		pt.theta = dist_theta(gen);
		pt.weight = 1.0;
		particles.push_back(pt);
		weights.push_back(pt.weight);
	}
	
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	
	for (int i=0; i < num_particles; i++){
		if (fabs(yaw_rate) < 1e-3){
			particles[i].x += velocity*delta_t*cos(particles[i].theta) + dist_x(gen);
			particles[i].y += velocity*delta_t*sin(particles[i].theta) + dist_y(gen); 	
		}
		else{
			particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) + dist_x(gen); 
			particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)) + dist_y(gen);
			particles[i].theta += yaw_rate*delta_t + dist_theta(gen); 	
		}
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (unsigned int i=0; i < observations.size(); i++){
		//maximum finite value of type double
		double min_dist = numeric_limits<double>::max();

		for (unsigned int j=0; j < predicted.size(); j++){
			double curr_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if (curr_dist < min_dist){
				min_dist = curr_dist;
				observations[i].id = predicted[j].id;
			} 
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	for (int i=0; i < num_particles; i++){
		double xp = particles[i].x;
		double yp = particles[i].y;
		double theta = particles[i].theta;
		
		//Vector of landmarks within sensor range of the particle
		vector<LandmarkObs> predictions;
		LandmarkObs pred;
		for (unsigned int j=0; j < map_landmarks.landmark_list.size(); j++){
			pred.x = map_landmarks.landmark_list[j].x_f;
			pred.y = map_landmarks.landmark_list[j].y_f;
			pred.id = map_landmarks.landmark_list[j].id_i;
			
			if(dist(pred.x, pred.y, xp, yp) <= sensor_range){
				predictions.push_back(pred);
			}		
		}
		

		//Transform observations from vehicle coordinate system to map coordinate system
		vector<LandmarkObs> obs_map;
		LandmarkObs obs;
		for (unsigned int k=0; k < observations.size(); k++){
			obs.x = xp + cos(theta)*observations[k].x - sin(theta)*observations[k].y;
			obs.y = yp + sin(theta)*observations[k].x + cos(theta)*observations[k].y;
			obs.id = observations[k].id;
			obs_map.push_back(obs);
		}

		dataAssociation(predictions, obs_map);	

		//Initialize particles weight to 1;
		particles[i].weight=1;
		
		//Calculate observation weight and update particle weight
		for (unsigned int j=0; j < obs_map.size(); j++){
		obs.x = obs_map[j].x;
		obs.y = obs_map[j].y;
		obs.id = obs_map[j].id;

		//vector<LandmarkObs>::iterator itr = find_if(predicted.begin(), predicted.end(), [] (LandmarkObs obs){ return obs.id == *itr.id);

			for(unsigned int k=0; k < predictions.size(); k++){
				if (obs.id == predictions[k].id){
					pred.x = predictions[k].x;
					pred.y = predictions[k].y;
					pred.id = predictions[k].id;
				}
			}
		
		double p_obs = (0.5/(M_PI*std_landmark[0]*std_landmark[1]))*exp(-0.5*(pow((obs.x - pred.x)/std_landmark[0], 2) + pow((obs.y - pred.y)/std_landmark[1], 2))); 
		
		particles[i].weight *= p_obs;
		}
		weights[i] = particles[i].weight;
	}
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	vector<Particle> res_particles (num_particles);
	
	discrete_distribution<int> dist_pt(weights.begin(), weights.end());
	
	for (int i=0; i < num_particles; i++){
		res_particles[i] = particles[dist_pt(gen)];
	}

	particles = res_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
