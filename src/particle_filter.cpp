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
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	random_device rd;
    mt19937 gen(rd());
	
	normal_distribution<double> dist_x(x,std[0]);
	normal_distribution<double> dist_y(y,std[1]);
	normal_distribution<double> dist_theta(theta,std[2]);
	
	num_particles = 10;
	
	Particle p;
	for (int i=0; i<num_particles; i++) {
		
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);

		p.weight = 1.0;
		particles.push_back(p);
		weights.push_back(p.weight);
		}
	
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	double x_del, y_del, theta_del;
	
	random_device rd;
    mt19937 gen(rd());
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0,std_pos[2]);
	
	for (int i=0; i<particles.size(); i++) {
		
		if (fabs(yaw_rate) > 0.00001) {
			x_del = velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t)\
									-sin(particles[i].theta));
									
			y_del = velocity/yaw_rate * (cos(particles[i].theta)\
									-cos(particles[i].theta + yaw_rate*delta_t));
		}
		else {
			x_del = velocity * delta_t * cos(particles[i].theta);
			y_del = velocity * delta_t * sin(particles[i].theta);
		}

		theta_del = yaw_rate * delta_t;

		particles[i].x += x_del + dist_x(gen);
		particles[i].y += y_del + dist_y(gen);
		particles[i].theta += theta_del + dist_theta(gen);	
		}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	double ob_x, ob_y, landmark_x, landmark_y;
	double min_dist, distance;
	
	for (int i=0; i<observations.size(); i++) {
				
		ob_x = observations[i].x;
		ob_y = observations[i].y;

		landmark_x = predicted[0].x;
		landmark_y = predicted[0].y;
					
		min_dist = dist(ob_x, ob_y, landmark_x, landmark_y);
		observations[i].id = predicted[0].id;
		
		//find nearest neighbor 
		for (int j=1; j<predicted.size(); j++) {
			
			landmark_x = predicted[j].x;
			landmark_y = predicted[j].y;
			
			distance = dist(ob_x, ob_y, landmark_x, landmark_y);
			
			if (distance < min_dist) {
				observations[i].id = predicted[j].id;
				min_dist = distance;
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
	
	const double std_x = std_landmark[0];
	const double std_y = std_landmark[1];
	
	vector<LandmarkObs> predicted;
	vector<LandmarkObs> observations_t;

	for (int k=0; k<particles.size(); k++) {
		
		const double p_x = particles[k].x;
		const double p_y = particles[k].y;
		const double p_theta = particles[k].theta;

		//only includes landmarks that are within the sensor range
		predicted.clear();
		for (int i=0; i<map_landmarks.landmark_list.size(); i++) {
			
			double landmark_x = map_landmarks.landmark_list[i].x_f;
			double landmark_y = map_landmarks.landmark_list[i].y_f;
			
			double distance = dist(p_x, p_y, landmark_x, landmark_y);
			
			if (distance <= sensor_range) {
				LandmarkObs lm_in_range;
				lm_in_range.x = map_landmarks.landmark_list[i].x_f;
				lm_in_range.y = map_landmarks.landmark_list[i].y_f;
				lm_in_range.id = map_landmarks.landmark_list[i].id_i;
				predicted.push_back(lm_in_range);
			}
		}
		
		//convert vehicle coordinates to map coordinates
		observations_t.clear();
		for (int i=0; i<observations.size(); i++) {		
			
			LandmarkObs ob_t;
			ob_t.id = -1;
			ob_t.x = cos(p_theta)* observations[i].x - sin(p_theta) * observations[i].y + p_x;
			ob_t.y = sin(p_theta)* observations[i].x + cos(p_theta) * observations[i].y + p_y;
			observations_t.push_back(ob_t);
		}
		
		//
		dataAssociation(predicted, observations_t);
		
		//
		particles[k].associations.clear();
		particles[k].sense_x.clear();
		particles[k].sense_y.clear();
		particles[k].weight = 1;
		for (int i=0; i<observations_t.size(); i++) {
			
			int association = observations_t[i].id;
			double ob_x = observations_t[i].x;
			double ob_y = observations_t[i].y;
			double diff_x = ob_x - map_landmarks.landmark_list[association-1].x_f;
			double diff_y = ob_y - map_landmarks.landmark_list[association-1].y_f;

			double norm = 1 / (2 * M_PI * std_x * std_y);
			double exponent = pow(diff_x, 2) / (2 * pow(std_x, 2))\
								  + pow(diff_y, 2) / (2 * pow(std_y, 2));
			particles[k].weight *= norm * exp(-exponent);
			
			particles[k].associations.push_back(association);
			particles[k].sense_x.push_back(ob_x);
			particles[k].sense_y.push_back(ob_y);
			}

			weights.push_back(particles[k].weight);		
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	vector<Particle> new_particles;
	random_device rd;
    mt19937 gen(rd());
	discrete_distribution<> resample_dist(weights.begin(), weights.end());
	
	for (int i=0; i<num_particles; i++) {
		new_particles.push_back(particles[resample_dist(gen)]);
	}
	
	particles = new_particles;
	weights.clear();
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
