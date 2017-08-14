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
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// DONE: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 100;
	particles.resize(num_particles);
	weights.resize(num_particles);

	std::random_device rd;
    std::mt19937 gen(rd());
	std::normal_distribution<double> d_x(x, std[0]);
	std::normal_distribution<double> d_y(y, std[1]);
	std::normal_distribution<double> d_theta(theta, std[2]);

	for (int i = 0; i < num_particles; i++) {
		particles[i].x = d_x(gen);
		particles[i].y = d_y(gen);
		particles[i].theta = d_theta(gen);
		particles[i].weight = 1.0;
		weights[i] = 0.;
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// DONE: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	std::random_device rd;
    std::mt19937 gen(rd());
	std::normal_distribution<double> d_x(0, std_pos[0]);
	std::normal_distribution<double> d_y(0, std_pos[1]);
	std::normal_distribution<double> d_theta(0, std_pos[2]);
	for (int i = 0; i < num_particles; i++) {
		//std::cout << "predict: pre " << particles[i].x << std::endl;

		if (fabs(yaw_rate) < 0.001) {
			particles[i].x += velocity*cos(particles[i].theta)*delta_t + d_x(gen);
			particles[i].y += velocity*sin(particles[i].theta)*delta_t + d_y(gen);
			particles[i].theta += d_theta(gen);
		} else {
			particles[i].x += velocity/yaw_rate 
				* (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) + d_x(gen);
			particles[i].y += velocity/yaw_rate 
				* (-cos(particles[i].theta + yaw_rate*delta_t) + cos(particles[i].theta)) + d_y(gen);
			particles[i].theta += yaw_rate*delta_t + d_theta(gen);
		}
		//std::cout << "predict: [" << i << "] x=" << particles[i].x << ", y=" << particles[i].y << std::endl;
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// I do not use it.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// DONE: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	double std_range = std_landmark[0];
	double std_bearing = std_landmark[1];
	double std_x = std_range;
	double std_y = std_range;
	double multi_gauss_const = 1./2./M_PI/std_x/std_y;

	double weights_sum = 0.0;
	for (int i = 0; i < num_particles; i++) {
		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;

		std::vector<Map::single_landmark_s> landmarks_in_range;
		//std::vector<double> std_ranges;

		double l_x, l_y;
		int o_id, l_id;
		double o_x, o_y;
		double d_x, d_y;
		double min_dist, the_dist;
		double range;

		struct obs2land_dist{
			int obs_index;
			double dist;
			double obs_x;
			double obs_y;
		};
		std::map<int,obs2land_dist> min_dist_for_landmark;


		for (int l = 0; l < map_landmarks.landmark_list.size(); l++) {
			l_id = map_landmarks.landmark_list[l].id_i;
			l_x = map_landmarks.landmark_list[l].x_f;
			l_y = map_landmarks.landmark_list[l].y_f;
			range = dist(particles[i].x, particles[i].y, l_x, l_y);
			if (range > sensor_range) {
				//std::cout << "update: [" << i << "] o=" << o << " l=" << l << " sensor_range out" << std::endl;
				continue;
			}
			landmarks_in_range.push_back(map_landmarks.landmark_list[l]);
			//std_ranges.push_back(range);
		}


		if (landmarks_in_range.size() == 0) {
			std::cout << "update: particle[" << i << "] has no close landmarks in range." << std::endl;
			particles[i].weight = 0.;
		} else {

			for (int o = 0; o < observations.size(); o++) {
				o_x = particles[i].x 
						+ observations[o].x * cos(particles[i].theta) - observations[o].y * sin(particles[i].theta);
				o_y = particles[i].y 
						+ observations[o].x * sin(particles[i].theta) + observations[o].y * cos(particles[i].theta);

				l_id = -1;
				min_dist = sensor_range;
				for (int l = 0; l < landmarks_in_range.size(); l++) {
					l_x = landmarks_in_range[l].x_f;
					l_y = landmarks_in_range[l].y_f;

					the_dist = dist(o_x, o_y, l_x, l_y);
					if (min_dist > the_dist) {
						min_dist = the_dist;
						l_id = landmarks_in_range[l].id_i;
						observations[o].id = l_id;

						obs2land_dist obs_dist;
						obs_dist.obs_index = o;
						obs_dist.dist = min_dist;
						obs_dist.obs_x = o_x;
						obs_dist.obs_y = o_y;
						
						pair<map<int,obs2land_dist>::iterator,bool> success = 
							min_dist_for_landmark.insert(pair<int,obs2land_dist>(l_id, obs_dist));
						if (false == success.second) {
							//std::cout << "update: particle[" << i << "] obs=" << o << " find the nearest observation." << std::endl;
							double old_dist = min_dist_for_landmark[l_id].dist;
							if (old_dist > min_dist) {
								min_dist_for_landmark[l_id] = obs_dist;
							}
						}
						
					}
				}
			}

			if (min_dist_for_landmark.size() == 0) {
				std::cout << "update: particle[" << i << "] has no assocations" << std::endl;
				particles[i].weight = 0.;
			} else {
				particles[i].weight = 1.;
				double prop;
				
				map<int,obs2land_dist>::iterator iter;
				for (iter = min_dist_for_landmark.begin(); iter != min_dist_for_landmark.end(); ++iter) {
					l_id = iter->first;
					Map::single_landmark_s landmark = map_landmarks.landmark_list[l_id - 1];
					obs2land_dist obs_dist = iter->second;

					d_x = obs_dist.obs_x - landmark.x_f;
					d_y = obs_dist.obs_y - landmark.y_f;
					prop = multi_gauss_const * exp(-0.5 * (pow(d_x/std_x, 2) + pow(d_y/std_y, 2)));
					particles[i].weight *= prop;

					associations.push_back(l_id);
					sense_x.push_back(obs_dist.obs_x);
					sense_y.push_back(obs_dist.obs_y);
				}
			}
		}
		weights[i] = particles[i].weight;
		weights_sum += weights[i];
		
		particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);
		//std::cout << "update: particle[" << i << "] has weights=" 
		//	<< weights[i] << ", associations=" << associations.size() << std::endl;
	}


	for (int i = 0; i < num_particles; i++) {
		if (weights_sum == 0.) {
			weights[i] = 1./double(num_particles);
		} else {
			weights[i] = weights[i]/weights_sum;
		}
	}

}

void ParticleFilter::resample() {
	// DONE: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::random_device rd;
    std::mt19937 gen(rd());
	std::vector<Particle> particles_new;
	std::discrete_distribution<> d(weights.begin(), weights.end());
	for (int i = 0; i < num_particles; i++) {
		Particle p = particles[d(gen)];
		particles_new.push_back(p);
	}
	particles = particles_new;
	/*
	for (int i = 0; i < num_particles; i++) {
		std::cout << "resampled: [" << i << "] weight=" << particles[i].weight 
				<< ", p(" << particles[i].x << "," << particles[i].y << ")" << std::endl;
	}
	*/

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
