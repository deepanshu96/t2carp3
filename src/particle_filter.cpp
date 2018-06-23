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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	if(is_initialized==true)
		return;
	num_particles = 200;

	default_random_engine gen;
	double std_x, std_y, std_theta;
	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
 	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for(int i=0;i<num_particles;i++)
		{
			Particle temp;
			temp.id = i;
			temp.x = dist_x(gen);
			temp.y = dist_y(gen);
			temp.theta = dist_theta(gen);
			temp.weight = 1.0;
			particles.push_back(temp);
		}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	double std_x, std_y, std_theta;
	std_x = std_pos[0];
	std_y = std_pos[1];
	std_theta = std_pos[2];

	normal_distribution<double> dist_x(0, std_x);
 	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);

	for(int i=0;i<num_particles;i++)
		{
			if(fabs(yaw_rate)<0.00001)
				{
					particles[i].x += velocity*delta_t*cos(particles[i].theta);
					particles[i].y += velocity*delta_t*sin(particles[i].theta);
				}
			else
				{
					particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      				particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      				particles[i].theta += yaw_rate*delta_t;
				}
			particles[i].x += dist_x(gen);
			particles[i].y += dist_y(gen);
			particles[i].theta += dist_theta(gen);
		}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	LandmarkObs a,b;
	int x1 = observations.size();
	int x2 = predicted.size();
	//numeric_limits<double>::max()
	double temp,cur;
	int tempid;
	for(int i=0;i<x1;i++)
	{
		temp = numeric_limits<double>::max();
		a = observations[i];
		for(int j=0;j<x2;j++)
		{
			b = predicted[j];
			cur = dist(a.x, a.y, b.x, b.y);
			if(cur<temp)
			{
				temp = cur;
				tempid = b.id;
			}
		}
		observations[i].id = tempid;
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
	double x,y,theta,landx,landy,landid;
	for(int i=0;i<num_particles;i++)
	{
		x = particles[i].x;
		y = particles[i].y;
		theta = particles[i].theta;

		vector<LandmarkObs> nearland;
		for(int j=0;j<map_landmarks.landmark_list.size();j++)
		{
			landx = map_landmarks.landmark_list[j].x_f;
			landy = map_landmarks.landmark_list[j].y_f;
			landid = map_landmarks.landmark_list[j].id_i;

			if(fabs(landx - x)<=sensor_range && fabs(landy - y)<=sensor_range)
			{
				LandmarkObs temp;
				temp.x = landx;
				temp.y = landy;
				temp.id = landid;
				nearland.push_back(temp);
			}
		}

		vector<LandmarkObs> transformed;

		for(int k=0;k<observations.size();k++)
		{
			double x2 = cos(theta)*observations[k].x - sin(theta)*observations[k].y + x;
			double y2 = sin(theta)*observations[k].x + cos(theta)*observations[k].y + y;
			LandmarkObs temp;
			temp.x = x2;
			temp.y = y2;
			temp.id = observations[k].id;
			transformed.push_back(temp);
		}

		dataAssociation(nearland, transformed);

		particles[i].weight = 1;
		double x4,y4;
		for(int j=0;j<transformed.size();j++)
		{
			double x3 = transformed[j].x;
			double y3 = transformed[j].y;
			int id3 = transformed[j].id;
			for(int k=0;k<nearland.size();k++)
			{
				if(nearland[k].id == id3)
				{
					x4 = nearland[k].x;
					y4 = nearland[k].y;
				}
			}


			double gauss_norm = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);
			
			double exponent = (pow (x4-x3,2)/(2* pow (std_landmark[0],2) )) + (pow (y4-y3,2)/(2* pow (std_landmark[1],2)));

			double weight = gauss_norm*exp(-exponent);

			particles[i].weight *= weight;
		}
	 //for particles
	}
 //updateWeights
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	vector<double> wts;
	vector<Particle> particles2;
	for(int i=0;i<num_particles;i++)
	{
		wts.push_back(particles[i].weight);
	}
	double mxwt = (*( max_element(wts.begin(), wts.end()) ) );
	uniform_real_distribution<double> dist1(0.0, mxwt);
	uniform_int_distribution<int> dist2(0.0, num_particles);

	int index = dist2(gen);
	double beta = 0;
	for(int i=0;i<num_particles;i++)
	{
		beta += dist1(gen)*2;
		while(beta>wts[index])
		{
			beta -= wts[index];
			index = (index+1)%num_particles;
		}
		particles2.push_back(particles[index]);
	}

	particles = particles2;

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
