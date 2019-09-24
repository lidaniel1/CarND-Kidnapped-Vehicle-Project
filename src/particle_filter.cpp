/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
 
  num_particles = 1000;  // Set the number of particles
  
  // define normal distribution with mean & sigma
  std::default_random_engine gen;
  std::normal_distribution<double> dtr_x(x, std[0]);
  std::normal_distribution<double> dtr_y(y, std[1]);
  std::normal_distribution<double> dtr_theta(theta, std[2]);
  
  Particle item;
  
  for (int i=0;i<num_particles;++i){
    item.x= dtr_x(gen);
    item.y = dtr_y(gen);
    item.theta = dtr_theta(gen);
    item.weight = 1;
    particles.push_back(item);
  }
  is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
    */
  
  // define random Gaussian noise
  std::default_random_engine gen;
  std::normal_distribution<double> dtr_x(0.0, std_pos[0]);
  std::normal_distribution<double> dtr_y(0.0, std_pos[1]);
  std::normal_distribution<double> dtr_theta(0.0, std_pos[2]);
  
  for (int i=0; i<num_particles;++i){
    //std::cout<<"i"<<i<<";x"<<particles[i].x<<std::endl;
    //std::cout<<"y"<<particles[i].y<<std::endl;
    //if(i==0){
    //std::cout<<"yaw_rate"<<yaw_rate<<std::endl;
    //}
    
    // prediction with two different sets of equation depending on whether yaw_rate is zero
    if(yaw_rate==0){
        particles[i].x += velocity*sin(particles[i].theta)*delta_t + dtr_x(gen);
    	particles[i].y += velocity*cos(particles[i].theta)*delta_t + dtr_y(gen);
    }
    else{
    particles[i].x += velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) + dtr_x(gen);
    particles[i].y += velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)) + dtr_y(gen);
    }
    particles[i].theta += yaw_rate*delta_t + dtr_theta(gen);
    //std::cout<<"x"<<particles[i].x<<std::endl;
    //std::cout<<"y"<<particles[i].y<<std::endl;
    //std::cout<<"theta"<<particles[i].theta<<std::endl;
      }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  // not used
  double min_dist;
  double calc_dist;
  double x;
  double y;
  int vlen1;
  int vlen2;
  int id;
  vlen1 = predicted.size();
  vlen2 = observations.size();

  for (int i=0; i<vlen1; ++i){
    min_dist = dist(predicted[i].x,predicted[i].y,observations[0].x,observations[0].y);
    id = observations[0].id;
    for (int j=1; j<vlen2; ++j){
      calc_dist = dist(predicted[i].x,predicted[i].y,observations[j].x,observations[j].y);
      if (calc_dist < min_dist){
        min_dist = calc_dist;
        id = observations[j].id;
        x = observations[j].x;
        y = observations[j].y;
      }
     predicted[i].id = id; 
     predicted[i].x = x;
     predicted[i].y = y;
    }
  }	
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  double x;
  double y;
  double theta;
  double calc_dist;
  double min_dist;
  double prob;
  double sum_weight;
  vector<LandmarkObs> tobserv; // transformed observation
  vector<LandmarkObs> aobserv; //associated observation
  LandmarkObs item;
  double gauss_norm;
  gauss_norm = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);
  //std::cout<<"gauss norm"<<gauss_norm<<std::endl;
  //std::cout<<"map landmark size"<<map_landmarks.landmark_list.size()<<std::endl;
  //std::cout<<"observation size"<<observations.size()<<std::endl;  
  
  //loop through each particle for transformation, association and update weight
  for(int i=0;i<num_particles;++i){
    /* transformation for each observation */
    tobserv.clear();
    aobserv.clear();
    x = particles[i].x;
    y = particles[i].y;
    theta = particles[i].theta;
    prob = 1;
   
    for (vector<LandmarkObs>::const_iterator it=observations.begin();it!=observations.end();++it){
           /* translate x y to map coordinate */
           item.x = it->x*cos(theta) - it->y*sin(theta) + x;
           item.y = it->x*sin(theta) + it->y*cos(theta) + y;
           item.id = it->id;
           tobserv.push_back(item); 
          if (i==0){
            std::cout<<"particle position"<<x<<";"<<y<<";"<<theta<<std::endl;
            std::cout<<"transformed observation"<<item.x<<";"<<item.y<<";"<<item.id<<std::endl;
          }
    }
    //std::cout<<"observation size"<<tobserv.size()<<std::endl;
    //associate observation
    //associate each observation to a landmark (closest observation to a landmark)
    for (unsigned int j=0; j<tobserv.size(); ++j){
      min_dist = dist(tobserv[j].x,tobserv[j].y,map_landmarks.landmark_list[0].x_f,
                      map_landmarks.landmark_list[0].y_f);
      if(i==0){
      std::cout<<"min dist"<<min_dist<<std::endl;
      }
      item.x = map_landmarks.landmark_list[0].x_f;
      item.y = map_landmarks.landmark_list[0].y_f;
      item.id =map_landmarks.landmark_list[0].id_i;
      for (unsigned int k=1; k<map_landmarks.landmark_list.size();++k){
        calc_dist = dist(tobserv[j].x,tobserv[j].y,map_landmarks.landmark_list[k].x_f,
                      map_landmarks.landmark_list[k].y_f);
        if (calc_dist<min_dist){
          min_dist = calc_dist;
          item.x =  map_landmarks.landmark_list[k].x_f;
          item.y =  map_landmarks.landmark_list[k].y_f;
          item.id = map_landmarks.landmark_list[k].id_i;
        }
      }
      aobserv.push_back(item);
      if(i==0){
      	std::cout<<"map land mark"<<item.x<<";"<<item.y<<std::endl;
        std::cout<<"observation land mark"<<tobserv[j].x<<";"<<tobserv[j].y<<std::endl;
      }
      prob*=  gauss_norm* exp(-((pow(item.x - tobserv[j].x,2)/(2*pow(std_landmark[0],2)))  + (pow(item.y - tobserv[j].y,2)/(2*pow(std_landmark[1],2))))); 
      //std::cout<<"prob"<<prob<<std::endl;
   }
   particles[i].weight = prob;
   // std::cout<<"ith weight"<<i<<";"<<prob<<std::endl;
  }
  // normalize weight
  sum_weight = 0;
  for (unsigned int i=0;i<particles.size();++i){
    sum_weight += particles[i].weight;
  }
  for (unsigned int i=0;i<particles.size();++i){
    particles[i].weight/=sum_weight;
  }
}                            

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  std::vector<double> w;
  vector<Particle> pr;
  int index;
  
  // resample using discrete_distribution function
  for (vector<Particle>::iterator it=particles.begin();it!=particles.end();++it){
    w.push_back(it->weight);
  }
  
  pr = particles;
  
  std::default_random_engine gen;
  std::discrete_distribution<int> dtr(w.begin(),w.end());
  
  for (int i=0; i<num_particles; ++i){
    index = dtr(gen);
  	particles[i] = pr[index];  
  }
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}