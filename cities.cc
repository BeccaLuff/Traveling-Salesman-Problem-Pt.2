#include "cities.hh"
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <assert.h>
#include <cmath>
#include <random>
#include <algorithm>

//Implement operator<< (std::istream&, Cities&) to allow for the outputting of a Cities object to a stream
std::ostream& operator<<(std::ostream& out, Cities &cities)
{
 Cities::city_list_t cityList = cities.get_list();	
 for(Cities::coord_t i: cityList){
     out << i.first << "\t" << i.second << "\n";
 }
 return out;
}

//Implement operator >> (std::istream&, Cities&) to allow for the reading  / filling of a Cities object from a stream.
std::istream& operator>>(std::istream& is, Cities &cities)
{
 int x;
 int y;
 Cities::city_list_t cityList; 
 while(is >> x >> y){
     Cities::coord_t loc = {x, y};
     cities.add_coord(loc);
 }
 return is;
}

// For a given permutation of the cities in this object,
// compute how long (distance) it would take to traverse all the cities in the
// order of the permutation, and then returning to the first city.
// The distance between any two cities is computed as the Euclidean
// distance on a plane between their coordinates.
double Cities::total_path_distance(const permutation_t& ordering) const{
  double distance = 0;
  int x1, y1, x2, y2;
  long unsigned int i = 0;
  while( i < ordering.size()){
	  x1 = city_list_[ordering[i]].first;
	  y1 = city_list_[ordering[i]].second;
	  x2 = city_list_[ordering[i+1]].first;
	  y2 = city_list_[ordering[i+1]].second;
	  distance += std::hypot(x2 - x1, y2 - y1);
	  i++;
  }
  //add the distance between the first and last cities
  int firstx = city_list_[ordering[0]].first;
  int lastx= city_list_[ordering[ordering.size() - 1]].first;
  int firsty = city_list_[ordering[0]].second;
  int lasty= city_list_[ordering[ordering.size() - 1]].second;
  distance += std::hypot(lastx - firstx, lasty - firsty);
  return distance;

}

// Given a permutation, return a new Cities object where the order of the
// cities reflects the original order of this class after reordering with
// the given ordering. So for example, the ordering { 1, 0 } simply swaps
// the first two elements (coordinates) in the new Cities object.
Cities Cities::reorder(const permutation_t& ordering) const{
  //change city list so that is reflects the permutation 
  Cities new_city_order;
  for(long unsigned int i = 0; i < ordering.size(); i++){
	  new_city_order.add_coord(city_list_[ordering[i]]);
  }
  return new_city_order;
}

//generates and returns a new permutation of all the numbers from 0 to len -1
Cities::permutation_t random_permutation(unsigned len){
   Cities::permutation_t randomized;
   //obtain a random number 
   std::random_device rd;
   //seed generator
   std::mt19937 eng(rd());
   //std::default_random_engine generator;
   //defining range
   std::uniform_int_distribution<int> distribution(0, len - 1);
   
   //while the randomized permutation is not full
   while(randomized.size() < len){
	//get random integer   
	int tmp = distribution(eng);
	//if the random integer is not already in the randomized permuation
	//add to randomized permutation
        if(std::find(randomized.begin(), randomized.end(), tmp) 
			== randomized.end()){
		randomized.push_back(tmp);
	}
   }
   return randomized;

}
