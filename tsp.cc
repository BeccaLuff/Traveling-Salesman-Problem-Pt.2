#include "cities.hh"
#include <iostream>
#include <fstream>
#include <ostream>
int main(int argc,char *argv[]){ 
  if(argc<2){ //Make sure filename was provided
 	  std::cerr<<"Need input filename\n";
	  return -1;
  }
  
  std::ifstream fin(std::string(argv[1]), std::ios::in);
 
  if(!fin.is_open()){ //Make sure file can be opened
    std::cerr<<"Can't open input file "<<argv[1]<<"\n";
    return -2; 
  } 
  Cities cities;
  fin>>cities;
  auto cities_size = cities.get_list().size();

  Cities::permutation_t permutation;
  for(unsigned int i=0; i<cities_size; i++){permutation.push_back(i);}
  auto shortest = cities.total_path_distance(permutation);
   
  for(int i=0; i<1000000; i++){ 
    auto new_permutation = random_permutation(cities_size);
    auto new_distance = cities.total_path_distance(new_permutation);
	
    if(new_distance<shortest){
      shortest = new_distance;
      cities = cities.reorder(new_permutation);
      std::cout<<i<<"\t"<<new_distance<<"\n";
	  }
  }
       
  std::ofstream ordering(std::string("shortest.tsv"), std::ios::out);
  ordering<<cities;
  return 0;}