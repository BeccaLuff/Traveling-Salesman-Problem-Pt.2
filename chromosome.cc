 /*
 * Implementation for Chromosome class
 */

#include <algorithm>
#include <cassert>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <chrono>

#include "chromosome.hh"

//////////////////////////////////////////////////////////////////////////////
// Generate a completely random permutation from a list of cities
Chromosome::Chromosome(const Cities* cities_ptr)
  : cities_ptr_(cities_ptr),
    order_(random_permutation(cities_ptr->size())),
    score_(get_fitness())
{
  assert(is_valid());
 // std::cout << "Creation success";
}

//////////////////////////////////////////////////////////////////////////////
// Clean up as necessary
Chromosome::~Chromosome()
{
  assert(is_valid());

}


/////////////////////////////////////////////////////////////////////////////
// Perform a single mutation on this chromosome
void
Chromosome::mutate()
{
  
  int i, j;

  i = rand() % order_.size();
  j = rand() % order_.size();

  /*  To make sure i!=j */
  while(i == j){
    j = rand() % order_.size();
  }

  std::swap(order_[i], order_[j]);

  assert(is_valid());
}

//////////////////////////////////////////////////////////////////////////////
// Return a pair of offsprings by recombining with another chromosome
// Note: this method allocates memory for the new offsprings

std::pair<Chromosome*, Chromosome*>
Chromosome::recombine(const Chromosome* other)
{
	assert(is_valid());
	assert(other->is_valid());

	int half_order_size = order_.size()/2;
	std::uniform_int_distribution<> distribution(0,half_order_size);
	
	generator_.seed(std::chrono::system_clock::now().time_since_epoch().count());
	auto begin = distribution(generator_);

//	std::cout<<"\nThe generated range is:"<<begin<<"\t size is :"<<order_.size()<<"\n";

	
	std::pair<Chromosome*, Chromosome*> offspring_pair;
	offspring_pair = 
		std::make_pair(
				create_crossover_child(this, other, begin, begin+half_order_size-1),
				create_crossover_child(other, this, begin, begin+half_order_size-1));
	return offspring_pair;
}

//////////////////////////////////////////////////////////////////////////////
// For an ordered set of parents, return a child using the ordered crossover.
// The child will have the same values as p1 in the range [b,e),
// and all the other values in the same order as in p2.
Chromosome*
Chromosome::create_crossover_child(const Chromosome* p1, const Chromosome* p2,
                                   unsigned b, unsigned e) const
{
  Chromosome* child = p1->clone();

  // We iterate over both parents separately, copying from parent1 if the
  // value is within [b,e) and from parent2 otherwise
  unsigned i = 0, j = 0;

  for ( ; i < p1->order_.size() && j < p2->order_.size(); ++i) {
    if (i >= b and i < e) {
      child->order_[i] = p1->order_[i];
    }
    else { // Increment j as long as its value is in the [b,e) range of p1
      while (p1->is_in_range(p2->order_[j], b, e)) {
        ++j;
      }
      assert(j < p2->order_.size());
      child->order_[i] = p2->order_[j];
      j++;
    }
  }

  assert(child->is_valid());
  return child;
}

// Return a positive fitness value, with higher numbers representing
// fitter solutions (shorter total-city traversal path).
double
Chromosome::get_fitness() const
{
  double score;

  /*Assuming a total_distance above 1500 miles for a TSP renders negative profit,
    giving each unit of distance a cost of -4, our score is our maximum subtracting
    our cost. Hence, shorter distances have a higher score */

  score = 6000 - 4 * (this->calculate_total_distance());
  return score;  
}

// A chromsome is valid if it has no repeated values in its permutation,
// as well as no indices above the range (length) of the chromosome.
// We implement this check with a sort, which is a bit inefficient, but simple
bool
Chromosome::is_valid() const
{
    Cities::permutation_t temp = order_;
    std::sort(temp.begin(), temp.end());

    for(unsigned long int i=0; i < order_.size(); i++){
      if(temp[i]!=i)  {return false;}
    }

    return true;
    
}

// Find whether a certain value appears in a given range of the chromosome.
// Returns true if value is within the specified the range specified
// [begin, end) and false otherwise.

bool
Chromosome::is_in_range(unsigned value, unsigned begin, unsigned end) const
{
  for(unsigned int i = begin; i < end; i++){
    if(order_[i] == value)  {return true;}
  }

  return false;
}
