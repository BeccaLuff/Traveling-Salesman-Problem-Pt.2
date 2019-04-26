/*
 * Declarations for Deme class to evolve a genetic algorithm for the
 * travelling-salesperson problem.  A deme is a population of individuals.
 */

#include "chromosome.hh"
#include "deme.hh"
#include <iostream>

// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].
Deme::Deme(const Cities* cities_ptr, unsigned pop_size, double mut_rate)
{
 pop_size_ = pop_size;
 mut_rate_ = mut_rate;
 for(unsigned i = 0; i < pop_size; i++){
         pop_.push_back(new Chromosome(cities_ptr));
 }

}

// Clean up as necessary
Deme::~Deme()
{

}

// Evolve a single generation of new chromosomes, as follows:
// We select pop_size/2 pairs of chromosomes (using the select() method below).
// Each chromosome in the pair can be randomly selected for mutation, with
// probability mut_rate, in which case it calls the chromosome mutate() method.
// Then, the pair is recombined once (using the recombine() method) to generate
// a new pair of chromosomes, which are stored in the Deme.
// After we've generated pop_size new chromosomes, we delete all the old ones.
void Deme::compute_next_generation()
{
   next_gen_.clear();
   for(unsigned i; i<pop_size_/2; i++){
      //randomly select parents
      auto mother = select_parent();
      auto father = select_parent();
      //generate random number
      std::uniform_real_distribution<double> dis(0.0, 1.0);
      double rand_num_fem = dis(generator_);
      double rand_num_mal = dis(generator_);
      //if rand num smaller than mutation probibility, then mutate
      if(rand_num_fem  < mut_rate_)     mother->mutate();
      if(rand_num_mal  < mut_rate_)     father->mutate();
      //recombine the two parents
      std::pair<Chromosome*, Chromosome*> children = mother->recombine(father);
      next_gen_.push_back(children.first);
      next_gen_.push_back(children.second);
   }
   pop_.clear();
   for(auto i: next_gen_)       pop_.push_back(i);

}

// Return a copy of the chromosome with the highest fitness.
const Chromosome* Deme::get_best() const
{ Chromosome* best_chrome = pop_[0];
  for(auto chrom: pop_){
          if(chrom->get_fitness() > best_chrome->get_fitness())
                  best_chrome = chrom;
  }
  return best_chrome;
}

// Randomly select a chromosome in the population based on fitness and
// return a pointer to that chromosome.
Chromosome* Deme::select_parent()
{
  //calculate sum of population fitness
  double s = 0;
  for(unsigned i = 0; i < pop_size_; i++){
          s += pop_[i]->get_fitness();
  }

  //generate a random number between 0 and s
  std::uniform_int_distribution<int> dis(0, s);
  int rand_num = dis(generator_);

  //starting from the top of the population keep adding the fitnesses
  //to the partial sum P until P < S
  //the individual for which P exceeds S is the chosen individual
  double partial_sum = 0;
  for(auto i : pop_){
          partial_sum += i ->get_fitness();
          if(rand_num < partial_sum)    return i;
  }
  return nullptr;
}

