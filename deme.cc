/*
 * Declarations for Deme class to evolve a genetic algorithm for the
 * travelling-salesperson problem.  A deme is a population of individuals.
 */

#include "chromosome.hh"
#include "deme.hh"
#include <chrono>

///////////////////////////////////////////////////////////////
// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].
Deme::Deme(const Cities* cities_ptr, unsigned pop_size, double mut_rate)
{
	mut_rate_ = mut_rate;
	for(unsigned i = 0; i < pop_size; i++)
	{
		pop_.push_back(new Chromosome (cities_ptr));
	}

}

///////////////////////////////////////////////////////////////
// Clean up as necessary
Deme::~Deme()
{
  // Add your implementation here
}

///////////////////////////////////////////////////////////////
// Evolve a single generation of new chromosomes, as follows:
// We select pop_size/2 pairs of chromosomes (using the select() method below).
// Each chromosome in the pair can be randomly selected for mutation, with
// probability mut_rate, in which case it calls the chromosome mutate() method.
// Then, the pair is recombined once (using the recombine() method) to generate
// a new pair of chromosomes, which are stored in the Deme.
// After we've generated pop_size new chromosomes, we delete all the old ones.
void Deme::compute_next_generation()
{
	child_pop_.clear();
	while (child_pop_.size() < pop_.size())
	{
		auto parent1 = select_parent();
		auto parent2 = select_parent();

		if (random_double() < mut_rate_)
		{
			parent1->mutate();
		}
		if (random_double() < mut_rate_)
		{
			parent2->mutate();
		}
		auto child_pair = parent1->recombine(parent2);
		child_pop_.push_back(child_pair.first);
		child_pop_.push_back(child_pair.second);
	}

	pop_.clear();
	for (auto i:child_pop_)
	{
		pop_.push_back(i);
	}
}

///////////////////////////////////////////////////////////////
// Return a copy of the chromosome with the highest fitness.
const Chromosome* Deme::get_best() const
{
	auto best = pop_[0];
	for (auto i:pop_)
	{
		if (i->get_fitness()>best->get_fitness())
		{
			best = i;
		}
	}
	return best;
}

///////////////////////////////////////////////////////////////
// Randomly select a chromosome in the population based on fitness and
// return a pointer to that chromosome.
Chromosome* Deme::select_parent()
{
	double rand_gen = random_double();

	double target_fitness = rand_gen * total_fitness();
	double current_total = 0;

	for(auto i:pop_)
	{
		current_total += i->get_fitness();
		if (current_total > target_fitness)
		{
			return i;
		}
	}

	// if the sum never exceeds target, return last chromosome in the vector
	return pop_[pop_.size()-1];

}

///////////////////////////////////////////////////////////////
// Returns the sum of fitness of all Chromosomes
double Deme::total_fitness() const
{
	double total = 0.;
	for (auto i:pop_)
	{
		total+=i->get_fitness();
	}
	return total;
}
//////////////////////////////////////////////////////////////
// Returns a random double between 0 and 1

double Deme::random_double() const
{

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator_(seed); // A random number generator for the various methods
	//generate a random double between 0 and 1
	std::uniform_real_distribution<double> distribution(0.0,1.0);
//	generator_.seed(seed);
	return distribution(generator_);
}
