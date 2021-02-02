#pragma once

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng *RandomNumberGenerator; /* Random number generator. */
unsigned RandomSeed; /* Random seed. */

unsigned GetRandomSeed()
{
  FILE *random_file;
  unsigned seed;
  int error;

  seed = 0;
  random_file = fopen("/dev/urandom", "rb");
  if (!random_file)
    return 0;
  do
  {
    error = fread(&seed, 1, sizeof(unsigned int), random_file);
    if (error != sizeof(unsigned int))
      return 0;
  }
  while (!seed);
  fclose(random_file);
  return seed;
}

/** InitialiseRandom.

    - Initialises the random number generator. **/

int InitialiseRandom(const gsl_rng_type *random_generator)
{
  RandomSeed = GetRandomSeed();
  if (!RandomSeed)
    return EXIT_FAILURE;
  RandomNumberGenerator = gsl_rng_alloc(random_generator);
  if (!RandomNumberGenerator)
    return EXIT_FAILURE;
  gsl_rng_set(RandomNumberGenerator, RandomSeed);
  return EXIT_SUCCESS;
}

/** FinaliseRandom.

    - Finalises the random number generator. **/

int FinaliseRandom()
{
  /* Release random number generator memory. */

  if (RandomNumberGenerator)
  {
    gsl_rng_free(RandomNumberGenerator);
    return EXIT_SUCCESS; // EXIT_SUCCESS = 0
  }
  return EXIT_FAILURE; // EXIT_FAILURE = 1
}

/////////////
