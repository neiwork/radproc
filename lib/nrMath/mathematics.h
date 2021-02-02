#pragma once

/*** mathematics.h:

     Mathematical and statistical functions. ***/

/** Avoid including this file twice. **/

#ifndef MATHEMATICS_H
#define MATHEMATICS_H

/** Libraries. **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/** Constants. **/

#define MATH_PI 3.1415926535897932385 /* Pi. */
#define MATH_2PI 6.2831853071795864770 /* 2 pi. */
#define MATH_PI_2 1.57079632679489661925 /* Pi / 2. */

/* Coordinates / axes. */

#define MATH_X 0 /* Castesian. */
#define MATH_Y 1
#define MATH_Z 2
#define MATH_RADIUS 0 /* Spherical. */
#define MATH_POLAR 1
#define MATH_AZIMUTH 2
#define MATH_MU 1 /* Pseudospherical (mu = cos(polar) replacing polar). */

/** Variables. **/

/* Random numbers. */

gsl_rng *RandomNumberGenerator; /* Random number generator. */
unsigned RandomSeed; /* Random seed. */
unsigned SortComponent = 0; /* Component used for sorting vectors. */

/** Function prototypes. **/

int SortAscending(const void *, const void *);
int SortIndexesAscending(const void *, const void *);
double ScalarProduct(double *, double *, unsigned);
double Norm(double *, unsigned);
double NormSquared(double *, unsigned);
double AngleCosine(double *, double *);

void SphericalToCartesian(double *, double *);
void PseudoSphericalToCartesian(double *, double *);
void CartesianToSpherical(double *, double *);
void CartesianToPseudoSpherical(double *, double *);
void Rotate(double *, double *, unsigned, double);
unsigned GetRandomSeed();
int InitialiseRandom(const gsl_rng_type *);
int FinaliseRandom();
double MetropolisHastings(double *, double *, double *, double *, unsigned, unsigned, unsigned, unsigned, double *, double *, double *,
                          void (*)(double *, double *, double *, double *, double *, unsigned), double (*)(double *, double *, double *, double *, double *, unsigned),
                          double (*)(double *));
double ParallelTempering(double *, double *, double *, double *, unsigned, unsigned, unsigned, unsigned, double *, double *, double *,
                         void (*)(double *, double *, double *, double *, double *, unsigned), double (*)(double *, double *, double *, double *, double *, unsigned),
                         double (*)(double *), unsigned, double *, double);
void TanWeights(double *, double *, double *, double *, unsigned, unsigned, unsigned, double *, double *, double *,
                double (*)(double *, double *, double *, double *, double *, unsigned), unsigned);
double TanNormalisation(double *, unsigned, double *);
double WeightAndSum(double *, double *, unsigned);
void UniformTransition(double *, double *, double *, double *, double *, unsigned);
double UniformTransitionPDF(double *, double *, double *, double *, double *, unsigned);
double PowerLawDeviate(double, double, double);
void IsotropicDirection(double, double *);
void StandardGaussianDeviates(double *);
void GaussianDeviates(double *, double, double);
complex ComplexGaussianDeviates(complex, complex);

/** SortAscending.

    - Function designed to be an argument of qsort().
    - Sorts the elements of a vector in ascending order. The vector elements are reorganised. **/

int SortAscending(const void *x, const void *y)
{
  double difference;

  difference = *((double *)x) - *((double *)y);
  if (difference > 0.0)
    return 1;
  else if (difference < 0.0)
    return -1;
  else
    return 0;
}

/** SortIndexesAscending.

    - Function designed to be an argument of qsort().
    - Sorts the indexes of a vector, in such a way that its elements are arranged in ascending order. The vector remains unchanged.
    - 'SortComponent' indicates the component used to sort, in the case that the vector elements have more than one dimension. **/

int SortIndexesAscending(const void *x, const void *y)
{
  double difference;

  difference = *(*(double **)x + SortComponent) - *(*(double **)y + SortComponent);
  if (difference > 0.0)
    return 1;
  else if (difference < 0.0)
    return -1;
  else
    return 0;
}

/** ScalarProduct.

    - Computes the scalar product between two vectors. **/

double ScalarProduct(double *vector_1, double *vector_2, unsigned vector_size)
{
  unsigned i;
  double product;

  /* Compute components product and sum. */

  product = 0.0;
  for (i = 0; i < vector_size; i++)
    product += vector_1[i] * vector_2[i];
  return product;
}

/** Norm.

    - Computes the norm of a vector. **/

double Norm(double *vector, unsigned vector_size)
{
  unsigned i;
  double *component_sq, max_component_sq, norm;

  /* Generate an auxiliary vector of component squares, in ascending order. */

  component_sq = malloc(vector_size * sizeof(double));
  for (i = 0; i < vector_size; i++)
    component_sq[i] = vector[i] * vector[i];
  qsort(component_sq, vector_size, sizeof(double), SortAscending);

  /* Divide component squares by the largest one to avoid overflow. */

  max_component_sq = component_sq[vector_size-1];
  if (max_component_sq <= 0.0) /* All components are null. */
    return 0.0;
  vector_size--;
  for (i = 0; i < vector_size; i++)
    component_sq[i] /= max_component_sq;

  /* Sum component squares in ascending order, to avoid loss of significance. */

  norm = 0.0;
  for (i = 0; i < vector_size; i++)
    norm += component_sq[i];
  norm += 1.0;

  /* Renormalise. */

  norm = sqrt(max_component_sq * norm);
  free(component_sq);
  return norm;
}

/** NormSquared.

    - Computes the square of the norm of a vector. **/

double NormSquared(double *vector, unsigned vector_size)
{
  unsigned i;
  double *component_sq, max_component_sq, norm;

  /* Generate an auxiliary vector of component squares, in ascending order. */

  component_sq = malloc(vector_size * sizeof(double));
  for (i = 0; i < vector_size; i++)
    component_sq[i] = vector[i] * vector[i];
  qsort(component_sq, vector_size, sizeof(double), SortAscending);

  /* Divide component squares by the largest one to avoid overflow. */

  max_component_sq = component_sq[vector_size-1];
  if (max_component_sq <= 0.0) /* All components are null. */
    return 0.0;
  vector_size--;
  for (i = 0; i < vector_size; i++)
    component_sq[i] /= max_component_sq;

  /* Sum component squares in ascending order, to avoid loss of significance. */

  norm = 0.0;
  for (i = 0; i < vector_size; i++)
    norm += component_sq[i];
  norm += 1.0;

  /* Renormalise. */

  norm *= max_component_sq;
  free(component_sq);
  return norm;
}

/** AngleCosine.

    - Computes the cosine of the angle between two 3D vectors. **/

double AngleCosine(double *vector_1, double *vector_2)
{
  double norm_1, norm_2, cosine;

  norm_1 = Norm(vector_1, 3);
  norm_2 = Norm(vector_2, 3);
  if (norm_1 <= 0.0 || norm_2 <= 0.0) /* At least one vector is null. */
    return 0.0;
  cosine = ScalarProduct(vector_1, vector_2, 3) / (norm_1 * norm_2);
  if (cosine > 1.0) /* Correct eventual round-off errors. */
    cosine = 1.0;
  else if (cosine < -1.0)
    cosine = -1.0;
  return cosine;
}

/** SphericalToCartesian.

    - Converts a 3D vector from spherical (radius, polar angle, azimuth) to Cartesian form.
    - Initial and final vectors may be the same. **/

void SphericalToCartesian(double *vec_in, double *vec_out)
{
  double sin_polar, cos_polar, radius, azimuth;

  radius = vec_in[MATH_RADIUS];
  sin_polar = sin(vec_in[MATH_POLAR]);
  cos_polar = cos(vec_in[MATH_POLAR]);
  azimuth = vec_in[MATH_AZIMUTH];
  vec_out[MATH_X] = radius * sin_polar * cos(azimuth);
  vec_out[MATH_Y] = radius * sin_polar * sin(azimuth);
  vec_out[MATH_Z] = radius * cos_polar;
}

/** PseudoSphericalToCartesian.

    - Converts a 3D vector from pseudo spherical (radius, polar angle cosine, azimuth) to Cartesian form.
    - Initial and final vectors may be the same. **/

void PseudoSphericalToCartesian(double *vec_in, double *vec_out)
{
  double cos_polar, sin_polar, radius, azimuth;

  radius = vec_in[MATH_RADIUS];
  azimuth = vec_in[MATH_AZIMUTH];
  cos_polar = vec_in[MATH_MU];
  sin_polar = 1.0 - cos_polar * cos_polar;
  if (sin_polar <= 0.0) /* Correct eventual round-off errors. */
    sin_polar = 0.0;
  else
    sin_polar = sqrt(sin_polar);
  vec_out[MATH_X] = radius * sin_polar * cos(azimuth);
  vec_out[MATH_Y] = radius * sin_polar * sin(azimuth);
  vec_out[MATH_Z] = radius * cos_polar;
}

/** CartesianToSpherical.

    - Converts a 3D vector from Cartesian to spherical (radius, polar angle, azimuth) form.
    - Initial and final vectors may be the same. **/

void CartesianToSpherical(double *vec_in, double *vec_out)
{
  double radius, azimuth, polar;

  radius = Norm(vec_in, 3);
  if (radius > 0.0)
    polar = acos(vec_in[MATH_Z] / radius);
  else
    polar = 0.0;
  azimuth = atan2(vec_in[MATH_Y], vec_in[MATH_X]);
  vec_out[MATH_RADIUS] = radius;
  vec_out[MATH_POLAR] = polar;
  vec_out[MATH_AZIMUTH] = azimuth;
}

/** CartesianToPseudoSpherical.

    - Converts a 3D vector from Cartesian to pseudo spherical (radius, polar angle cosine, azimuth) form.
    - Initial and final vectors may be the same. **/

void CartesianToPseudoSpherical(double *vec_in, double *vec_out)
{
  double radius, cos_polar, azimuth;

  radius = Norm(vec_in, 3);
  if (radius > 0.0)
    cos_polar = vec_in[MATH_Z] / radius;
  else
    cos_polar = 1.0;
  azimuth = atan2(vec_in[MATH_Y], vec_in[MATH_X]);
  vec_out[MATH_RADIUS] = radius;
  vec_out[MATH_MU] = cos_polar;
  vec_out[MATH_AZIMUTH] = azimuth;
}

/** Rotate.

    - Rotates a vector.
    - Vector is rotated counterclockwise (or axes clockwise) as seen from the tip of the rotation axis specified by argument 'axis'.
    - Initial and final vectors may be the same. **/

void Rotate(double *vec_in, double *vec_out, unsigned axis, double angle)
{
  unsigned i, j;
  double cos_angle, sin_angle, x, y;

  cos_angle = cos(angle);
  sin_angle = sin(angle);
  axis %= 3;
  vec_out[axis] = vec_in[axis];
  i = (axis + 1) % 3;
  j = (axis + 2) % 3;
  x = vec_in[i];
  y = vec_in[j];
  vec_out[i] = x * cos_angle - y * sin_angle;
  vec_out[j] = x * sin_angle + y * cos_angle;
}

/** GetRandomSeed.

    - Provides a random seed by reading garbage from /dev/urandom. **/

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
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

/** MetropolisHastings.

    - Constructs a Markov chain using Metropolis-Hastings sampling from a given PDF ('sampling_pdf').
    - 'states' contains the actual Markov chain, comprising 'chain_size' elements. Each element contains 'state_size' doubles. The initial state must be stored in
      the first element of 'states'. 'states' must be allocated by the calling function.
    - 'proposals' is the chain of proposals derived from each state (same length and element size as 'states'). It must be allocated by the calling function.
    - 'lower_bound' and 'upper_bound' must contain the absolute bounds for the variables defining the state.
    - 'step' must contain the parameter defining the size of the transition step.
    - 'sampling_pdf' must point to a function computing the PDF sampled by the Markov chain. Its only argument is a state of the chain; any other parameter must be
      passed through global variables.
    - 'transition_pdf' must point to a function computing the PDF of a Metropolis-Hastings transition, given a Markov chain state. Its arguments are the final and
      initial states of the transition, in this order. Any other parameter must be passed through global variables.
    - 'transition' must point to a function computing the proposed final state of a Metropolis-Hastings transition, given a Markov chain state. Its arguments are the
      final and initial states of the transition, in this order. Any other parameter must be passed through global variables.
    - 'proposals_PDF' is an array of doubles of the same length of the Markov chain. It must be allocated by the calling function if storage of sampling PDF values
      for the proposals is desired, or set to NULL otherwise.
    - 'states_PDF' is an array of doubles of the same length of the Markov chain. It must be allocated by the calling function if storage of sampling PDF values
      for the states is desired, or set to NULL otherwise.
    - If the initial state set in 'states' has a null value of the sampling PDF, the first transition probability is set to unity. Once a state with a non-null value
      of the sampling PDF is reached, the sampling PDF values do not vanish over the rest of the Markov chain by construction.
    - The order of the state and proposal chains is the following: any element of 'proposals' is sampled from the corresponding element of 'states', whereas any
      element of 'states' (except the first one) is the previous element of either 'states' or 'proposals'.
    - If 'save_states' is set to null, only the final state and proposal are given. Otherwise, the full chains are saved.
    - If 'log_flag' is non null, 'sampling_pdf' is assumed to give the natural logarithm of the PDF.
    - Returns the acceptance ratio of the chain. **/

double MetropolisHastings(double *states, double *proposals, double *proposals_pdf, double *states_pdf, unsigned state_size, unsigned chain_size, unsigned save_states,
                          unsigned log_flag, double *lower_bound, double *upper_bound, double *step,
			  void (*transition)(double *, double *, double *, double *, double *, unsigned),
                          double (*transition_pdf)(double *, double *, double *, double *, double *, unsigned), double (*sampling_pdf)(double *))
{
  unsigned i, accepted_transitions, state_size_bytes;
  double direct_transition_pdf, inverse_transition_pdf, pdf_ratio, state_pdf, proposal_pdf, *proposal, *state, *old_state, acceptance_ratio;

  /* Set initial stuff. */

  state_size_bytes = state_size * sizeof(double);
  proposal = proposals;
  state = states;
  state_pdf = sampling_pdf(state);
  if (states_pdf)
    states_pdf[0] = state_pdf;

  /* Construct the Markov chain. */

  accepted_transitions = 0;
  chain_size--;
  for (i = 0; i < chain_size; i++)
  {
    /* Set new proposal. */

    transition(proposal, state, lower_bound, upper_bound, step, state_size);
    proposal_pdf = sampling_pdf(proposal);
    if (proposals_pdf && save_states)
      proposals_pdf[i] = proposal_pdf;

    /* Compute transition probabilities. */

    direct_transition_pdf = transition_pdf(proposal, state, lower_bound, upper_bound, step, state_size);
    inverse_transition_pdf = transition_pdf(state, proposal, lower_bound, upper_bound, step, state_size);
    if (!log_flag && state_pdf <= 0.0)
      pdf_ratio = 1.1;
    else if (direct_transition_pdf <= 0.0)
      pdf_ratio = -0.1;
    else if (log_flag)
      pdf_ratio = exp(proposal_pdf - state_pdf) * inverse_transition_pdf / direct_transition_pdf;
    else
      pdf_ratio = proposal_pdf * inverse_transition_pdf / (state_pdf * direct_transition_pdf);

    /* Accept of reject transition. */

    if (save_states)
    {
      old_state = state;
      state += state_size;
    }
    if (gsl_rng_uniform(RandomNumberGenerator) < pdf_ratio)
    {
      memcpy(state, proposal, state_size_bytes);
      state_pdf = proposal_pdf;
      accepted_transitions++;
    }
    else if (save_states)
      memcpy(state, old_state, state_size_bytes);
    if (save_states)
    {
      proposal += state_size;
      if (states_pdf)
        states_pdf[i+1] = state_pdf;
    }
  }

  /* Compute final stuff. */

  transition(proposal, state, lower_bound, upper_bound, step, state_size);
  acceptance_ratio = accepted_transitions / (double) chain_size;
  if (proposals_pdf && save_states)
    proposals_pdf[chain_size] = sampling_pdf(proposal);
  return acceptance_ratio;
}

/** ParallelTempering.

    - Constructs a Markov chain using Metropolis-Hastings sampling from a given PDF ('sampling_pdf'), with parallel tempering.
    - 'states' contains the cold Markov chain, comprising 'chain_size' elements. Each element contains 'state_size' doubles. The initial state must be stored in
      the first element of 'states'. 'states' must be allocated by the calling function.
    - 'proposals' is the chain of proposals derived from each state (same length and element size as 'states'). It must be allocated by the calling function.
    - 'lower_bound' and 'upper_bound' must contain the absolute bounds for the variables defining the state.
    - 'step' must contain the parameter defining the size of the transition step.
    - 'sampling_pdf' must point to a function computing the PDF sampled by the cold Markov chain. Its only argument is a state of the chain; any other parameter must
      be passed through global variables.
    - 'transition_pdf' must point to a function computing the PDF of a Metropolis-Hastings transition, given a Markov chain state. Its arguments are the final and
      initial states of the transition, in this order. Any other parameter must be passed through global variables.
    - 'transition' must point to a function computing the proposed final state of a Metropolis-Hastings transition, given a Markov chain state. Its arguments are the
      final and initial states of the transition, in this order. Any other parameter must be passed through global variables.
    - 'proposals_PDF' is an array of doubles of the same length of the Markov chain. It must be allocated by the calling function if storage of sampling PDF values
      for the proposals is desired, or set to NULL otherwise.
    - 'states_PDF' is an array of doubles of the same length of the Markov chain. It must be allocated by the calling function if storage of sampling PDF values
      for the states is desired, or set to NULL otherwise.
    - If the initial state set in 'states' has a null value of the sampling PDF, the first transition probability is set to unity. Once a state with a non-null value
      of the sampling PDF is reached, the sampling PDF values do not vanish over the rest of the Markov chain by construction.
    - The order of the state and proposal chains is the following: any element of 'proposals' is sampled from the corresponding element of 'states', whereas any
      element of 'states' (except the first one) is the previous element of either 'states' or 'proposals'.
    - If 'save_states' is set to null, only the final state and proposal are given. Otherwise, the full chains are saved.
    - If 'log_flag' is non null, 'sampling_pdf' is assumed to give the natural logarithm of the PDF.
    - The number of hot chains in the parallel-tempering scheme is 'hot_chains_count', which have tempering factors given by 'tempering_factors'.
    - 'exchange_probability' gives the probability that two chains are interchanged in a given step.
    - Returns the acceptance ratio of the chain. **/

double ParallelTempering(double *states, double *proposals, double *proposals_pdf, double *states_pdf, unsigned state_size, unsigned chain_size, unsigned save_states,
                         unsigned log_flag, double *lower_bound, double *upper_bound, double *step,
                         void (*transition)(double *, double *, double *, double *, double *, unsigned),
                         double (*transition_pdf)(double *, double *, double *, double *, double *, unsigned), double (*sampling_pdf)(double *),
                         unsigned hot_chains_count, double *tempering_factors, double exchange_probability)
{
  unsigned i, j, accepted_transitions, state_size_bytes;
  double direct_transition_pdf, inverse_transition_pdf, pdf_ratio, state_pdf, proposal_pdf, *proposal, *state, *old_state, acceptance_ratio, x, *aux_state;
  double *cool_state, hot_proposals_pdf, factor, *hot_proposals, *hot_states, *hot_states_pdf, *cool_state_pdf;

  /* Set initial stuff. */

  state_size_bytes = state_size * sizeof(double);
  proposal = proposals;
  state = states;
  state_pdf = sampling_pdf(state);
  if (states_pdf)
    states_pdf[0] = state_pdf;
  hot_proposals = malloc(hot_chains_count * state_size_bytes);
  hot_states = malloc(hot_chains_count * state_size_bytes);
  aux_state = malloc(state_size_bytes);
  for (i = 0; i < hot_chains_count; i++)
    memcpy(hot_states + i * state_size, state, state_size_bytes); /* Copy initial state of cold chain to all hot chains. */
  hot_states_pdf = malloc(hot_chains_count * sizeof(double));
  for (i = 0; i < hot_chains_count; i++)
  {
    if (log_flag)
      hot_states_pdf[i] = state_pdf * tempering_factors[i];
    else
      hot_states_pdf[i] = pow(state_pdf, tempering_factors[i]);
  }

  /* Construct the Markov chains. */

  accepted_transitions = 0;
  chain_size--;
  for (i = 0; i < chain_size; i++)
  {
    /* Set new proposal for cold chain. */

    transition(proposal, state, lower_bound, upper_bound, step, state_size);
    proposal_pdf = sampling_pdf(proposal);
    if (proposals_pdf && save_states)
      proposals_pdf[i] = proposal_pdf;

    /* Compute transition probabilities for cold chain. */

    direct_transition_pdf = transition_pdf(proposal, state, lower_bound, upper_bound, step, state_size);
    inverse_transition_pdf = transition_pdf(state, proposal, lower_bound, upper_bound, step, state_size);
    if (!log_flag && state_pdf <= 0.0)
      pdf_ratio = 1.1;
    else if (direct_transition_pdf <= 0.0)
      pdf_ratio = -0.1;
    else if (log_flag)
      pdf_ratio = exp(proposal_pdf - state_pdf) * inverse_transition_pdf / direct_transition_pdf;
    else
      pdf_ratio = proposal_pdf * inverse_transition_pdf / (state_pdf * direct_transition_pdf);

    /* Accept of reject transition for cold chain. */

    if (save_states)
    {
      old_state = state;
      state += state_size;
    }
    if (gsl_rng_uniform(RandomNumberGenerator) < pdf_ratio)
    {
      memcpy(state, proposal, state_size_bytes);
      state_pdf = proposal_pdf;
      accepted_transitions++;
    }
    else if (save_states)
      memcpy(state, old_state, state_size_bytes);
    if (save_states)
    {
      proposal += state_size;
      if (states_pdf)
        states_pdf[i+1] = state_pdf;
    }

    /* Hot chains. */

    for (j = 0; j < hot_chains_count; j++)
    {
      /* Set new proposal. */

      transition(hot_proposals + j * state_size, hot_states + j * state_size, lower_bound, upper_bound, step, state_size);
      hot_proposals_pdf = sampling_pdf(hot_proposals + j * state_size);
      if (log_flag)
        hot_proposals_pdf *= tempering_factors[j];
      else
        hot_proposals_pdf = pow(hot_proposals_pdf, tempering_factors[j]);

      /* Compute transition probabilities. */

      direct_transition_pdf = transition_pdf(hot_proposals + j * state_size, hot_states + j * state_size, lower_bound, upper_bound, step, state_size);
      inverse_transition_pdf = transition_pdf(hot_states + j * state_size, hot_proposals + j * state_size, lower_bound, upper_bound, step, state_size);
      if (!log_flag && hot_states_pdf[j] <= 0.0)
        pdf_ratio = 1.1;
      else if (direct_transition_pdf <= 0.0)
        pdf_ratio = -0.1;
      else if (log_flag)
        pdf_ratio = exp(hot_proposals_pdf - hot_states_pdf[j]) * inverse_transition_pdf / direct_transition_pdf;
      else
        pdf_ratio = hot_proposals_pdf * inverse_transition_pdf / (hot_states_pdf[j] * direct_transition_pdf);

      /* Accept of reject transition. */

      if (gsl_rng_uniform(RandomNumberGenerator) < pdf_ratio)
      {
        memcpy(hot_states + j * state_size, hot_proposals + j * state_size, state_size_bytes);
        hot_states_pdf[j] = hot_proposals_pdf;
      }
    }

    /* Exchange chain states. */

    if (gsl_rng_uniform(RandomNumberGenerator) < exchange_probability)
    {
      j = gsl_rng_uniform_int(RandomNumberGenerator, hot_chains_count);
      factor = tempering_factors[j];
      if (j)
      {
        factor /= tempering_factors[j-1];
        cool_state_pdf = hot_states_pdf + j - 1;
        cool_state = hot_states + (j - 1) * state_size;
      }
      else
      {
        cool_state_pdf = &state_pdf;
        cool_state = state;
      }
      if (log_flag)
        pdf_ratio = exp(*cool_state_pdf * (factor - 1.0) + hot_states_pdf[j] * (1.0 / factor - 1.0));
      else
        pdf_ratio = pow(*cool_state_pdf, factor - 1.0) * pow(hot_states_pdf[j], 1.0 / factor - 1.0);
      if (gsl_rng_uniform(RandomNumberGenerator) < pdf_ratio)
      {
        if (log_flag)
        {
          x = hot_states_pdf[j];
          hot_states_pdf[j] = *cool_state_pdf * factor;
          *cool_state_pdf = x / factor;
        }
        else
        {
          x = hot_states_pdf[j];
          hot_states_pdf[j] = pow(*cool_state_pdf, factor);
          *cool_state_pdf = pow(x, 1.0 / factor);
        }
        memcpy(aux_state, hot_states + j * state_size, state_size_bytes);
        memcpy(hot_states + j * state_size, cool_state, state_size_bytes);
        memcpy(cool_state, aux_state, state_size_bytes);
      }
    }
  }
  transition(proposal, state, lower_bound, upper_bound, step, state_size);
  acceptance_ratio = accepted_transitions / (double) chain_size;
  if (proposals_pdf && save_states)
    proposals_pdf[chain_size] = sampling_pdf(proposal);
  free(hot_proposals);
  free(hot_states);
  free(hot_states_pdf);
  free(aux_state);
  return acceptance_ratio;
}

/** TanWeights.

    - Computes the weights proposed by Tan (2006, J.Comp.Graph.Stat., 15, 735) to use Metropolis-Hastings sampling for the evaluation of integrals via the Monte
      Carlo method.
    - 'states' contains the actual Markov chain, comprising 'chain_size' elements. Each element contains 'state_size' doubles.
    - 'proposals' contains the chain of proposals derived from each state (same length and element size as 'states').
    - 'proposals_pdf' is an array of doubles of the same length of the Markov chain, containing the sampling PDF values for the proposals.
    - If the present chain is an extension of an old one, 'old_chain_size' is the size of the latter.
    - 'weights' is an array of 'chain_size' elements for storing the resulting weights. It must be allocated by the calling function. Elements of 'weights' between
      'old_chain_size' and 'chain_size' are erased before computing.
    - 'lower_bound' and 'upper_bound' must contain the absolute bounds for the variables defining the state.
    - 'step' must contain the parameter defining the size of the transition step.
    - 'sort_component' indicates the component of the state used for sorting them, to reduce the computational cost if the PDF is bound.
    - 'transition_pdf' must point to a function computing the PDF of a Metropolis-Hastings transition, given a Markov chain state. Its arguments are the final and
      initial states of the transition, in this order. Any other parameter must be passed through global variables. **/

void TanWeights(double *states, double *proposals, double *proposals_pdf, double *weights, unsigned state_size, unsigned chain_size, unsigned old_chain_size,
                double *lower_bound, double *upper_bound, double *step, double (*transition_pdf)(double *, double *, double *, double *, double *, unsigned),
                unsigned sort_component)
{
  unsigned i, j, k, m, n;
  double *proposal, **pointers, pdf;

  if (sort_component >= state_size)
    SortComponent = 0;
  else
    SortComponent = sort_component;
  memset(weights + old_chain_size, 0, (chain_size - old_chain_size) * sizeof(double));
  pointers = malloc(chain_size * sizeof(double *));
  for (i = 0; i < chain_size; i++)
    pointers[i] = states + i * state_size;
  qsort(pointers, chain_size, sizeof(double *), SortIndexesAscending);
  for (i = 0; i < chain_size; i++)
  {
    k = pointers[i] - states;
    m = k / state_size;
    proposal = proposals + k;
    for (j = i; j < chain_size; j++)
    {
      n = (pointers[j] - states) / state_size;
      if (m < old_chain_size && n < old_chain_size)
        continue; /* Avoid re-computation of weights already computed. */
      pdf = transition_pdf(proposal, pointers[j], lower_bound, upper_bound, step, state_size);
      if (pdf <= 0.0)
        break;
      weights[m] += pdf;
    }
    for (j = 1; j <= i; j++)
    {
      n = (pointers[i-j] - states) / state_size;
      if (m < old_chain_size && n < old_chain_size)
        continue;
      pdf = transition_pdf(proposal, pointers[i-j], lower_bound, upper_bound, step, state_size);
      if (pdf <= 0.0)
        break;
      weights[m] += pdf;
    }
    if (weights[m] > 0.0)
      weights[m] = proposals_pdf[m] / weights[m];
  }
  free(pointers);
}

/** TanNormalisation.

    - Computes the PDF normalisation proposed by Tan (2006, J.Comp.Graph.Stat., 15, 735), together with an estimate of its undertainty ('error').
    - 'weights' is the array of weights (of size 'element_count') computed by 'TanWeights'. **/

double TanNormalisation(double *weights, unsigned element_count, double *error)
{
  unsigned i;
  double norm, **pointers, x;

  /* Perform weighted sum. */

  pointers = malloc(element_count * sizeof(double *));
  for (i = 0; i < element_count; i++)
    pointers[i] = weights + i;
  qsort(pointers, element_count, sizeof(double *), SortIndexesAscending);
  norm = 0.0;
  *error = 0.0;
  for (i = 0; i < element_count; i++)
  {
    x = pointers[i][0];
    norm += x;
    *error += x * x;
  }
  *error = sqrt(*error - norm * norm / element_count);
  free(pointers);
  return norm;
}

/** TanMean.

    - Computes the mean of a function as proposed by Tan (2006, J.Comp.Graph.Stat., 15, 735), together with an estimate of its undertainty ('error').
    - 'weights' is the array of weights (of size 'element_count') computed by 'TanWeights'.
    - 'norm' is the PDF normalisation computed by 'TanNormalisation'.
    - 'function' is the array of function values (of size 'element_count') at the proposal states. **/

double TanMean(double *weights, double *function, double norm, unsigned element_count, double *error)
{
  unsigned i;
  double mean, x;

  /* Perform weighted sum. */

  mean = 0.0;
  *error = 0.0;
  for (i = 0; i < element_count; i++)
    mean += weights[i] * function[i];
  mean /= norm;
  for (i = 0; i < element_count; i++)
  {
    x = (function[i] - mean) * weights[i];
    *error += x * x;
  }
  *error = sqrt(*error) / norm;
  return mean;
}

/** UniformTransition.

    - Samples an n-dimensional transition in a Markov chain, from a uniform PDF.
    - 'state' contains the present state of the Markov chain, comprising 'state_size' doubles.
    - 'proposal' will store the resulting proposed state (same size as 'state').
    - 'lower_bound' and 'upper_bound' must contain the absolute bounds for the variables defining the state.
    - 'step' must contain the parameter defining the size of the transition step. **/

void UniformTransition(double *proposal, double *state, double *lower_bound, double *upper_bound, double *step, unsigned state_size)
{
  unsigned i;
  double min, max;

  for (i = 0; i < state_size; i++)
  {
    if (step[i] <= 0.0)
      proposal[i] = state[i];
    else
    {
      min = state[i] - 0.5 * step[i];
      max = min + step[i];
      if (min < lower_bound[i])
        min = lower_bound[i];
      if (max > upper_bound[i])
        max = upper_bound[i];
      proposal[i] = min + (max - min) * gsl_rng_uniform(RandomNumberGenerator);
    }
  }
}

/** UniformTransitionPDF.

    - Computes the n-dimensional PDF value of a Markov chain transition, using a uniform PDF.
    - 'state' contains the state of the Markov chain from which the transition occurs, comprising 'state_size' doubles.
    - 'proposal' contains the proposed state of the chain after the transition (same size as 'state').
    - 'lower_bound' and 'upper_bound' must contain the absolute bounds for the variables defining the state.
    - 'step' must contain the parameter defining the size of the transition step. **/

double UniformTransitionPDF(double *proposal, double *state, double *lower_bound, double *upper_bound, double *step, unsigned state_size)
{
  unsigned i;
  double min, max, density;

  density = 1.0;
  for (i = 0; i < state_size; i++)
  {
    if (step[i] <= 0.0)
      continue;
    min = state[i] - 0.5 * step[i];
    max = min + step[i];
    if (min < lower_bound[i])
      min = lower_bound[i];
    if (max > upper_bound[i])
      max = upper_bound[i];
    if (proposal[i] < min || proposal[i] > max)
      return 0.0;
    density /= max - min;
  }
  return density;
}

/** PowerLawDeviate.

    - Generates a power-law deviate with an index 'index' in any closed interval ['min', 'max']. **/

double PowerLawDeviate(double min, double max, double index)
{
  double deviate, std_deviate;

  std_deviate = gsl_rng_uniform(RandomNumberGenerator);
  if (index == -1.0)
    deviate = min * exp(std_deviate * log(max / min));
  else if (index == 0.0)
    deviate = min + (max - min) * std_deviate;
  else
  {
    index += 1.0;
    min = pow(min, index);
    max = pow(max, index);
    deviate = pow(min + (max - min) * std_deviate, 1.0 / index);
  }
  return deviate;
}

/** IsotropicDirection.

    - Generates a random direction isotropically distributed within a cone of a given semi-aperture 'semi_aperture'.
    - The output direction has unit norm, by construction.
    - Relies on the fact that, for an isotropic distribution of directions in three dimensions, the projection onto any axis is distributed uniformly in [-1,1]. **/

void IsotropicDirection(double semi_aperture, double *direction)
{
  double cos_polar, sin_polar, azimuth, cos_polar_max;

  /* Sample polar angle cosine with uniform PDF in [cos_polar_max, 1]. */

  cos_polar_max = cos(semi_aperture);
  cos_polar = cos_polar_max + gsl_rng_uniform(RandomNumberGenerator) * (1.0 - cos_polar_max);
  sin_polar = 1.0 - cos_polar * cos_polar;
  if (sin_polar > 0.0) /* Correct eventual round-off errors. */
    sin_polar = sqrt(sin_polar);
  else
    sin_polar = 0.0;

  /* Sample azimuth with uniform PDF in [0, 2 pi). */

  azimuth = MATH_2PI * gsl_rng_uniform(RandomNumberGenerator);

  /* Compute direction. */

  direction[MATH_X] = sin_polar * cos(azimuth);
  direction[MATH_Y] = sin_polar * sin(azimuth);
  direction[MATH_Z] = cos_polar;
}

/** StandardGaussianDeviates.

    - Generates a couple of standard gaussian deviates.
    - Deviates are stored in 'deviates', which must be allocated by the calling program. **/

void StandardGaussianDeviates(double *deviates)
{
  double x, y, r;

  y = 2.0 * gsl_rng_uniform(RandomNumberGenerator) - 1.0;
  do
  {
    x = y;
    y = 2.0 * gsl_rng_uniform(RandomNumberGenerator) - 1.0;
    r = x * x + y * y;
  }
  while (r == 0.0 || r >= 1.0);
  r = sqrt(-2.0 * log(r) / r);
  deviates[0] = x * r;
  deviates[1] = y * r;
}

/** GaussianDeviates.

    - Generates a couple of gaussian deviates with mean 'mean' and standard deviation 'std_dev'.
    - Deviates are stored in 'deviates', which must be allocated by the calling program. **/

void GaussianDeviates(double *deviates, double mean, double std_dev)
{
  double x, y, r;

  y = 2.0 * gsl_rng_uniform(RandomNumberGenerator) - 1.0;
  do
  {
    x = y;
    y = 2.0 * gsl_rng_uniform(RandomNumberGenerator) - 1.0;
    r = x * x + y * y;
  }
  while (r == 0.0 || r >= 1.0);
  r = sqrt(-2.0 * log(r) / r) * std_dev;
  deviates[0] = x * r + mean;
  deviates[1] = y * r + mean;
}

/** ComplexGaussianDeviates.

    - Generates a gaussian deviate with mean 'mean' and standard deviation 'std_dev' in the complex plane. **/

complex ComplexGaussianDeviates(complex mean, complex std_dev)
{
  double x, y, r;
  complex deviate;

  y = 2.0 * gsl_rng_uniform(RandomNumberGenerator) - 1.0;
  do
  {
    x = y;
    y = 2.0 * gsl_rng_uniform(RandomNumberGenerator) - 1.0;
    r = x * x + y * y;
  }
  while (r == 0.0 || r >= 1.0);
  r = sqrt(-2.0 * log(r) / r);
  deviate = mean + r * (x * creal(std_dev) + I * y * cimag(std_dev));
  return deviate;
}

#endif

/*** End of 'mathematics.h'. ***/
