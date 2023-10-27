#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#include <shmem.h>
#include "wctimer.h"

// #define DEBUG
// int DEBUG = 1;
#ifdef DEBUG
#define dbg_printf(...) dbg_printf_real(__VA_ARGS__)
#else
#define dbg_printf(...)
#endif
#define dprintf(...) dbg_printf_real(__VA_ARGS__)

#define N 10000
#define G 6.67e-11
#define TIMESTEP 0.25
#define NSTEPS 10

// types
struct body_s
{
  double x;
  double y;
  double z;
  double dx;
  double dy;
  double dz;
  double fx;
  double fy;
  double fz;
  double mass;
};
typedef struct body_s body_t;

struct global_s
{
  int rank;
  int nproc;         // # of parallel processes
  int n;             // N = number of bodies in simulation
  int nsteps;        // # of timesteps to run simulation
  int bodiesPerProc; // # bodies per proc
  int remainder;     // # remainder
};
typedef struct global_s global_t;

/*
 *  global data structure
 */
global_t g;

/*
 *  prototypes
 */
int dbg_printf_real(const char *format, ...);
int eprintf(const char *format, ...);
void print_body(body_t *b);
static inline double dist(double dx, double dy, double dz);

/**
 * dist - computes distance between two bodies
 *  @param dx x-coordinate difference
 *  @param dy y-coordinate difference
 *  @param dz z-coordinate difference
 *  @returns distance magnitude
 */
static inline double dist(double dx, double dy, double dz)
{
  return sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

/**
 *  dbg_printf - debug printing wrapper, only enabled if DEBUG defined
 *    @param format printf-like format string
 *    @param ... arguments to printf format string
 *    @return number of bytes written to stderr
 */
int dbg_printf_real(const char *format, ...)
{
  va_list ap;
  int ret, len;
  char buf[1024], obuf[1024];

  va_start(ap, format);
  ret = vsprintf(buf, format, ap);
  va_end(ap);
  len = sprintf(obuf, "%4d: %s", g.rank, buf);
  write(STDOUT_FILENO, obuf, len);
  return ret;
}

/**
 *  eprintf - error printing wrapper
 *    @param format printf-like format string
 *    @param ... arguments to printf format string
 *    @return number of bytes written to stderr
 */
int eprintf(const char *format, ...)
{
  va_list ap;
  int ret;
  char buf[1024];

  if (g.rank == 0)
  {
    va_start(ap, format);
    ret = vsprintf(buf, format, ap);
    va_end(ap);
    write(STDOUT_FILENO, buf, ret);
    return ret;
  }
  else
    return 0;
}

/**
 * print_body - prints the contents of a body
 *  @param b - body to print
 */
void print_body(body_t *b)
{
  dprintf("x: %7.3f y: %7.3f z: %7.3f dx: %7.3f dy: %7.3f dz: %7.3f\n",
          b->x, b->y, b->z, b->dx, b->dy, b->dz);
}

/**
 *  Init initial value for bodies
 *  @param bodies - ptr to bodies needed initialization
 */
static inline void init(body_t *bodies)
{
  for (int i = 0; i < g.bodiesPerProc; i++)
  {
    int index = i + g.bodiesPerProc * g.rank;
    bodies[i].x = 100.0 * (index + 0.1);
    bodies[i].y = 200.0 * (index + 0.1);
    bodies[i].z = 300.0 * (index + 0.1);
    bodies[i].dx = index + 400.0;
    bodies[i].dy = index + 500.0;
    bodies[i].dz = index + 600.0;
    bodies[i].fx = 0;
    bodies[i].fy = 0;
    bodies[i].fz = 0;
    bodies[i].mass = 10e6 * (index + 100.2);
  }
}

/**
 *  Update fx fy fz of a body
 *  @param me - ptr to affected body
 *  @param calculating_body - ptr to affecting body
 */
static inline void updateForces(body_t *me, body_t *calculating_body)
{
  double d, f;       // distance, force
  double dx, dy, dz; // position deltas

  // compute the distances in each dimension
  dx = me->x - calculating_body->x;
  dy = me->y - calculating_body->y;
  dz = me->z - calculating_body->z;

  // compute the distance magnitude
  d = dist(dx, dy, dz);

  // skip over ourselves (d==0)
  if (d != 0)
  {
    // F = G m1 m2 / r^2
    f = (G * me->mass * calculating_body->mass) / (d * d);

    // compute force components in each dimension
    me->fx += (f * dx) / d;
    me->fy += (f * dy) / d;
    me->fz += (f * dz) / d;
  }
}

/**
 *  Update the state at t+1 of a body
 *  @param me - ptr to t body
 *  @param nextme - ptr to t+1 body
 */
static inline void updateState(body_t *me, body_t *nextme)
{
  double ax, ay, az; // acceleration components
  // acc = force / mass (F=ma)
  ax = me->fx / me->mass;
  ay = me->fy / me->mass;
  az = me->fz / me->mass;

  // update the body velocity at time t+1
  nextme->dx = me->dx + (TIMESTEP * ax);
  nextme->dy = me->dy + (TIMESTEP * ay);
  nextme->dz = me->dz + (TIMESTEP * az);

  // update the body position at t+1
  nextme->x = me->x + (TIMESTEP * me->dx);
  nextme->y = me->y + (TIMESTEP * me->dy);
  nextme->z = me->z + (TIMESTEP * me->dz);

  // update fx, fy, fz at t+1
  nextme->fx = 0;
  nextme->fy = 0;
  nextme->fz = 0;

  // copy over the mass
  nextme->mass = me->mass;
}

int main(int argc, char **argv)
{
  char c;
  wc_timer_t ttimer; // total time
  wc_timer_t itimer; // per-iteration timer

  memset(&g, 0, sizeof(g)); // zero out global data structure
  shmem_init();             // initialize OpenSHMEM

  g.nproc = shmem_n_pes();
  g.rank = shmem_my_pe();

  wc_tsc_calibrate();

  g.n = N;
  g.nsteps = NSTEPS;

  while ((c = getopt(argc, argv, "hn:t:")) != -1)
  {
    switch (c)
    {
    case 'h':
      eprintf("usage: lab3 [-n #bodies] [-t #timesteps]\n");
      shmem_finalize();
      exit(0);
      break;
    case 'n':
      g.n = atoi(optarg);
      break;
    case 't':
      g.nsteps = atoi(optarg);
      break;
    }
  }

  eprintf("beginning N-body simulation of %d bodies with %d processes over %d timesteps\n",
          g.n, g.nproc, g.nsteps);

  // calculate bodies per processor
  g.bodiesPerProc = g.n / g.nproc;
  g.remainder = g.n % g.nproc;

  // initialize currect bodies, next bodies' state, communicating bodies
  body_t *communicating_bodies = shmem_calloc(g.bodiesPerProc + 1, sizeof(body_t));
  body_t *bodies = shmem_calloc(g.bodiesPerProc + 1, sizeof(body_t));
  body_t *nexts = shmem_calloc(g.bodiesPerProc + 1, sizeof(body_t));
  body_t *calculating_bodies = communicating_bodies;

  // init bodies and remainder
  init(bodies);
  // if has remainder, distribute 1 to each remainder first processor
  if (g.rank < g.remainder)
  {
    int index = g.rank + g.bodiesPerProc * g.nproc;
    bodies[g.bodiesPerProc].x = 100.0 * (index + 0.1);
    bodies[g.bodiesPerProc].y = 200.0 * (index + 0.1);
    bodies[g.bodiesPerProc].z = 300.0 * (index + 0.1);
    bodies[g.bodiesPerProc].dx = index + 400.0;
    bodies[g.bodiesPerProc].dy = index + 500.0;
    bodies[g.bodiesPerProc].dz = index + 600.0;
    bodies[g.bodiesPerProc].mass = 10e6 * (index + 100.2);
  }

  WC_INIT_TIMER(ttimer);
  WC_START_TIMER(ttimer);

  // main simulation loop
  for (int ts = 0; ts < g.nsteps; ts++)
  {
    WC_INIT_TIMER(itimer);
    WC_START_TIMER(itimer);

    // loop through all bodies from other processors
    for (int pe = 0; pe < g.nproc; pe++)
    {

      // if calculate with this processor's bodies
      if (pe == g.rank)
      {
        calculating_bodies = bodies;
      }
      // if calculating with other processor's bodies
      else
      {
        shmem_getmem(communicating_bodies, bodies, (g.bodiesPerProc + 1) * sizeof(body_t), pe);
        calculating_bodies = communicating_bodies;
      }

      // loop through this processor's bodies
      for (int i = 0; i < g.bodiesPerProc + (g.rank < g.remainder); i++)
      {
        // loop through the communicating bodies
        for (int j = 0; j < g.bodiesPerProc + (pe < g.remainder); j++)
        {
          updateForces(&bodies[i], &calculating_bodies[j]);
        }
      }
    }

    // update state to t+1 of this processor's bodies
    for (int i = 0; i < g.bodiesPerProc + (g.rank < g.remainder); i++)
    {
      updateState(&bodies[i], &nexts[i]);
    }

    shmem_barrier_all();

    // copy the t+1 state to be the new time t
    memcpy(bodies, nexts, (g.bodiesPerProc + 1) * sizeof(body_t));

    WC_STOP_TIMER(itimer);
    eprintf("timestep %d complete: %7.4f ms\n", ts, WC_READ_TIMER_MSEC(itimer));
  }

  WC_STOP_TIMER(ttimer);
  eprintf("execution time: %7.4f ms\n", WC_READ_TIMER_MSEC(ttimer));

  shmem_finalize(); // finalize OpenSHMEM
}
