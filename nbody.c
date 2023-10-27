#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#define N 10000
#define NSTEPS 100
// #define N 300000
// #define NSTEPS 5
#define G 6.67e-11
#define TIMESTEP 0.25

/*
 * body data structure
 */
struct body_s
{
	double x;
	double y;
	double z;
	double dx;
	double dy;
	double dz;
	double mass;
};
typedef struct body_s body_t;

/*
 * function prototypes
 */
void init(void);
double dist(double dx, double dy, double dz);

body_t bodies[N]; // array of N-bodies at timestep t
body_t next[N];	  // array of N-bodies at timestep t+1

/**
 * init - give the planets initial values for position, velocity, mass
 */
void init(void)
{
	for (int i = 0; i < N; i++)
	{
		bodies[i].x = 100.0 * (i + 0.1);
		bodies[i].y = 200.0 * (i + 0.1);
		bodies[i].z = 300.0 * (i + 0.1);
		bodies[i].dx = i + 400.0;
		bodies[i].dy = i + 500.0;
		bodies[i].dz = i + 600.0;
		bodies[i].mass = 10e6 * (i + 100.2);
	}
}

/**
 * dist - determine the distance between two bodies
 *    @param dx - distance in the x dimension
 *    @param dy - distance in the y dimension
 *    @param dz - distance in the z dimension
 *    @return distance
 */
double dist(double dx, double dy, double dz)
{
	return sqrt((dx * dx) + (dy * dy) + (dz * dz));
	;
}

/**
 * computeforce - compute the superposed forces on one body
 *   @param me     - the body to compute forces on at time t
 *   @param nextme - the body at time t+1
 */
void computeforce(body_t *me, body_t *nextme)
{
	double d, f;	   // distance, force
	double dx, dy, dz; // position deltas
	double fx, fy, fz; // force components
	double ax, ay, az; // acceleration components

	fx = fy = fz = 0.0;

	// for every other body relative to me
	for (int i = 0; i < N; i++)
	{
		// compute the distances in each dimension
		dx = me->x - bodies[i].x;
		dy = me->y - bodies[i].y;
		dz = me->z - bodies[i].z;

		// compute the distance magnitude
		d = dist(dx, dy, dz);

		// skip over ourselves (d==0)
		if (d != 0)
		{

			// F = G m1 m2 / r^2
			f = (G * me->mass * bodies[i].mass) / (d * d);

			// compute force components in each dimension
			fx += (f * dx) / d;
			fy += (f * dy) / d;
			fz += (f * dz) / d;
		}
	}

	// acc = force / mass (F=ma)
	ax = fx / me->mass;
	ay = fy / me->mass;
	az = fz / me->mass;

	// update the body velocity at time t+1
	nextme->dx = me->dx + (TIMESTEP * ax);
	nextme->dy = me->dy + (TIMESTEP * ay);
	nextme->dz = me->dz + (TIMESTEP * az);

	// update the body position at t+1
	nextme->x = me->x + (TIMESTEP * me->dx);
	nextme->y = me->y + (TIMESTEP * me->dy);
	nextme->z = me->z + (TIMESTEP * me->dz);

	// copy over the mass
	nextme->mass = me->mass;
}

/**
 *  print_body - prints a body for debugging
 *    @param b - body to print
 */
void print_body(body_t *b)
{
	printf("x: %7.3f y: %7.3f z: %7.3f dx: %7.3f dy: %7.3f dz: %7.3f\n",
		   b->x, b->y, b->z, b->dx, b->dy, b->dz);
}

/*
 * get_wctime - returns wall clock time as double
 *   @return double representation of wall clock time
 */
double get_wctime(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

/**
 * main
 */
int main(int argc, char **argv)
{
	double start, tsstart;

	// setup initial conditions
	init();

	printf("beginning N-body simulation of %d bodies.\n", N);

	start = get_wctime();
	// for each timestep in the simulation
	for (int ts = 0; ts < NSTEPS; ts++)
	{
		tsstart = get_wctime();
		// for every body in the universe
		for (int i = 0; i < N; i++)
		{
			computeforce(&bodies[i], &next[i]);
		}

		// copy the t+1 state to be the new time t
		for (int i = 0; i < N; i++)
		{
			memcpy(&bodies[i], &next[i], sizeof(body_t));
		}
		printf("timestep %d complete: %7.3f ms\n", ts, (get_wctime() - tsstart) * 1000);
	}
	printf("simulation complete: %9.3f ms\n", (get_wctime() - start) * 1000);
	printf("velocity: %f\n", bodies[0].dx);
	printf("bodies[%d] mass: %f\n", N - 1, bodies[N - 1].mass);
	printf("bodies[0] x - y - z: %f, %f, %f\n", bodies[0].x, bodies[0].y, bodies[0].z);
	printf("bodies[%d] x - y - z: %f, %f, %f\n", N - 1, bodies[N - 1].x, bodies[N - 1].y, bodies[N - 1].z);
	return 0;
}