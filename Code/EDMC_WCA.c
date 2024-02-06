#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#define MAXNEIGH 24


//Include random number generator 
#include "xoshiro256plus.c"

//Include algorithm for exponential random numbers 
#include "ran_exp.c"

//Header file
#include "EDMC_WCA.h"

//Number of extra events (e.g. write) to allocate space for
#define EXTRAEVENTS 12

//Pi (if not already defined)
#ifndef M_PI
#define M_PI 3.1415926535897932
#endif

// Event driven Monte Carlo for the WCA model.
// Based on the methodology of [Peters & De With, Phys. Rev. E 85, 026703 (2012)]
// Uses efficient event queue for large systems: see [G. Paul, J. Comput. Phys. 221, 615 (2007)] (but using a binary search tree instead of a complete binary heap... which should not have a significant effect)
// Events and particles are the same struct in this code



double maxtime = 100;                   //Simulation stops at this time
int makesnapshots = 1;                  //Whether to make snapshots during the run (yes = 1, no = 0)
double writeinterval = 1;               //Time between output to screen / data file
double snapshotinterval = 10;           //Time between snapshots (should be a multiple of writeinterval)
double epsilon = 20;                    //Interaction strength in units of k_B T.

int initialconfig = 0;                  //= 0 load from file, 1 = FCC crystal, 2 = simple cubic crystal, 3 = random
char inputfilename[100] = "init.sph";   //File to read as input snapshot (for initialconfig = 0)
double density = 0.8;                   //Number density (for initialconfig = 1 or 2)
int N = 1000;                           //Number of particles (for FCC/simple cubic/random)
double maxrange, maxrange2;
double sizeratio = 0.78;                //Used for random placement, generates a binary mixture with this size ratio if composition < 1.
double composition = 1.0 / 3.0;         //For random placement

uint64_t randomseed = 42;               //Seed for random number generator.

//Variables related to the event queueing system. These can affect efficiency.
//The system schedules only events in the current block of time with length "eventlisttime" into a sorted binary search tree. 
//The rest are scheduled in unordered linked lists associated with the "numeventlists" next blocks.
//"numeventlists" is roughly equal to maxscheduletime / eventlisttime
//Any events occurring even later are put into an overflow list
//After every time block with length "eventlisttime", the set of events in the next linear list is moved into the binary search try.
//All events in the overflow list are also rescheduled.

//After every "writeinterval", the code will output two listsizes to screen. 
//The first is the average number of events in the first that gets moved into the event tree after each block.
//The second is the average length of the overflow list.
//Ideally, we set maxscheduletime large enough that the average overflow list size is negligible (i.e. <10 events)
//Also, there is some optimum value for the number of events per block (scales approximately linearly with "eventlisttime").
//Choices below should be reasonable for most cases.
double maxscheduletime = 0.5;          //All events further in the future than this go in the overflow list. (Could be adjusted to tweak performance.)
double eventlisttimemultiplier = 0.2;  //event list time will be this / N (Could be adjusted to tweak performance.)
int numeventlists;
double eventlisttime;
double reftime = 0; //Reference time marking the start of the first block in the event queue. (Must = 0 at beginning.)
particle* root;     //Root of event tree


//Neighbor lists
double shellthickness = 0.2;            //Shell thickness for the neighbor list
double shellsize;



//Internal variables
double time = 0;
int currentlist = 0;
int totalevents;

int listcounter1 = 0, listcounter2 = 0, mergecounter = 0;

particle** eventlists; //Last one is overflow list


particle* particles;    //Main particle array (allocated in initparticles())
particle** celllist;    //Cell list (allocated in initcelllist())
particle* writeevent;   //Pointer to event keeping track of next write to screen/disk

double xsize, ysize, zsize;         //Box size
double hx, hy, hz;                  //Half box size
double icxsize, icysize, iczsize;   //Inverse cell size
int    cx, cy, cz;                  //Number of cells

double dvtot = 0;                   //Cumulative momentum transfer (for calculating pressure)
unsigned int colcounter = 0;        //Collision counter (could overflow in a long run, but not really important)


//Options for making snapshots on a logarithmic time scale
//This makes a number of batches of snapshots over the course of the simulations, 
//   where the time between snapshots grows exponentially.
//Useful for calculating e.g. F(q,t).
#define MAXLOGSNAPSHOTS 10000
int makelogsnapshots = 0;       //Toggle on or off
logsnaptime logsnapshottimes[MAXLOGSNAPSHOTS];
int nextlogsnapshot;
int numsnapshots = 40;		    //Number of snapshots per batch
int numbatches = 10;			//Number of batches
double minimumtime = 0.01;		//Time difference between the first two snapshots
particle* logsnapshotevent;

//The internal clock is periodically reset to improve numerical stability.
//Probably not crucial for spherical particles, but time costs are negligble.
double timeoffset = 0;                          //Offset between actual time and internal time
double resetinterval = 10000;                   //Interval after which the internal clock is reset
particle* resetevent;                           //Event keeping track of resetting

int main(int argc, char** argv)
{
    init();
    printf("Starting\n");


    while (time <= maxtime - timeoffset)
    {
      step();
    }
    time = maxtime - timeoffset;
    printstuff();
    outputsnapshot();

    free(celllist);
    free(particles);
    free(eventlists);
    return 0;
}

/**************************************************
**                 PRINTSTUFF
** Some data at the end of the simulation
**************************************************/
void printstuff()
{
    int i;
    particle* p;
    double v2tot = 0;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        v2tot += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
    }
    double ekin = 0.5*v2tot/N;
    double temp = v2tot/(3*(N-1));
    printf("Average kinetic energy: %lf\n", ekin);
    printf("Temperature           : %lf\n", temp);
    printf("Total time simulated  : %lf\n", time);

}


/**************************************************
**                    INIT
**************************************************/
void init()
{
    int i;
    //Pre-calculate maximum interaction range and the resulting shell size for neighbor list
    maxrange = pow(2.0,1.0/6.0);
    maxrange2 = maxrange*maxrange;
    shellsize = maxrange + shellthickness;

    //Initialize random number generator
    printf("Seed: %" PRIu64 "\n", randomseed);
    init_genrand(randomseed);
    random_exponential_init();


    //Initialize particles
    if      (initialconfig == 0) loadparticles();
    else if (initialconfig == 1) fcc();
    else if (initialconfig == 2) simplecubic();
    else                         randomparticles();

    for (i = 0; i < N; i++)
    {
        particle* p = particles + i;
        p->boxestraveledx = 0;
        p->boxestraveledy = 0;
        p->boxestraveledz = 0;
        p->nneigh = 0;
        p->counter = 0;
        p->t = 0;
        p->xn = p->x;
        p->yn = p->y;
        p->zn = p->z;
    }



    //Initialize velocities
    randommovement();

    //Pre-calculate half box lengths, used for periodic boundary conditions
    hx = 0.5 * xsize; hy = 0.5 * ysize; hz = 0.5 * zsize;

    //Initialize event queue
    initevents();

    //Initialize cell list
    initcelllist();

    //Initialize neighbor list
    for (i = 0; i < N; i++)
    {
        makeneighborlist(particles + i);
    }

    //Initialize logarithmic snapshots (if enabled)
    if (makelogsnapshots) initlogsnapshots();


    printf("Done with initialization\n");
}

/******************************************************
**               INITPARTICLES
** Allocates particle array
** Should only be called once the number of particles
**    is known.
** Note that an event and a particle are the same thing.
** We allocate some extra space for events that are
**    not associated with a particle.
******************************************************/
void initparticles(int n)
{
    N = n;
    totalevents = N + EXTRAEVENTS;
    particles = (particle*)malloc((totalevents) * sizeof(particle));
    if (!particles)
    {
        printf("Failed to allocate memory for particles\n");
        exit(3);
    }
    else memset(particles, 0, (totalevents) * sizeof(particle));
}


/******************************************************
**               MYGETLINE
** Reads a single line, skipping over lines
** commented out with #
******************************************************/
int mygetline(char* str, FILE* f)
{
    int comment = 1;
    while (comment)
    {
        if (!fgets(str, 255, f)) return -1;
        if (str[0] != '#') comment = 0;
    }
    return 0;
}


/**************************************************
**                    FCC
** Puts particles on an FCC lattice
** N should be 4 * a perfect cube
** Monodisperse
**************************************************/
void fcc()
{
    int i, j, k;
    particle* p;

    int ncell = cbrt(N / 4) + 0.0001;

    if (ncell * ncell * ncell * 4 != N)
    {
        printf("N should be 4 * a perfect cube! (e.g. %d)\n", ncell * ncell * ncell * 4);
        exit(3);
    }

    double volume = N / (density);
    xsize = cbrt(volume);
    ysize = xsize;
    zsize = xsize;

    double step = xsize / ncell;

    printf("step: %lf\n", step);
    initparticles(N);
    initcelllist();
    printf("Placing particles\n");

    p = particles;
    for (i = 0; i < ncell; i++) for (j = 0; j < ncell; j++) for (k = 0; k < ncell; k++)
    {
        p->x = (i + 0.25) * step;
        p->y = (j + 0.25) * step;
        p->z = (k + 0.25) * step;
        p++;
        p->x = (i + 0.75) * step;
        p->y = (j + 0.75) * step;
        p->z = (k + 0.25) * step;
        p++;
        p->x = (i + 0.75) * step;
        p->y = (j + 0.25) * step;
        p->z = (k + 0.75) * step;
        p++;
        p->x = (i + 0.25) * step;
        p->y = (j + 0.75) * step;
        p->z = (k + 0.75) * step;
        p++;
    }

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        p->radius = 0.5;
        p->type = 0;
        p->mass = 1;
    }

    printf("Packing fraction: %lf\n", M_PI / (6.0 * xsize * ysize * zsize) * N);
    printf("Starting configuration from fcc crystal\n");
}

/**************************************************
**                    SIMPLECUBIC
** Puts particles on a simple cubic lattice
** N should be a perfect cube
** Monodisperse
**************************************************/
void simplecubic()
{
    int i, j, k;
    particle* p;

    int ncell = cbrt(N) + 0.0001;

    if (ncell * ncell * ncell != N)
    {
        printf("N should be a perfect cube! (e.g. %d)\n", ncell * ncell * ncell);
        exit(3);
    }

    double volume = N / (density);
    xsize = cbrt(volume);
    ysize = xsize;
    zsize = xsize;

    double step = xsize / ncell;

    printf("step: %lf\n", step);
    initparticles(N);
    initcelllist();
    printf("Placing particles\n");

    p = particles;
    for (i = 0; i < ncell; i++) for (j = 0; j < ncell; j++) for (k = 0; k < ncell; k++)
    {
        p->x = (i + 0.25) * step;
        p->y = (j + 0.25) * step;
        p->z = (k + 0.25) * step;
        p++;
    }

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        p->radius = 0.5;
        p->type = 0;
        p->mass = 1;
    }

    printf("Packing fraction: %lf\n", M_PI / (6.0 * xsize * ysize * zsize) * N);
    printf("Starting configuration from simple cubic crystal\n");
}

/**************************************************
**                RANDOMPARTICLES
** Places particles randomly
** May cause heavy slowdown at the beginning if
** the density is too high.
** Allows for binary systems, based on 
**    sizeratio and composition.
**************************************************/
void randomparticles()
{
    int i;
    particle* p;

    double volume = N / (density);
    xsize = cbrt(volume);
    ysize = xsize;
    zsize = xsize;

    initparticles(N);
    initcelllist();
    printf("Placing particles\n");

    p = particles;
    for (i = 0; i < N; i++)
    {
        p->x = genrand_real2()*xsize;
        p->y = genrand_real2()*ysize;
        p->z = genrand_real2()*zsize;
        if (i < N * composition)
        {
            p->radius = 0.5;
            p->type = 0;
        }
        else
        {
            p->radius = sizeratio * 0.5;
            p->type = 1;
        }
        p->mass = 1;
        p++;
    }

    printf("Density: %lf\n", N / (xsize * ysize * zsize));
    printf("Starting configuration from random placement\n");
}


/**************************************************
**                    LOADPARTICLES
** Loads particles from a snapshot
** Only mono- and bidisperse
** Largest size should be diameter 1 (radius 0.5)
**************************************************/
void loadparticles()
{
    char tmp;
    int i, npart;
    particle* p;
    char buffer[255];

    FILE* file;
    file = fopen(inputfilename, "r");
    if (!file)
    {
        printf("File not found!\n");
        exit(3);
    }
    mygetline(buffer, file);
    int ftmp = sscanf(buffer, "%d", &npart);
    if (ftmp != 1) { printf("Read error (n or box)\n"); exit(3); }
    mygetline(buffer, file);
    ftmp = sscanf(buffer, "%lf %lf %lf\n", &xsize, &ysize, &zsize);
    if (ftmp != 3) { printf("Read error (n or box)\n"); exit(3); }


    initparticles(npart);
    printf("Placing particles\n");
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        mygetline(buffer, file);
        ftmp = sscanf(buffer, "%c %lf  %lf  %lf %lf\n", &tmp, &(p->x), &(p->y), &(p->z), &(p->radius));
        backinbox(p);
        if (ftmp != 5) { printf("Read error (particle) %d \n String: %s\n", ftmp, buffer); exit(3); }
        p->type = tmp - 'a';
        p->mass = 1;
    }
    fclose(file);
    initcelllist();

    if (particles[N-1].type != 0 && particles[0].type == 0)//binary?
    {
        sizeratio = particles[N-1].radius / particles[0].radius;
        printf ("Sizeratio: %lf\n", sizeratio);
    }


    printf("Density: %lf\n", N / ( xsize * ysize * zsize));

}

/**************************************************
**                RANDOMMOVEMENT
** Assigns particles random velocities (Gaussian)
** Zero center-of-mass motion
**************************************************/
void randommovement()
{
    particle* p;
    double v2tot = 0, vxtot = 0, vytot = 0, vztot = 0;
    double mtot = 0;
    int i;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        double imsq = 1.0 / sqrt(p->mass);

        p->vx = imsq * random_gaussian();
        p->vy = imsq * random_gaussian();
        p->vz = imsq * random_gaussian();
        vxtot += p->mass * p->vx;					//Keep track of total v
        vytot += p->mass * p->vy;
        vztot += p->mass * p->vz;
        mtot += p->mass;
    }


    vxtot /= mtot; vytot /= mtot; vztot /= mtot;
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->vx -= vxtot;					//Make sure v_cm = 0
        p->vy -= vytot;
        p->vz -= vztot;
        v2tot += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
    }
    double fac = sqrt(3.0 / (v2tot / (N-1)));
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->vx *= fac;					//Fix energy
        p->vy *= fac;
        p->vz *= fac;
    }
    printf("Starting configuration read from %s\n", inputfilename);
}

/**************************************************
**                UPDATE
** Update a particle to the current time.
**************************************************/
void update(particle* p1)
{
    double dt = time - p1->t;
    p1->t = time;
    p1->x += dt * p1->vx;
    p1->y += dt * p1->vy;
    p1->z += dt * p1->vz;
}

/**************************************************
**                 INITCELLLIST
** Initialize the cell list and add all particles
**************************************************/
void initcelllist()
{
    int i;
    cx = (int)(xsize - 0.0001) / shellsize;				//Set number of cells
    cy = (int)(ysize - 0.0001) / shellsize;
    cz = (int)(zsize - 0.0001) / shellsize;

    if (cx < 3 || cy < 3 || cz < 3)
    {
        printf("WARNING: Box too small for cell list: using 3 cells.");
    }
    if (cx < 3) cx = 3;
    if (cy < 3) cy = 3;
    if (cz < 3) cz = 3;
    while (cx*cy*cz > 8*N)     //Very dilute, let's make the cells larger
    {
        cx *= 0.9;
        cy *= 0.9;
        cz *= 0.9;
    }

    printf("Cells: %d, %d, %d\n", cx, cy, cz);
    celllist = (particle**) malloc(cx*cy*cz*sizeof(particle*));

    icxsize = cx / xsize;						//Set inverse cell size
    icysize = cy / ysize;
    iczsize = cz / zsize;
    memset(celllist, 0, cx*cy*cz*sizeof(particle*));
    for (i = 0; i < N; i++) 
    {
        particle* p = particles + i;
        addtocelllist(p, p->x * icxsize, p->y * icysize, p->z * iczsize);
    }
}

/**************************************************
**               CELLOFFSET
** Convert from three cell-indices to a single
**  index used for the array of cells.
**************************************************/
int celloffset(int a, int b, int c)
{
    return a + b*cx + c*cy*cy;
}


/**************************************************
**               REMOVEFROMCELLLIST
** Remove particle from cell list
**************************************************/
void removefromcelllist(particle* p1)
{
    if (p1->prev) p1->prev->next = p1->next;    //Remove particle from celllist
    else          celllist[p1->cell] = p1->next;
    if (p1->next) p1->next->prev = p1->prev;
}

/**************************************************
**                    ADDTOCELLLIST
** Add particle to cell list
** Also flag whether the particle is close to the
** box edge.
**************************************************/
void addtocelllist(particle* p, int cellx, int celly, int cellz)
{
    p->cell = celloffset(cellx, celly, cellz);
    p->next = celllist[p->cell];	//Add particle to celllist
    if (p->next) p->next->prev = p;			//Link up list
    celllist[p->cell] = p;
    p->prev = NULL;
    p->nearboxedge = (cellx == 0 || celly == 0 || cellz == 0 || cellx == cx - 1 || celly == cy - 1 || cellz == cz - 1);

}

/**************************************************
**                     STEP
** Handle the next event.
**************************************************/
void step()
{
    particle* ev;
    ev = root->right;
    while (ev == NULL) //If there are no events in the tree, move on to the next bucket of events.
    {
        addnexteventlist();
        ev = root->right;
    }

    while (ev->left) ev = ev->left;		//Find first event in tree
    if (ev->eventtime < time)
    {
        printf ("Negative time\n");
    }

    time = ev->eventtime;       //Update the time
    removeevent(ev);            //De-schedule the event we're about to handle
    switch(ev->eventtype)       //Which type?
    {
        case 0:
            collision(ev);
            break;
        case 8:
            makeneighborlist(ev);
            break;
        case 50:
            resettime();
            break;
        case 100:
            write();
            break;
        case 102:
            writelogsnapshot();
            break;   
    }
}



/**************************************************
**                MAKENEIGHBORLIST
** Update the neighbor list of p1
**************************************************/
void makeneighborlist(particle* p1)
{
    p1->counter++;
    int cdx, cdy, cdz;
    particle* p2;
    double dx, dy, dz, r2, rm;

    update(p1);

    //Check periodic boundaries
    if (p1->x >= xsize) { p1->x -= xsize; p1->boxestraveledx++; }
    else if (p1->x < 0) { p1->x += xsize; p1->boxestraveledx--; }
    if (p1->y >= ysize) { p1->y -= ysize; p1->boxestraveledy++; }
    else if (p1->y < 0) { p1->y += ysize; p1->boxestraveledy--; }
    if (p1->z >= zsize) { p1->z -= zsize; p1->boxestraveledz++; }
    else if (p1->z < 0) { p1->z += zsize; p1->boxestraveledz--; }
    
    p1->xn = p1->x;     //Update neighbor list position
    p1->yn = p1->y;
    p1->zn = p1->z;

    //Remove particle from cell list
    removefromcelllist(p1);

    //Erase particle from neighbor lists
    int i, j;
    for (i = 0; i < p1->nneigh; i++)
    {
        p2 = p1->neighbors[i];
        for (j = 0; j < p2->nneigh; j++)
        {
            if (p2->neighbors[j] == p1)
            {
                p2->nneigh--;
                p2->neighbors[j] = p2->neighbors[p2->nneigh];
                break;
            }
        }
    }

    //New position in cell list
    int cellx = p1->x * icxsize, celly = p1->y * icysize, cellz = p1->z * iczsize;
    addtocelllist(p1, cellx, celly, cellz);

    cellx += cx;        //For easier modulo calculations, ensures that (cellx-1)%cx is non-negative
    celly += cy;
    cellz += cz;
    p1->nneigh = 0;

    for (cdx = cellx - 1; cdx < cellx + 2; cdx++)
        for (cdy = celly - 1; cdy < celly + 2; cdy++)
            for (cdz = cellz - 1; cdz < cellz + 2; cdz++)
            {
                p2 = celllist[celloffset(cdx % cx, cdy % cy, cdz % cz)];
                while (p2)
                {
                    if (p2 != p1)
                    {
                        dx = p1->xn - p2->xn;   //Distance between neighbor list positions
                        dy = p1->yn - p2->yn;
                        dz = p1->zn - p2->zn;
                        if (p1->nearboxedge)
                        {
                            if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
                            if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
                            if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
                        }
                        r2 = dx * dx + dy * dy + dz * dz;
                        rm = shellsize;
                        if (r2 < rm * rm)
                        {
                            p1->neighbors[p1->nneigh++] = p2;
                            p2->neighbors[p2->nneigh++] = p1;
                        }
                    }
                    p2 = p2->next;
                }
            }

    findcollisions(p1);


}


/**************************************************
**                FINDNEIGHBORLISTUPDATE
** Assumes p1 is up to date
** Note that the particle is always in the same
** box as its neighborlist position (p->xn)
** (so no need to worry about PBCs)
**************************************************/
double findneighborlistupdate(particle* p1)
{
    double dx = p1->x - p1->xn;
    double dy = p1->y - p1->yn;
    double dz = p1->z - p1->zn;

    double dvx = p1->vx, dvy = p1->vy, dvz = p1->vz;

    double b = dx * dvx + dy * dvy + dz * dvz;                  //dr.dv

    double dv2 = dvx * dvx + dvy * dvy + dvz * dvz;
    double dr2 = dx * dx + dy * dy + dz * dz;
    double md = (shellsize - maxrange) * 0.5;

    double disc = b * b - dv2 * (dr2 - md * md);       //Should always be positive, as dr2 < md*md
    double t = (-b + sqrt(disc)) / dv2;
    return t;
}

/**************************************************
**                POTENTIAL
** Weeks-Chandler-Andersen potential
**************************************************/
double potential(double r2)
{
    if (r2 > maxrange2) return 0;
    double ir6 = 1.0/(r2*r2*r2);
    return 4*epsilon*(ir6*ir6 - ir6 + 0.25);
}

/**************************************************
**                FINDCOLLISION
** Detect the next collision for two particles
** Note that p1 is always up to date in
** findcollision
** tmin is the time of the earliest predicted
      collision for p1 so far
**************************************************/
int findcollision(particle* p1, particle* p2, double* tmin)
{
    double dt2 = time - p2->t;
    double dx = p1->x - p2->x - dt2 * p2->vx;    //relative distance at current time
    double dy = p1->y - p2->y - dt2 * p2->vy;
    double dz = p1->z - p2->z - dt2 * p2->vz;
    if (p1->nearboxedge)
    {
        if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
        if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
        if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
    }

    double dvx = p1->vx - p2->vx;                               //relative velocity
    double dvy = p1->vy - p2->vy;
    double dvz = p1->vz - p2->vz;

    double b = dx * dvx + dy * dvy + dz * dvz;                  //dr.dv
    if (b>0) return 0;      //Particles flying apart

    double dr2 = dx * dx + dy * dy + dz * dz;
    double dv2 = dvx * dvx + dvy * dvy + dvz * dvz;
    double sigmaij2 = p1->radius + p2->radius;
    sigmaij2 *= sigmaij2;
    double ucurrent = potential(dr2/sigmaij2);          //Current energy
    double unew = ucurrent + random_exponential();      //Energy at collision
    double md = cbrt(2.0/(1.0+sqrt(unew / epsilon))) * sigmaij2; //Square of collision distance

    double A = md-dr2;
    double disc = b * b + dv2 * A;
    if (disc < 0) return 0;
    double t = (-b - sqrt(disc)) / dv2;     //Well collision
    if (t < -0.00000000001)
    {
        printf ("Error\n");     //Should not happen
    }
    if (t < *tmin)              //Is this our earliest predicted collision for p1 so far?
    {
        *tmin = t;
        return 1;
    }
    return 0;        


}


/**************************************************
**                FINDCOLLISIONS
** Find all collisions for particle p1.
**************************************************/
void findcollisions(particle* p1)    //All collisions of particle p1
{
    int i;
    double tmin = findneighborlistupdate(p1);   //Time to update neighbor list
    int type = 8;
    particle* partner = p1;
    particle* p2;
    for (i = 0; i < p1->nneigh; i++)
    {
        p2 = p1->neighbors[i];
        if(findcollision(p1, p2, &tmin))
        {
            partner = p2;
            type = 0;
        }
    }
    createevent(tmin + time, p1, partner, type);
    p1->counter2 = partner->counter;
}



/**************************************************
**                FINDALLCOLLISIONS
** All collisions of all particle pairs
**************************************************/
void findallcollisions()      
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        particle* p1 = particles + i;
        particle* partner = p1;
        double tmin = findneighborlistupdate(p1);
        int type = 8;
        for (j = 0; j < p1->nneigh; j++)
        {
            particle* p2 = p1->neighbors[j];
            if (p2 > p1)
            {
                if(findcollision(p1, p2, &tmin))
                {
                    partner = p2;
                    type = 0;
                }
            }
        }
        createevent(tmin, p1, partner, type);
        p1->counter2 = partner->counter;
    }
}







/**************************************************
**                  COLLISION
** Process a single collision event
**************************************************/
void collision(particle* p1)
{
    particle* p2 = p1->p2;
    update(p1);
    if (p1->counter2 != p2->counter)    //Event already invalidated?
    {
        p1->counter++;
        findcollisions(p1);
        return;
    }

    update(p2);
    p1->counter++;
    p2->counter++;

    double m1 = p1->mass;
    double m2 = p2->mass;

    double dx = (p1->x - p2->x);		//Normalized distance vector
    double dy = (p1->y - p2->y);
    double dz = (p1->z - p2->z);
    if (p1->nearboxedge)
    {
        if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
        if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
        if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
    }
    double r = sqrt(dx*dx+dy*dy+dz*dz);
    double rinv = 1.0 / r;
    dx *= rinv;  dy *= rinv;  dz *= rinv;

    double dvx = p1->vx - p2->vx;                               //relative velocity
    double dvy = p1->vy - p2->vy;
    double dvz = p1->vz - p2->vz;

    double b = dx * dvx + dy * dvy + dz * dvz;                  //dr.dv
    b *= 2.0 / (m1 + m2);
    double dv1 = b * m2, dv2 = b * m1;

    p1->vx -= dv1 * dx;         //Change velocities after collision
    p1->vy -= dv1 * dy;         //delta v = (-) dx2.dv2
    p1->vz -= dv1 * dz;
    p2->vx += dv2 * dx;
    p2->vy += dv2 * dy;
    p2->vz += dv2 * dz;


    dvtot += b * r;             //Keep track of momentum transfer to calculate pressure
    colcounter++;

    removeevent(p2);


    findcollisions(p1);
    findcollisions(p2);
}






/**************************************************
**                 INITEVENTS
** Initialize event queue
**************************************************/
void initevents()
{
    eventlisttime = eventlisttimemultiplier / N;
    numeventlists = ceil(maxscheduletime / eventlisttime);
    maxscheduletime = numeventlists * eventlisttime;
    printf("number of event lists: %d\n", numeventlists);

    eventlists = (particle**)malloc((numeventlists + 1) * sizeof(particle*));
    if (!eventlists)
    {
        printf("Failed to allocate memory for eventlists\n");
        exit(3);
    }
    else memset(eventlists, 0, (numeventlists + 1) * sizeof(particle*));


    root = particles + N;				//Create root event for the binary tree
    root->eventtime = -99999999999.99;	//Root event is empty, but makes sure every other event has a parent
    root->eventtype = 127;				//This makes sure we don't have to keep checking this when adding/removing events
    root->parent = NULL;

    writeevent = particles + N + 1;		//Pick first unused event
    writeevent->eventtime = 0;
    writeevent->eventtype = 100;
    root->right = writeevent;
    writeevent->parent = root;
    writeevent->p2 = NULL;
    printf("Event tree initialized.\n");


    logsnapshotevent = particles + N + 4;
    logsnapshotevent->eventtype = 102;
    logsnapshotevent->p2 = NULL;
    //(First logsnap-event gets added in initlogsnapshots())

    resetevent = particles + N + 3;
    resetevent->eventtype = 50;
    resetevent->p2 = NULL;
    resetevent->eventtime = resetinterval;
    addevent(resetevent);
}

/**************************************************
**                  ADDEVENTTOTREE
**************************************************/
void addeventtotree(particle* newevent)
{
    double time = newevent->eventtime;
    particle* loc = root;
    int busy = 1;
    while (busy)						//Find location to add event into tree (loc)
    {
        if (time < loc->eventtime)				//Go left
        {
            if (loc->left) loc = loc->left;
            else
            {
                loc->left = newevent;
                busy = 0;
            }
        }
        else						//Go right
        {
            if (loc->right) loc = loc->right;
            else
            {
                loc->right = newevent;
                busy = 0;
            }
        }
    }
    newevent->parent = loc;
    newevent->left = NULL;
    newevent->right = NULL;

}

/**************************************************
**                  ADDEVENT
**************************************************/
void addevent(particle* newevent)
{
    double dt = newevent->eventtime - reftime;

    if (dt < eventlisttime) //Put it in the tree
    {
        newevent->queue = currentlist;
        addeventtotree(newevent);
    }
    else
    {
        int list_id = currentlist + dt / eventlisttime;
        if (list_id >= numeventlists)
        {
            list_id -= numeventlists;
            if (list_id > currentlist - 1) list_id = numeventlists; //Overflow
        }

        newevent->queue = list_id;
        newevent->right = eventlists[list_id]; //Add to linear list
        newevent->left = NULL;
        if (newevent->right) newevent->right->left = newevent;
        eventlists[list_id] = newevent;
    }
}
/**************************************************
**                  CREATEEVENT
**************************************************/
void createevent(double time, particle* p1, particle* p2, int type)
{
    p1->eventtime = time;
    p1->eventtype = type;
    p1->p2 = p2;
    addevent(p1);
}

/**************************************************
**                     ADDNEXTEVENTLIST
**************************************************/
void addnexteventlist()
{
    currentlist++;
    reftime += eventlisttime;
    if (currentlist == numeventlists) 
    {
        currentlist = 0;
        particle* ev = eventlists[numeventlists];//Overflow queue
        eventlists[numeventlists] = NULL;
        while (ev)
        {
            particle* nextev = ev->right;
            addevent(ev);
            ev = nextev;
            listcounter2 += 1;
        }            
    }

    //   printf("Currentlist is now %d (%lf)\n", currentlist, reftime);

    particle* ev = eventlists[currentlist];
    while (ev)
    {
        particle* nextev = ev->right;
        addeventtotree(ev);
        ev = nextev;
        listcounter1++;
    }
    eventlists[currentlist] = NULL;
    mergecounter++;
}

/**************************************************
**                  REMOVEEVENT
**************************************************/
void removeevent(particle* oldevent)
{
    if (oldevent->queue != currentlist)
    {
        if (oldevent->right) oldevent->right->left = oldevent->left;
        if (oldevent->left) oldevent->left->right = oldevent->right;
        else
        {
            eventlists[oldevent->queue] = oldevent->right;
        }
        return;
    }

    particle* parent = oldevent->parent;
    particle* node;					//This node will be attached to parent in the end


    if (oldevent->left == NULL)			//Only one child: easy to delete
    {
        node = oldevent->right;			//Child2 is attached to parent
        if (node)
        {
            node->parent = parent;
        }
    }
    else if (oldevent->right == NULL)		//Only one child again
    {
        node = oldevent->left;			//Child1 is attached to parent
        node->parent = parent;
    }
    else		  //Node to delete has 2 children
    {               //In this case: a) Find first node after oldevent     (This node will have no left)
                    //              b) Remove this node from the tree     (Attach node->right to node->parent)
                    //              c) Put this node in place of oldevent (Oldevent's children are adopted by node)
        node = oldevent->right;
        while (node->left) node = node->left;	//Find first node of right tree of descendants of oldevent
        particle* pnode = node->parent;
        if (pnode != oldevent)			//node is not a child of oldevent
        {						//Both of oldevent's children should be connected to node
            pnode->left = node->right;		//Remove node from right tree
            if (node->right) node->right->parent = pnode;
            oldevent->left->parent = node;
            node->left = oldevent->left;
            oldevent->right->parent = node;
            node->right = oldevent->right;
        }
        else					//This means node == oldevent->right
        {						//Only left has to be attached to node
            oldevent->left->parent = node;
            node->left = oldevent->left;
        }
        node->parent = parent;
    }
    if (parent->left == oldevent) parent->left = node;
    else                          parent->right = node;
}


/**************************************************
**                  OUTPUTSNAPSHOT
**************************************************/
void outputsnapshot()
{
    char filename[200];
    sprintf(filename, "snapshot_end.sph");
    FILE* file = fopen(filename, "w");
    int i;
    particle* p;
    fprintf(file, "%d\n%.12lf %.12lf %.12lf\n", (int)N, xsize, ysize, zsize);
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        double dt = time - p->t;
        p->x += p->vx * dt;
        p->y += p->vy * dt;
        p->z += p->vz * dt;
        p->t = time;


        fprintf(file, "%c %.12lf  %.12lf  %.12lf  %lf\n", 'a' + p->type, p->x + xsize * p->boxestraveledx, p->y + ysize * p->boxestraveledy, p->z + zsize * p->boxestraveledz, p->radius);
    }
    fclose(file);

}
/**************************************************
**                    GETEN
** Calculate potential energy of a particle
**************************************************/
double geten(particle* p1)
{
    int i;
    particle* p2;
    double en = 0;
    for (i = 0; i < p1->nneigh; i++)
    {
        p2 = p1->neighbors[i];
        double dx = p1->x - p2->x;
        double dy = p1->y - p2->y;
        double dz = p1->z - p2->z;
        if (p1->nearboxedge)
        {
            if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
            if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
            if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
        }
        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > maxrange2) continue;
        double sigmaij2 = p1->radius + p2->radius;
        sigmaij2 *= sigmaij2;
        en += potential(r2/sigmaij2);        
    }
    return en;
}

/**************************************************
**                    WRITE
** Writes information to the screen and makes snapshots
**************************************************/
void write()
{
    static int counter = 0;
    static int first = 1;
    static double lastsnapshottime = -999999999.9;
    static double timelast = 0;
    static double dvtotlast = 0;
    int i;
    particle* p;
    FILE* file;
    double realtime = time  + timeoffset;

    double en = 0, kinen = 0;
    int maxneigh = 0, minneigh = MAXNEIGH;
    for (i = 0; i < N; i++) update(particles+i);

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        kinen += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
        en += geten(p);
        if (p->nneigh > maxneigh) maxneigh = p->nneigh;
        if (p->nneigh < minneigh) minneigh = p->nneigh;
    }
    double temperature = 0.5 * kinen / (double)(N-1) / 1.5;
    en *= 0.5/N;

    double volume = xsize * ysize * zsize;
    double dens = N / volume;
    double pressid = dens;
    double pressnow = pressid -(dvtot - dvtotlast) / (3.0 * volume * (realtime - timelast));
    if (realtime == timelast) pressnow = 0;
    dvtotlast = dvtot;
    timelast = realtime;


    double listsize1 = (double)listcounter1 / mergecounter;
    double listsize2 = (double)listcounter2 / ((double) mergecounter / numeventlists);
    if (mergecounter == 0) listsize1 = listsize2 = 0.0;
    listcounter1 = listcounter2 = mergecounter = 0;

    printf("Simtime: %lf, Collisions: %u, Press: %lf, PotEn: %lf, Temperature: %lf, Listsizes: (%lf, %lf), Neigh: %d - %d\n", 
                realtime, colcounter, pressnow, en, temperature, listsize1, listsize2, minneigh, maxneigh);
    char filename[200];
    if (makesnapshots && realtime - lastsnapshottime > snapshotinterval - 0.001)
    {
        sprintf(filename, "mov.n%d.v%.4lf.sph", N, xsize * ysize * zsize);
        if (first) { first = 0; file = fopen(filename, "w"); }
        else                     file = fopen(filename, "a");
        fprintf(file, "%d\n%.12lf %.12lf %.12lf\n", (int)N, xsize, ysize, zsize);
        for (i = 0; i < N; i++)
        {
            p = &(particles[i]);
            update(p);

            fprintf(file, "%c %.12lf  %.12lf  %.12lf  %lf\n", 'a' + p->type, p->x + xsize * p->boxestraveledx, p->y + ysize * p->boxestraveledy, p->z + zsize * p->boxestraveledz, p->radius);
        }
        fclose(file);
        lastsnapshottime = realtime;
    }

    sprintf(filename, "press.n%d.v%.4lf.sph", N, xsize * ysize * zsize);
    if (counter == 0) file = fopen(filename, "w");
    else              file = fopen(filename, "a");
    fprintf(file, "%lf %lf %lf\n", realtime, pressnow, en);
    fclose(file);

    counter++;

    writeevent->eventtime = time + writeinterval;
    addevent(writeevent);
}



/**************************************************
**                    BACKINBOX
** Put particle back in the box.
** Just for initialization
**************************************************/
void backinbox(particle* p)
{
    p->x -= xsize * floor(p->x / xsize);
    p->y -= ysize * floor(p->y / ysize);
    p->z -= zsize * floor(p->z / zsize);
}



/**************************************************
**                 RANDOM_GAUSSIAN
** Generates a random number following a 
** normal distribution (mean 0, std. dev. 1)
**************************************************/
double random_gaussian()
{
    static int have_deviate = 0;
    static double u1, u2;
    double  x1, x2, w;

    if (have_deviate)
    {
        have_deviate = 0;
        return u2;
    }
    else
    {
        do
        {
            x1 = 2.0 * genrand_real2() - 1.0;
            x2 = 2.0 * genrand_real2() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);
        w = sqrt((-2.0 * log(w)) / w);
        u1 = x1 * w;
        u2 = x2 * w;
        have_deviate = 1;
        return u1;
    }
}



/**************************************************
**                COMPARE
** For sorting logarithmic snapshot times
**************************************************/
int compare (const void * a, const void * b)
{
    double d = ((logsnaptime*)a)->t - ((logsnaptime*)b)->t;

    return ( (0<d) - (d<0) );
}

/**************************************************
**                INITLOGSNAPSHOTS
**************************************************/
void initlogsnapshots()
{
  double shifttime = maxtime / numbatches;
  int i,j, counter= 0;
  int numshifts = maxtime / shifttime;
  double fac = pow(maxtime/minimumtime, 1.0/numsnapshots);
  double t,dt;
  for (i = 0; i < numshifts; i++)
  {
        if (counter == MAXLOGSNAPSHOTS - 1)
        {
            printf ("Warning: trying to make too many snapshots (reduce the number of snapshots to make, or increase MAXLOGSNAPSHOTS). (%d, %d, %d)\n",
                numshifts, numsnapshots, (int)MAXLOGSNAPSHOTS);
            break;
        }
        dt = minimumtime;
        for (j = 0; j <= numsnapshots; j++)
        {
            if (counter == MAXLOGSNAPSHOTS) break;
            if (j == 0)
            {
                t = i*shifttime;
            }
            else
            {
                t = i * shifttime + dt;
                dt*=fac;
            }
            if (t < maxtime)
            {
                logsnapshottimes[counter].t = t;
                logsnapshottimes[counter].batch = i;
                logsnapshottimes[counter].index = j;
                counter++;
            }
        }
  }


  qsort(logsnapshottimes, counter, sizeof(logsnaptime), compare);
  logsnapshottimes[counter].t   = maxtime - 0.000000000001;		//Last snapshot
  logsnapshottimes[counter].batch = i;		//Last snapshot
  logsnapshottimes[counter].index = 0;
  logsnapshottimes[counter+1].t = maxtime +1;			//Just to make sure that the next snapshot won't happen
  logsnapshottimes[counter+1].batch = i;		//Last snapshot
  logsnapshottimes[counter+1].index = 0;

  nextlogsnapshot = 0;
  logsnapshotevent->eventtime = logsnapshottimes[nextlogsnapshot].t;
  addevent(logsnapshotevent);


//   for (i = 0; i < counter; i++)
//   {
//     printf ("%d   %lf\n", i, logsnapshottimes[i]);
//   }
//   exit(3);

}


/**************************************************
**                    WRITELOGSNAPSHOT
** Writes a snapshot
**************************************************/
void writelogsnapshot()
{
  int batch = logsnapshottimes[nextlogsnapshot].batch;
  int index = logsnapshottimes[nextlogsnapshot].index;
  particle* p;
  int i;
  FILE* file;
  char filename[200];
  sprintf (filename, "logsnap.b%03d.s%05d", batch, index);
  file = fopen (filename, "w");
  if (!file)
  {
    printf ("Failed to open file %s\n", filename);
    return;
  }
  fprintf (file, "%d\n%lf %lf %lf    Time: %lf\n", (int) N, xsize, ysize, zsize, time+timeoffset);
  for ( i = 0; i < N; i++)
  {
    p = &(particles[i]);
    update(p);
    fprintf(file, "%c %.12lf  %.12lf  %.12lf  %lf\n", 'a' + p->type,   p->x + xsize * p->boxestraveledx,  p->y + ysize * p->boxestraveledy,  p->z + zsize * p->boxestraveledz, p->radius);
  }
  fclose(file);
  nextlogsnapshot++;
  logsnapshotevent->eventtime = logsnapshottimes[nextlogsnapshot].t - timeoffset;
  addevent (logsnapshotevent);     //Add next write interval

}


/**************************************************
**                     RESETTIME
** Resets simulation time to zero for additional
** precision.
**************************************************/
void resettime()
{
    int i;
    particle* p;

    for (i = 0; i < N; i++)     //Adjust time for all particles
    {
        p = particles + i;
        update(p);
        removeevent(p);
        p->t = 0;
        p->counter = 0;
        p->counter2 = 0;
    }
    time = 0;
    timeoffset += resetinterval;        //Keep track of the total offset
    removeevent(writeevent);                //Fix the write event timing
    if (makelogsnapshots) removeevent(logsnapshotevent);



    reftime = 0;
    currentlist = 0;

    for (i = 0; i < N; i++)
    {
        makeneighborlist(particles + i);        //This also reschedules all events
    }

    writeevent->eventtime -= resetinterval;
    addevent(writeevent);

    if(makelogsnapshots)
    {
        logsnapshotevent->eventtime -= resetinterval;
        addevent(logsnapshotevent);
    }




    //Schedule next reset event    
    resetevent->eventtime = resetinterval;
    addevent(resetevent);
}
