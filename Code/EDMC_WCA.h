// Event driven MD, headers.

typedef struct sparticle
{
	double x, y, z;
	double vx, vy, vz;
	double xn, yn, zn;
	struct sparticle* neighbors[MAXNEIGH];
	uint8_t nneigh;
	double t;
	double radius;
	double mass;
	uint8_t nearboxedge;
	int cell;
	int boxestraveledx, boxestraveledy, boxestraveledz;
	unsigned int counter;
	unsigned int counter2;
	struct sparticle* prev, * next;
	uint8_t type;
	double eventtime;
	struct sparticle* left;
	struct sparticle* right;
	struct sparticle* parent;
	struct sparticle* p2;
	int queue;
	unsigned char eventtype;
} particle;

typedef struct slogsnaptime
{
  double t;
  int batch;
  int index;
} logsnaptime;


int main();
void printstuff();
void init();


void initevents();
void fcc();
void simplecubic();
void loadparticles();
void randomparticles();
void randommovement();
void initcelllist();
void addtocelllist(particle* p, int cellx, int celly, int cellz);
int celloffset(int a, int b, int c);

void step();
int findcollision(particle*, particle*, double*);
void findallcollisions();
void findcollisions(particle*);
void collision(particle*);

void addevent(particle*);
void removeevent(particle*);
void createevent(double time, particle* p1, particle* p2, int type);
void addnexteventlist();
double findneighborlistupdate(particle* p1);
void makeneighborlist(particle* p1);

void outputsnapshot();
void write();
double random_gaussian();
void backinbox(particle* p);

double potential(double r2);
void resettime();
void initlogsnapshots();
void writelogsnapshot();

