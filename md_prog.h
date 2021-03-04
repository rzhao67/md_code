#include <stdio.h>

typedef double REAL;  // Real number
typedef struct {
	REAL x, y;
} VecR;  // 2D vector of doubles
typedef struct {
	VecR r, rv, ra;
} Mol;  // Position, velocity, and acceleration vectors associated with an atom or molecule
typedef struct {
    int x, y;
} VecI; // 2D vector of integers
typedef struct {
    REAL val, sum, sum2;
} Prop; // Property

typedef enum {N_I, N_R} VType;

typedef struct {
	char *vName;
	void *vPtr;
	VType vType;
	int vLen, vStatus;
} NameList;

/* Function definitions */
int GetNameList(int argc, char **argv);
void PrintNameList(FILE *fp);
void SetParams();
int SetupJob();
int AllocArrays();
void InitCoords();
void InitVels();
void InitAccels();
void InitRand(int randSeedI);
void VRand(VecR *p);
REAL RandR();
void AccumProps(int icode);
void SingleStep();
void LeapfrogStep(int part);
void ApplyBoundaryCond();
void ComputeForces();
void EvalProps();
void PrintSummary(FILE *fp);
void EvalVelDist();
void PrintVelDist(FILE *fp);