#include "md_prog.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <sys/time.h>

/* Vector operations */
#define VAdd(v1, v2, v3) \
	(v1).x = (v2).x + (v3).x, \
    (v1).y = (v2).y + (v3).y
#define VSub(v1, v2, v3) \
    (v1).x = (v2).x - (v3).x, \
    (v1).y = (v2).y - (v3).y
#define VDot(v1, v2) \
    ((v1).x * (v2).x + (v1).y * (v2).y)
#define VSAdd(v1, v2, s3, v3) \
    (v1).x = (v2).x + (s3) * (v3).x, \
    (v1).y = (v2).y + (s3) * (v3).y
#define VSet(v, sx, sy) \
    (v).x = sx, \
    (v).y = sy
#define VSetAll(v, s) VSet(v, s, s)
#define VZero(v) VSetAll(v, 0)
#define VVSAdd(v1, s2, v2) VSAdd(v1, v1, s2, v2)
#define VLen(v) sqrt(VDot(v, v))
#define VLenSq(v) VDot(v, v)
#define VMul(v1, v2, v3) \
    (v1).x = (v2).x * (v3).x, \
    (v1).y = (v2).y * (v3).y
#define VDiv(v1, v2, v3) \
    (v1).x = (v2).x / (v3).x, \
    (v1).y = (v2).y / (v3).y
#define VScale(v, s) \
    (v).x *= s, \
    (v).y *= s
#define VVAdd(v1, v2) VAdd(v1, v1, v2)
#define VSCopy(v2, s1, v1) \
    (v2).x = (s1) * (v1).x, \
    (v2).y = (s1) * (v1).y
#define VProd(v) ((v).x * (v).y)
#define VCSum(v) \
    ((v).x + (v).y) // + (v).z) for 3D

/* Periodic boundary functions */
#define VWrap(v, t) \
    if (v.t >= 0.5 * region.t) v.t -= region.t; \
    else if (v.t < -0.5 * region.t) v.t += region.t
#define VWrapAll(v) \
    {VWrap(v, x); \
     VWrap(v, y);}

/* Prop functions */
#define PropZero(v) \
    v.sum = 0., \
    v.sum2 = 0.
#define PropAccum(v) \
    v.sum += v.val, \
    v.sum2 += Sqr(v.val)
#define PropAvg(v, n) \
    v.sum /= n, \
    v.sum2 = sqrt(Max(v.sum2 / n - Sqr(v.sum), 0.))
#define PropEst(v) \
    v.sum, v.sum2

/* Input processing functions */
#define NameI(x) {#x, &x, N_I, sizeof(x) / sizeof(int)}
#define NameR(x) {#x, &x, N_R, sizeof(x) / sizeof(REAL)}
#define NP_I ((int *) (nameList[k].vPtr) + j)
#define NP_R ((REAL *) (nameList[k].vPtr) + j)

/* Constants */
/* Random number generation */
#define IADD 453806245
#define IMUL 314159269
#define MASK 2147483647
#define SCALE 0.4656612873e-9
int randSeedP = 17;
/* Other constants */
#define NDIM 2

/* Other macros */
#define Sqr(x) ((x) * (x))
#define Cube(x) ((x) * (x) * (x))
#define Max(x1, x2) \
    (((x1) > (x2)) ? (x1) : (x2))
#define Min(x1, x2) \
    (((x1) < (x2)) ? (x1) : (x2))
#define DO_MOL for (n = 0; n < nMol; n++)
#define AllocMem(a, n, t) a = (t *) malloc((n) * sizeof(t))
// #define NameVal(x) \
//     if (!strncmp(bp, #x, strlen(#x))) { \
//         bp += strlen(#x); \
//         x = strtod(bp, &bp); \
//     }

/* Initialize variables */
Mol *mol;
VecR region, vSum;
VecI initUcell;
Prop kinEnergy, pressure, totEnergy;
REAL deltaT, density, temperature, rCut, velMag, timeNow, uSum, virSum, vvSum, *histVel, rangeVel;
int stepAvg, stepEquil, stepLimit, nMol, stepCount, moreCycles, stepVel, countVel, limitVel, sizeHistVel;
REAL hFunction;
int randSeed;

NameList nameList[] = {
    NameR(deltaT),
    NameR(density),
    NameI(initUcell),
    NameI(stepAvg),
    NameI(stepEquil),
    NameI(stepLimit),
    NameR(temperature),
    NameI(limitVel),
    NameR(rangeVel),
    NameI(sizeHistVel),
    NameI(stepVel),
    NameI(randSeed),
};

int GetNameList(int argc, char **argv) {
    int id, j, k, match, ok;
    char buff[80], *token;
    FILE *fp;

    strcpy(buff, argv[0]);  // Copies the executable name into buff
    strcat(buff, ".in");
    ok = -1;
    fp = fopen(buff, "r");
    if (fp == NULL) {
        return ok;
    }
    for (k = 0; k < sizeof(nameList) / sizeof(NameList); k++)
        nameList[k].vStatus = 0;  // Initialize nameList
    while (1) {
        fgets(buff, 80, fp);
        if (feof(fp)) break;
        token = strtok(buff, " \t\n");
        if (!token) break;
        match = 0;
        for (k = 0; k < sizeof(nameList) / sizeof(NameList); k++) {
            if (strcmp(token, nameList[k].vName) == 0) {
                match = 1;
                if (nameList[k].vStatus == 0) {
                    nameList[k].vStatus = 1;
                    for (j = 0; j < nameList[k].vLen; j++) {
                        token = strtok(NULL, ", \t\n");
                        if (token) {
                            switch (nameList[k].vType) {
                                case N_I:
                                    *NP_I = atol(token);
                                    break;
                                case N_R:
                                    *NP_R = atof(token);
                                    break;
                            }
                        } else {
                            nameList[k].vStatus = 2;
                            ok = 0;
                        }
                    }
                    token = strtok(NULL, ", \t\n");
                    if (token) {
                        nameList[k].vStatus = 3;
                        ok = 0;
                    }
                    break;
                } else {
                    nameList[k].vStatus = 4;
                    ok = 0;
                }
            }
        }
        if (!match) ok = 0;
    }
    fclose(fp);
    for (k = 0; k < sizeof(nameList) / sizeof(NameList); k++) {
        if (nameList[k].vStatus != 1) ok = 0;
    }
    return ok;
}

void PrintNameList(FILE *fp) {
    int j, k;

    fprintf(fp, "NameList -- data\n");
    for (k = 0; k < sizeof(nameList) / sizeof(NameList); k++) {
        fprintf(fp, "%s\t", nameList[k].vName);
        if (strlen(nameList[k].vName) < 8) fprintf(fp, "\t");
        if (nameList[k].vStatus > 0) {
            for (j = 0; j < nameList[k].vLen; j++) {
                switch (nameList[k].vType) {
                    case N_I:
                        fprintf(fp, "%d ", *NP_I);
                        break;
                    case N_R:
                        fprintf(fp, "%#g ", *NP_R);
                        break;
                }
            }
        }
        switch (nameList[k].vStatus) {
            case 0:
                fprintf(fp, "** no data");
                break;
            case 1:
                break;
            case 2:
                fprintf(fp, "** missing data");
                break;
            case 3:
                fprintf(fp, "** extra data");
                break;
            case 4:
                fprintf(fp, "** multiply defined");
                break;
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "----\n");
}

void SetParams() {
    rCut = pow(2., 1./6.);
    VSCopy(region, 1./sqrt(density), initUcell);
    nMol = VProd(initUcell);
    velMag = sqrt(NDIM * (1. - 1./nMol) * temperature);
    printf("rCut:%f\nregion_x:%f\nregion_y:%f\nnMol:%i\nvelMag:%f\n", rCut, region.x, region.y, nMol, velMag);
    printf("----\n");
}

int SetupJob() {
    printf("Setting up job...\n");
    int failure = AllocArrays();
    if (failure) {
        printf("Failed to allocate memory.\n");
        return -1;
    }
    stepCount = 0;
    InitCoords();
    InitVels();
    InitAccels();
    printf("Initializing total energy, kinetic energy, and pressure...\n");
    AccumProps(0);
    countVel = 0;
    InitRand(randSeed);
    printf("Setup complete!\n");
    printf("----\n");
    return 0;
}

int AllocArrays() {
    printf("Allocating memory...\n");
    AllocMem(mol, nMol, Mol);
    if (mol == NULL) {
        return -1;
    } else {
        printf("Successfully allocated memory to molecules!\n");
    }
    AllocMem(histVel, sizeHistVel, REAL);
    return 0;
}

void InitCoords() {
    VecR c, gap;
    int n, nx, ny;

    printf("Initializing coordinates...\n");
    VDiv(gap, region, initUcell);
    n = 0;
    for (ny = 0; ny < initUcell.y; ny++) {
        for (nx = 0; nx < initUcell.x; nx++) {
            VSet(c, nx + 0.5, ny + 0.5);
            VMul(c, c, gap);
            VVSAdd(c, -0.5, region);
            mol[n].r = c;
            ++n;
        }
    }
}

void InitVels() {
    int n;

    printf("Initializing velocities...\n");
    VZero(vSum);
    DO_MOL {
        VRand(&mol[n].rv);
        VScale(mol[n].rv, velMag);
        VVAdd(vSum, mol[n].rv);
    }
    DO_MOL VVSAdd(mol[n].rv, -1. / nMol, vSum);
}

void InitAccels() {
    int n;

    printf("Initializing accelerations...\n");
    DO_MOL VZero(mol[n].ra);
}

void InitRand(int randSeedI) {
    struct timeval tv;

    if (randSeedI != 0) randSeedP = randSeedI;
    else {
        gettimeofday(&tv, 0);
        randSeedP = tv.tv_usec;
    }
}

// void VRand(VecR *p) { // For 3D vectors
//     REAL s, x, y;

//     s = 2.;
//     while (s > 1.) {
//         x = 2. * RandR() - 1.;
//         y = 2. * RandR() - 1.;
//         s = Sqr(x) + Sqr(y);
//     }
//     p->z = 1. - 2. * s;
//     s = 2. * sqrt(1. - s);
//     p->x = s * x;
//     p->y = s * y;
// }

void VRand(VecR *p) {
    REAL s;

    s = 2. * M_PI * RandR();
    p->x = cos(s);
    p->y = sin(s);
}

REAL RandR() {
    randSeedP = (randSeedP * IMUL + IADD) & MASK;
    return randSeedP * SCALE;
}

void AccumProps(int icode) {
    if (icode == 0) {
        // printf("Initializing total energy, kinetic energy, and pressure...\n");
        PropZero(totEnergy);
        PropZero(kinEnergy);
        PropZero(pressure);
    } else if (icode == 1) {
        PropAccum(totEnergy);
        PropAccum(kinEnergy);
        PropAccum(pressure);
    } else if (icode == 2) {
        PropAvg(totEnergy, stepAvg);
        PropAvg(kinEnergy, stepAvg);
        PropAvg(pressure, stepAvg);
    }
}

void SingleStep() {
	++stepCount;
	timeNow = stepCount * deltaT;
	LeapfrogStep(1);
	ApplyBoundaryCond();
	ComputeForces();
	LeapfrogStep(2);
	EvalProps();
	AccumProps(1);
	if (stepCount % stepAvg == 0) {
		AccumProps(2);
		PrintSummary(stdout);
		AccumProps(0);
	}
    // if (stepCount >= stepEquil &&
    //     (stepCount - stepEquil) % stepVel == 0) EvalVelDist();
}

void LeapfrogStep(int part) {
    int n;

    if (part == 1) {
        DO_MOL {
            VVSAdd(mol[n].rv, 0.5 * deltaT, mol[n].ra);
            VVSAdd(mol[n].r, deltaT, mol[n].rv);
        }
    } else {
        DO_MOL VVSAdd(mol[n].rv, 0.5 * deltaT, mol[n].ra);
    }
}

void ApplyBoundaryCond() {
    int n;

    DO_MOL VWrapAll(mol[n].r);
}

void ComputeForces() {
	VecR dr;
	REAL fcVal, rr, rrCut, rri, rri3;
	int j1, j2, n;

	rrCut = Sqr(rCut);
	DO_MOL VZero(mol[n].ra);
	uSum = 0.;
    virSum = 0.;
	for (j1 = 0; j1 < nMol - 1; j1++) {
		for (j2 = j1 + 1; j2 < nMol; j2++) {
			VSub(dr, mol[j1].r, mol[j2].r);
            VWrapAll(dr);
            rr = VLenSq(dr);
            if (rr < rrCut) {
                rri = 1./ rr;
                rri3 = Cube(rri);
                fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
                VVSAdd(mol[j1].ra, fcVal, dr);
                VVSAdd(mol[j2].ra, -fcVal, dr);
                uSum += 4. * rri3 * (rri3 - 1.) + 1.;
                virSum += fcVal * rr;
			}
		}
	}
}

void EvalProps() {
    REAL vv;
    int n;

    VZero(vSum);
    vvSum = 0.;
    DO_MOL {
        VVAdd(vSum, mol[n].rv);
        vv = VLenSq(mol[n].rv);
        vvSum += vv;
        kinEnergy.val = 0.5 * vvSum / nMol;
        totEnergy.val = kinEnergy.val + uSum / nMol;
        pressure.val = density * (vvSum + virSum) / (nMol * NDIM);
    }
}

void PrintSummary(FILE *fp) {
    fprintf(fp,
        "%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
        stepCount, timeNow, VCSum(vSum) / nMol, PropEst(totEnergy),
        PropEst(kinEnergy), PropEst(pressure));
}

void EvalVelDist() {
    REAL deltaV, histSum;
    int j, n;

    if (countVel == 0) {
        for (j = 0; j < sizeHistVel; j++) histVel[j] = 0.;
    }
    deltaV = rangeVel / sizeHistVel;
    DO_MOL {
        j = VLen(mol[n].rv) / deltaV;
        ++histVel[Min(j, sizeHistVel - 1)];
    }
    ++countVel;
    if (countVel == limitVel) {
        histSum = 0.;
        for (j = 0; j < sizeHistVel; j++) histSum += histVel[j];
        for (j = 0; j < sizeHistVel; j++) histVel[j] /= histSum;
        /* This part of code may be in the wrong spot */
        hFunction = 0.;
        for (j = 0; j < sizeHistVel; j++) {
            if (histVel[j] > 0.) hFunction += histVel[j] * log(histVel[j] / ((j + 0.5) * deltaV));
        }
        /* The above piece of code may be in the wrong spot */
        PrintVelDist(stdout);
        countVel = 0;
    }

}

void PrintVelDist(FILE *fp) {
    REAL vBin;
    int n;

    printf("vdist (%.3f)\n", timeNow);
    for (n = 0; n < sizeHistVel; n++) {
        vBin = (n + 0.5) * rangeVel / sizeHistVel;
        fprintf(fp, "%8.3f %7.3f\n", vBin, histVel[n]);
    }
    fprintf(fp, "hfun: %8.3f %8.3f\n", timeNow, hFunction);
}

int main(int argc, char **argv) {
	int failure = GetNameList(argc, argv);
    if (!failure) {
        printf("Error in input file.\n");
        return -1;
    }
	PrintNameList(stdout);
	SetParams();
	failure = SetupJob();
    if (failure) {
        return -1;
    }
	moreCycles = 1;
	while (moreCycles) {
		SingleStep();
		if (stepCount >= stepLimit) moreCycles = 0;
	}
}
