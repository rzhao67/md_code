#include <stdio.h>

#define Vadd(v1, v2, v3) \
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
#define VLenSq(v) VDot(v, v)
#define VWrap(v, t) \
    if (v.t >= 0.5 * region.t) v.t -= region.t; \
    else if (v.t < -0.5 * region.t) v.t += region.t
#define VWrapAll(v) \
    {VWrap(v, x); \
     VWrap(v, y);}
#define Sqr(x) ((x) * (x))
#define Cube(x) ((x) * (x) * (x))
#define DO_MOL for (n = 0; n < nMol; n++)
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
#define NDIM 2
#define AllocMem(a, n, t) a = (t *) malloc((n) * sizeof(t))
#define VSCopy(v2, s1, v1) \
    (v2).x = (s1) * (v1).x, \
    (v2).y = (s1) * (v1).y
#define VProd(v) ((v).x * (v).y)
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
#define VLen(v) sqrt(VDot(v, v))

typedef double real;
typedef struct {
	real x, y;
} VecR;
typedef struct {
	VecR r, rv, ra;
} Mol;
typedef struct {
    int x, y;
} VecI;
typedef struct {
    real val, sum, sum2;
} Prop;

Mol *mol;
VecR region, vSum;
VecI initUCell;
Prop kinEnergy, pressure, totEnergy;
real deltaT, density, rCut, temperature, timeNow, uSum, velMag, virSum, vvSum;
int moreCycles, nMol, stepAvg, stepCount, stepEquil, stepLimit;
real *histVel, rangeVel;
int countVel, limitVel, sizeHistVel, stepVel;
real hFunction;

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
    if (stepCount >= stepEquil &&
        (stepCount - stepEquil) % stepVel == 0) EvalVelDist();
}

void SetupJob() {
	AllocArrays();
	stepCount = 0;
	InitCoords();
	InitVels();
	InitAccels();
	AccumProps(0);
    countVel = 0;
    InitRand(randSeed);
}

void ComputeForces() {
	VecR dr;
	real fcVal, rr, rrCut, rri, rri3;
	int j1, j2, n;

	rrCut = Sqr(rcut);
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

void InitCoords() {
    VecR c, gap;
    int n, nx, ny;

    VDiv(gap, region, unitUcell);
    n = 0;
    for (ny = 0; ny < initUCell.y; ny++) {
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

    DO_MOL VZero(mol[n].ra);
}

void AllocArrays() {
    AllocMem(mol, nMol, Mol);
    AllocMem(histVel, sizeHistVel, real);
}

void SetParams() {
    rCut = pow(2., 1./6.);
    VSCopy(region, 1./sqrt(density), initUcell);
    nMol = VProd(initUcell);
    velMag = sqrt(NDIM * (1. - 1./nMol) * temperature);
}

void EvalProps() {
    real vv;
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

void AccumProps(int icode) {
    if (icode == 0) {
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

void PrintSummary(FILE *fp) {
    fprintf(fp,
        "%5d %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
        stepCount, timeNow, VCSum(vSum) / nMol, PropEst(totEnergy),
        PropEst(kinEnergy), PropEst(pressure));
}

void EvalVelDist() {
    real deltaV, histSum;
    int j, n;

    if (countVel == 0) {
        for (j = 0; j < sizeHistVel; j++) histVel[j] = 0.;
    }
    deltaV = rangeVel / sizeHistVel;
    DO_MOL {
        j = VLen(mol[n].rv) / deltaV;
        ++histVel[Min(j, sizeHistVel - i)];
    }
    ++countVel;
    if (countVel == limitVel) {
        histSum = 0.;
        for (j = 0; j < sizeHistVel; j++) histSum += histVel[j];
        for (j = 0; j < sizeHistVel; j++) histVel[j] /= histSum;
        PrintVelDist(stdout);
        countVel = 0;
    }
}

void PrintVelDist(FILE *fp) {
    real vBin;
    int n;

    printf("vdist (%.3f)\n", timeNow);
    for (n = 0; n < sizeHistVel; n++) {
        vBin = (n + 0.5) * rangeVel / sizeHistVel;
        fprintf(fp, "%8.3f %7.3f\n", vBin, histVel[n]);
    }
}

int main(int argc, char **argv) {
	GetNameList(argc, argv);
	PrintNameList(stdout);
	SetParams();
	SetupJob();
	moreCycles = 1;
	while (moreCycles) {
		SingleStep();
		if (stepCount >= stepLimit) moreCycles = 0;
	}
}