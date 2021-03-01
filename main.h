typedef double REAL;  // Real number
typedef struct {
	REAL x, y;
} VecR;  // 2D vector of doubles
typedef struct {
	VecR r, rv, ra;
} Mol;  // Position, velocity, and acceleration vectors associated with an atom or molecule
// typedef struct {
//     int x, y;
// } VecI;
// typedef struct {
//     REAL val, sum, sum2;
// } Prop;