#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Redefine printf to only be on when compiled with DEBUG defined
#ifdef DEBUG
#define dbprintf(...) printf(__VA_ARGS__)
#else
#define dbprintf(...) { int no_op; }
#endif

// Multinomial choose(p,q) - probably a better way to do this in C
inline int polychoose(int* p, int* q)
{
  int ret = 1;
  // tgamma(p+1) = p factorial
  for (int d=0; d < 3; d++)
    ret *= tgamma(p[d]+1) / (tgamma(q[d]+1)*tgamma(p[d]-q[d]+1));
  return ret;
}

#if 0
void testPchooseQ()
{
  int p[3], q[3];
  printf("Test p choose q:");

  p[0] = 0; p[1] = 0; p[2] = 0;
  q[0] = 0; q[1] = 0; q[2] = 0; 
  int val = polychoose(p,q);
  printf("  (%d, %d, %d) choose (%d, %d, %d) = %d \n",
      p[0],p[1],p[2],q[0],q[1],q[2],val);

  p[0] = 1; p[1] = 0; p[2] = 0;
  q[0] = 0; q[1] = 0; q[2] = 0; 
  val = polychoose(p,q);
  printf("  (%d, %d, %d) choose (%d, %d, %d) = %d \n",
      p[0],p[1],p[2],q[0],q[1],q[2],val);

  p[0] = 1; p[1] = 1; p[2] = 2;
  q[0] = 0; q[1] = 1; q[2] = 1; 
  val = polychoose(p,q);
  printf("  (%d, %d, %d) choose (%d, %d, %d) = %d \n",
      p[0],p[1],p[2],q[0],q[1],q[2],val);

  p[0] = 1; p[1] = 1; p[2] = 4;
  q[0] = 0; q[1] = 1; q[2] = 2; 
  val = polychoose(p,q);
  printf("  (%d, %d, %d) choose (%d, %d, %d) = %d \n",
      p[0],p[1],p[2],q[0],q[1],q[2],val);

  p[0] = 3; p[1] = 2; p[2] = 4;
  q[0] = 2; q[1] = 1; q[2] = 2; 
  val = polychoose(p,q);
  printf("  (%d, %d, %d) choose (%d, %d, %d) = %d \n",
      p[0],p[1],p[2],q[0],q[1],q[2],val);
}
#endif

// Create matrices A for Pth-order interpolation in 3D
// --> A * x \approx b , A is Nrow x Ncol, Nrow > Ncol
// inputs:
// radius - determines order of polynomial, values 1-5 are reasonable in 3D
// batch - number of matrices to generate
// outputs:
// nrow, ncol - size of each matrix
// outmat - pointer to memory for matrix entries, col-major, batch index:
//    row1, col1, matrix1
//    row2, col1, matrix1
//    ...
//    rowN-1, colN, matrix1
//    rowN, colN, matrix1
//    row1, col1, matrix2
//    row2, col1, matrix2
//    ...
int genmat(int radius, int batch, int* outrow, int* outcol, double** outmat)
{
    int P = radius; // order of polynomial interpolation
    int N = 1 + (P*(8 + P*(6 + P*4)))/3; // star-shaped interpolation points
    // Size of the matrix, Nrow x Ncol, Nrow > Ncol for WLS
    int Nrow = N; 
    *outrow = Nrow;
    int Ncol = (P+1)*(P+2)*(P+3)/6;
    *outcol = Ncol;
    dbprintf("Radius = %d --> (Nrow, Ncol) = (%d, %d)\n\n", radius, Nrow, Ncol);
    int Nmat = batch; // number of matrices to generate
    dbprintf("batch size Nmat = %d \n", Nmat);
    // Allocate memory for return variable
    double* mat = (double*) malloc(Nmat * Nrow * Ncol * sizeof(double));
    *outmat = mat;
    
    // Determine what size problem you'll need to generate the batch
    int Nx = (2*P+1); // size of each cube we'll loop over
    int Npts = Nx*Nx*Nx; // total points of each cube we'll loop over
    dbprintf("size of cube Nx = %d \n", Nx);
    int Ncube = 1+Nmat/Npts; // number of cubes
    int Ntot = Npts*Ncube; // number of total points to generate moments for
    dbprintf("num cubes Ncube = %d, Npts = %d, Ntot = %d \n", Ncube, Npts, Ntot);

    int row, col, i, j, k, x, y, z;
    double x0, y0, z0;
    int* powers = (int*) malloc(3 * Ncol * sizeof(int));
    double* mom = (double*) malloc(Ntot * Ncol * sizeof(double));

    // set the power values in 2D array, [Ncol][3]
    int** pvals = (int**) malloc(Ncol * sizeof(int*));
    int ix=0;
    for (int pk=0; pk <= P; pk++)
    for (int pj=0; pj <= P; pj++)
    for (int pi=0; pi <= P; pi++)
    {
      if (pi+pj+pk <= P) // ignore otherwise
      {
        pvals[ix] = (powers+3*ix);
        pvals[ix][0]=pi;
        pvals[ix][1]=pj;
        pvals[ix][2]=pk;
        ix++;
      }
    }

    dbprintf("Power list:\n");
    for (ix=0; ix < Ncol; ix++)
      dbprintf("  pvals[%d]=(%d,%d,%d)\n", ix, 
        pvals[ix][0], pvals[ix][1], pvals[ix][2]);
    dbprintf("\n");
  
    // Moment indexes for centroids, treated differently
    int ploc[3] = {1, P+1, (P+1)*(P+2)/2};

    // Loop through the pts, creating random center and calculating moments
    double** mvals = (double**) malloc(Ntot * sizeof(double*));
    for (int ip=0; ip < Ntot; ip++)
    {
      mvals[ip] = (mom+Ncol*ip);
      int ix = ip % Npts; // loop around after each unit of cube Npts
      i = ix % Nx;
      j = ix/Nx % Nx;
      k = ix/(Nx*Nx);

      //  Calculate random centroid location in [0,1]
      double xc[3], dx[3];
      for (int d=0; d < 3; d++)
      {
        xc[d] = (double) rand() / (double) RAND_MAX;
        // For debugging use this line instead for regular cell moments
        // xc[d] = 0.5;
        //  Assume dx is the min of distance to nearest grid line
        dx[d] = 2.*fabs( (xc[d] > (1-xc[d])) ? 1-xc[d] : xc[d] );
      }
      xc[0] += i; // shift centroid to this index
      xc[1] += j; // y centroid
      xc[2] += k; // z centroid
      // Copy into their correct moment indexes
      mvals[ip][0] = 1.;
      mvals[ip][ploc[0]] = xc[0]; // x centroid
      mvals[ip][ploc[1]] = xc[1]; // y centroid
      mvals[ip][ploc[2]] = xc[2]; // z centroid
      // Fill in the others - relative to centroid
      for (int ic=2; ic < Ncol; ic++) // skip first 2 moments we have set
      {
        if (ic==ploc[1] || ic==ploc[2]) // skip others we have set
          continue;
        // Set all the other moments as if a slab dimensions dx[d]
        double tmp = 1;
        int even = 1;
        for (int d=0; d < 3; d++)
        {
          int pd = pvals[ic][d];
          even = even && !(pd % 2);
          tmp *= pow(dx[d]/2., pd) / (double) (pd+1);
        }
        mvals[ip][ic] = (even) ? tmp : 0; // Any odd moments are 0 about xc
      }
      double tmp;
      tmp = xc[0]+1;
      dbprintf("  xc[%4d]=( %1.2e, %1.2e, %1.2e)\n", ip, xc[0], xc[1], xc[2]);
      dbprintf("  dx[%4d]=( %1.2e, %1.2e, %1.2e)\n", ip, dx[0], dx[1], dx[2]);
      dbprintf("  mom[%4d] for (%d,%d,%d) = ", ip, i, j, k);
      for (int ic=0; ic < Ncol; ic++)
        dbprintf(" %1.2e, ", mvals[ip][ic]);
      dbprintf("\n");
    }
    dbprintf("\n");

    // For the star-shaped stencil, we'll need offsets
    int** offset = (int**) malloc(Nrow * sizeof(int*));
    int* offvals = (int*) malloc(3*Nrow * sizeof(int));
    ix=0;
    int offsetix = 0;
    dbprintf("Calculating stencil offsets:\n");
    for (int k=-P; k <= P; k++)
    for (int j=-P; j <= P; j++)
    for (int i=-P; i <= P; i++)
    {
      if (abs(i)+abs(j)+abs(k) <= P) // ignore otherwise
      {
        offset[ix] = (offvals+3*ix);
        offset[ix][0]=i;
        offset[ix][1]=j;
        offset[ix][2]=k;
        if (i==0 && j==0 && k==0)
          offsetix = ix;
        ix++;
      }
    }
    for (int ix=0; ix < Nrow; ix++)
      dbprintf("  offset[%d]=(%d,%d,%d)\n", ix, 
        offset[ix][0], offset[ix][1], offset[ix][2]);
    dbprintf("\n");

    // Build up the moment matrices
    double** avals = (double**) malloc(Nmat * sizeof(double*));
    for (int ip=0; ip < Nmat; ip++)
    {
      avals[ip] = (mat+Ncol*Nrow*ip);
      int ix = ip % Npts; // loop around after each unit of cube Npts
      int ic = ip / Npts; // which "cube" we're on
      // Indices of this point location
      i = ix % Nx;
      j = ix/Nx % Nx;
      k = ix/(Nx*Nx);

      // Matrix entries are shifted to this cell's xc
      double xc0[3], dx0[3];
      for (int d=0; d < 3; d++)
        xc0[d] = mvals[ip][ploc[d]]; // save centroid to shift to

      // Copy other cell's moments into the matrix
      dbprintf("For matrix[%4d]:\n", ip);
      for (int r=0; r < Nrow; r++)
      {
        if (r == offsetix) // Do special thing for this cell's moments
        {
          // Copy center cell's moments into the matrix
          for (int ic=0; ic < Ncol; ic++)
            avals[ip][r + ic*Nrow] = mvals[ip][ic];
          // Set the first moments (centroids) to be 0 
          for (int d=0; d < 3; d++)
            avals[ip][r + ploc[d]*Nrow] = 0; 
          dbprintf("  row[%d] (%d,%d,%d) center\n",r,i,j,k);
          continue; 
        }

        // 1. Loop over the offset
        int ir = offset[r][0] + i;
        int jr = offset[r][1] + j;
        int kr = offset[r][2] + k;

        // 2. If any of the offset indicess are <0 or >Nx, zero-fill
        if ((ir < 0 || ir >= Nx) ||
            (jr < 0 || jr >= Nx) ||
            (kr < 0 || kr >= Nx))
        {
          for (int ic=0; ic < Ncol; ic++)
            avals[ip][r + ic*Nrow] = 0;
          dbprintf("  row[%d] (%d,%d,%d) offset, =0\n",r,ir,jr,kr);
          continue;
        }

        // 3. Otherwise, calculate its index and copy moments into
        int ixr = ic*Npts + ir + Nx*(jr + kr*Nx);

        // Calculate the shift then set original centroids to 0
        for (int d=0; d < 3; d++)
          dx0[d] = mvals[ixr][ploc[d]] - xc0[d]; // d-th centroid

        dbprintf("  row[%d] index = %d, offset=(%1.2e,%1.2e,%1.2e)\n",
            r,ixr,dx0[0],dx0[1],dx0[2]);

        // Calculate the shift-based weight, dist^(-P-1)
        double dist = sqrt(dx0[0]*dx0[0]+dx0[1]*dx0[1]+dx0[2]*dx0[2]);
        double weight = pow(dist,-P-1);

        // 4. shift it to the current centroid using dx0 = old - new shift
        for (int ic=0; ic < Ncol; ic++)
        {
          avals[ip][r + ic*Nrow] = 0; // zero to accumulate

          for (int is=0; is < Ncol; is++) // shifted location
          {
            int* p = pvals[ic]; // new moment
            int* q = pvals[is]; // shifted moment
            if ((q[0] > p[0]) || (q[1] > p[1]) || (q[2] > p[2]))
              continue; // skip if not a valid q

            int pCq = polychoose(p,q);
            /*
            dbprintf("  (%d, %d, %d) choose (%d, %d, %d) = %d \n",
                p[0],p[1],p[2],q[0],q[1],q[2],pCq);
            */
            double shift = pCq*pow(dx0[0],p[0]-q[0])
                              *pow(dx0[1],p[1]-q[1]) 
                              *pow(dx0[2],p[2]-q[2]);
            // Fix for centroid entries which should be 0
            double momq = ((q[0]+q[1]+q[2])==1) ? 0 : mvals[ixr][is];
            // accumulate the shifted moment info
            avals[ip][r + ic*Nrow] += shift*momq;
          }
        }

        dbprintf("  row[%d] index = %d, shifted moments=\n",r,ixr);
        for (int ic=0; ic < Ncol; ic++)
          dbprintf("    %1.2e, ", avals[ip][r + ic*Nrow]);
        dbprintf("\n");

        // multiply by the weights
        for (int ic=0; ic < Ncol; ic++)
          avals[ip][r + ic*Nrow] *= weight;
      }
    }
    // Free local memory
    free(avals);
    free(offset);
    free(offvals);
    free(mvals);
    free(pvals);
    free(mom);
    free(powers);

    return 0;
}


int usage(int argc, char* argv[]) 
{
    printf("Usage: \"%s radius batch seed\"\n\n", argv[0]);
    printf("    radius - int radius of the stencil, [1,5]\n" 
           "    batch - int, number of matrices generated\n"
           "    seed - random number seed for repeatability \n\n");
    printf("This generates rank-deficient test matrices for \n" 
           "weighted least squares polynomial interpolation. \n"
           "Output is a # header, # line separated matrices, \n" 
           "with each matrix in col-major (row index first) ordering.\n");
    return 1; // exit
}


int main(int argc, char* argv[])
{
    if (argc != 4)
      return usage(argc, argv);

    int radius, batch, type, seed;
    // radius = 1-5 are good values to vary the size from 27x10 to 165x84
    radius = atoi(argv[1]);
    if ((radius > 6) || (radius < 1))
      return usage(argc, argv);

    // batch > 0
    batch = atoi(argv[2]);
    if (batch <= 0)
      return usage(argc, argv);

    // seed = for srand
    seed = abs(atoi(argv[3]));
    srand(seed);
    type = 1; // maybe provide more options later

    // Call the matrix generation routine
    double* mat = NULL;
    int row, col;
    genmat(radius, batch, &row, &col, &mat);

    // Write the output
    printf("# Matrix output:\n");
    printf("# batch=%d, Nrow=%d, Ncol=%d\n", batch, row, col);
    for (int b = 0; b < batch; b++)
    {
      printf("# Matrix %d \n", b);
      for (int j = 0; j < col; j++)
        for (int i = 0; i < row; i++)
        {
            int ix = i + row*j + b*row*col;
            printf("%1.16e\n", mat[ix]);
        }
    }
    // Release memory
    free(mat);

    printf("#  %d matrices written\n# %s done, exiting\n", batch, argv[0]);

    return 0;
}

