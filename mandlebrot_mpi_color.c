#include<stdio.h>
#include<assert.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>

typedef struct complex{
  double real;
  double imag;
} Complex;

typedef unsigned char Color[3];

const double MAX_LENGTH_SQ = 4;
const int MAX_ITER = 10000;
const double MAX_DISTANCE = 1 << 18;


#pragma acc routine(hsv_to_rgb) seq
void hsv_to_rgb(const double *hsv, double *rgb) {
  rgb[0] = 0; rgb[1] = 0; rgb[2] = 0;
  double c = hsv[1]*hsv[2];
  double h_prime = hsv[0]/60;
  double x = c*(1 - abs(fmod(h_prime,2.0) - 1));
  if(h_prime >= 0 && h_prime <= 1) {
    rgb[0] = c; rgb[1] = x;
  }
  else if(h_prime >= 1 && h_prime <= 2) {
    rgb[0] = x; rgb[1] = c;
  }
  else if(h_prime >= 2 && h_prime <= 3) {
    rgb[1] = c; rgb[2] = x;
  }
  else if(h_prime >= 3 && h_prime <= 4) {
    rgb[1] = x; rgb[2] = c;
  }
  else if(h_prime >= 4 && h_prime <= 5) {
    rgb[2] = c; rgb[0] = x;
  }
  else if(h_prime >= 5 && h_prime <= 6) {
    rgb[2] = x; rgb[0] = c;
  }
  double m = hsv[2] - c;
  for(int i = 0; i < 3; i++) {
    rgb[i] += m;
  }
} 

#pragma acc routine(cal_color) seq
Color cal_color(double pixel_size, double distance, int iter) {
  Color c;
  if(iter >= MAX_ITER) {
    c[0] = 255; c[1] = 255; c[2] = 255;
    return c;
  }
  // compute hsv
  double temp;
  temp = log(iter) / log(MAX_ITER);
  temp *= 10; temp -= floor(temp);
  double hsv[3], rgb[3];
  hsv[0] = temp;
  hsv[1] = 0.7;
  if(distance < 0.5*pixel_size) {
    temp = pow(2*distance/pixel_size,1./3);
    hsv[2] = temp;
  }
  else {
    hsv[2] = 1;
  }
  // convert to hsv
  hsv_to_rgb(hsv,rgb);
  // convert to char
  for(int i = 0; i < 3; ++i) {
    c[i] = (char)255*rgb[i];
  }
  return c;
}

#pragma acc routine(cal_pixel) seq
Color cal_pixel(double pixel_size, Complex pt){
  int iter;
  Complex z, dz;
  double temp, length_sq, distance;
  z.real = 0; z.imag = 0;
  dz.real = 0; dz.real = 0;
  iter = 0;

  do{
    temp = z.real * dz.real - z.imag * dz.imag;
    dz.imag = 4 * dz.real * z.imag;
    dz.real = 2 * temp + 1; 
    temp = z.real * z.real - z.imag * z.imag + pt.real;
    z.imag = 2 * z.real * z.imag + pt.imag;
    z.real = temp;
    length_sq = z.real * z.real + z.imag * z.imag;
    iter++;
  }
  while ((length_sq < MAX_LENGTH_SQ) && (iter < MAX_ITER));
  temp = dz.real * dz.real + dz.imag + dz.imag;
  distance = 2*log(length_sq) * length_sq / temp;
  return cal_color(pixel_size,distance,iter);
}

#define MASTERPE 0
int main(int argc, char **argv){
  FILE *file;
  int i, j;
  Complex pt;
  int tmp;
  // for writing out integers
  Color *data;
  int nx, ny;
  int mystrt, myend;
  int nrows_l;
  int nprocs, mype;
  double pixel_size;

  MPI_Status status;

  /* regular MPI initialization stuff */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  if (argc != 3){
    int err = 0;
    printf("argc %d\n", argc);
    if (mype == MASTERPE){
      printf("usage: mandelbrot nx ny");
      MPI_Abort(MPI_COMM_WORLD,err );
    }
  }

  /* get command line args */
  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  pixel_size = 4.0/((double)nx);

  /* assume divides equally */
  nrows_l = nx/nprocs;

  /* create buffer for local work only */
  data = (Color*)malloc(nrows_l * ny * sizeof(Color));
#pragma acc enter data create(data[0:nrows_l*ny])
  //  data_l_tmp = data_l;

  /* calculate each processor's region of work */
  mystrt = mype*nrows_l;
  myend  = mystrt + nrows_l - 1;

  /* calc each procs coordinates and call local mandelbrot set function */
#pragma acc parallel loop private(c, tmp) collapse(2) present(data[0:nrows_l*ny])
  for (i = mystrt; i <= myend; ++i){
    for (j = 0; j < ny; ++j){
      pt.real = i/((double) nx) * 4. - 2.;
      pt.imag = j/((double) ny) * 4. - 2.;
      data[i*nx + j] = cal_pixel(pixel_size,pt);
    }
  }
#pragma acc exit data copyout(data[0:nrows_l*ny])
  if (mype == MASTERPE){
    file = fopen("mandelbrot.bin", "w");
    printf("nrows_l, ny  %d %d\n", nrows_l, ny);
    fprintf(file,"P5\n");
    fprintf(file,"%d %d\n",nrows_l,ny);
    fprintf(file,"%d\n",255);
    fwrite(data, nrows_l*ny, sizeof(Color), file);
    fclose(file);

    for (i = 1; i < nprocs; ++i){
      //MPI_Recv(data_l, nrows_l * ny, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      //printf("received message from proc %d\n", i);
      //file = fopen("mandelbrot.bin_0000", "a");
      //fwrite(data_l, nrows_l*ny, sizeof(double), file);
      //fclose(file);
    }
  }

  else{
    //MPI_Send(data_l, nrows_l * ny, MPI_DOUBLE, MASTERPE, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
}
