#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define FFT_SIZE 16384
#define DFT_SIZE  8        // When we DFT rather than FFT
#define I_SCALE      2     // Scale factor for the integer FFT data (e.g 24 bits -> 24+I_SCALE bits)
#define I_TRIG_SCALE 22    // Scale factor for the integer TRIG constants (e.g cos(0) = 1<<I_TRIG_SCALE)

static double    d_sines[FFT_SIZE], d_cosines[FFT_SIZE];
static float     f_sines[FFT_SIZE], f_cosines[FFT_SIZE];
static long long i_sines[FFT_SIZE], i_cosines[FFT_SIZE];

static double    d_din_r[FFT_SIZE], d_din_i[FFT_SIZE], d_dout_r[FFT_SIZE], d_dout_i[FFT_SIZE];
static float     f_din_r[FFT_SIZE], f_din_i[FFT_SIZE], f_dout_r[FFT_SIZE], f_dout_i[FFT_SIZE];
long long        i_din_r[FFT_SIZE], i_din_i[FFT_SIZE], i_dout_r[FFT_SIZE], i_dout_i[FFT_SIZE];

double           d_dout_r_ref[FFT_SIZE], d_dout_i_ref[FFT_SIZE];

//=================================================================
// Save calling sin() and cos() all the time
//=================================================================
static void init_tables(void) {
   for(int i = 0; i < FFT_SIZE; i++) {
     double phase = 2.0*M_PI*i/FFT_SIZE;
     d_cosines[i] = cos(phase);
     d_sines[i]   = sin(phase);
     f_cosines[i] = cos(phase);
     f_sines[i]   = sin(phase);
     i_cosines[i] = cos(phase)*(1<<I_TRIG_SCALE);
     i_sines[i]   = sin(phase)*(1<<I_TRIG_SCALE);
  }
}

static void print_out( double    *d_r, double    *d_i, 
                       float     *f_r, float     *f_i, 
                       long long *i_r, long long *i_i, 
                       int size) {
   printf("bin  ====  Doubles ========== "
              " ==== Floats ===========  == err (double-float) == "
              " ====== Fixed point ====  == err (double-fixed) ==\n");
   for(int i = 0; i < size; i++) {
      if(i < 128 || i > FFT_SIZE-129) {
         printf("%3i, (%10.2f, %10.2f), "
                  "(%10.2f,%10.2f), (%10.2f,%10.2f), "
                  "(%10.2f,%10.2f), (%10.2f, %10.2f)\n", 
              i,
              d_r[i], d_i[i],
              f_r[i], f_i[i],
              fabs(d_r[i]-f_r[i]),                  fabs(d_i[i]-f_i[i]),
              i_r[i]*1.0/(1<<I_SCALE), i_i[i]*1.0/(1<<I_SCALE),
              fabs(d_r[i]-i_r[i]*1.0/(1<<I_SCALE)), fabs(d_i[i]-i_i[i]*1.0/(1<<I_SCALE))
        );
     }
   }
}


static int reverse_bits(int i, int max) {
   int o = 0;
   i |= max;
   while(i != 1) {
      o <<= 1;
      if(i&1) 
        o |= 1; 
      i >>= 1;
   }
   return o;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
static void d_dft(double *r_i, double *i_i, double *r_o, double *i_o, int size, int stride) {
   for(int bin = 0; bin < size; bin++) {
      int i = 0;
      double total_i = 0.0;
      double total_q = 0.0;
      for(int s = 0; s < size; s++) {
         total_i +=  r_i[s*stride] * d_cosines[i] - i_i[s*stride] *   d_sines[i];
         total_q +=  r_i[s*stride] *   d_sines[i] + i_i[s*stride] * d_cosines[i];
         i += bin*(FFT_SIZE/size);
         if(i >= FFT_SIZE) i -= FFT_SIZE;
      }
      r_o[bin] = total_i/FFT_SIZE;
      i_o[bin] = total_q/FFT_SIZE;
   }
}

static void d_fft(double *r_i, double *i_i, double *r_o, double *i_o, int size) {
    int i, stride, step;

    stride = size/DFT_SIZE;
    for(i = 0; i < stride ; i++) {
        int out_offset = reverse_bits(i,stride) * DFT_SIZE;
        d_dft(r_i+i, i_i+i,  r_o+out_offset, i_o+out_offset, DFT_SIZE, stride);
    }

    stride = DFT_SIZE*2;
    step = FFT_SIZE/stride;
    while(stride <= FFT_SIZE) {
       for(i = 0; i < FFT_SIZE; i+= stride) {
          double *real = r_o+i;
          double *imag = i_o+i;
          size = stride/2;
          for(int j = 0; j < size; j++) {
             double c = d_cosines[j*step], s = d_sines[j*step];
             double rotated_r =  real[size] * c - imag[size] * s;
             double rotated_i =  real[size] * s + imag[size] * c;
             real[size] = real[0] - rotated_r;
             imag[size] = imag[0] - rotated_i;
             real[0]    = real[0] + rotated_r;
             imag[0]    = imag[0] + rotated_i;
             real++;
             imag++;
          }
       }
       stride *= 2;
       step   /= 2;
    }
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
static void f_dft(float *r_i, float *i_i, float *r_o, float *i_o, int size, int stride) {
   for(int bin = 0; bin < size; bin++) {
      int i = 0;
      float total_i = 0.0;
      float total_q = 0.0;
      for(int s = 0; s < size; s++) {
         total_i +=  r_i[s*stride] * d_cosines[i] - i_i[s*stride] *   d_sines[i];
         total_q +=  r_i[s*stride] *   d_sines[i] + i_i[s*stride] * d_cosines[i];
         i += bin*(FFT_SIZE/size);
         if(i >= FFT_SIZE) i -= FFT_SIZE;
      }
      r_o[bin] = total_i/FFT_SIZE;
      i_o[bin] = total_q/FFT_SIZE;
   }
}

static void f_fft(float *r_i, float *i_i, float *r_o, float *i_o, int size) {
    int i, stride, step;

    stride = size/DFT_SIZE;
    for(i = 0; i < stride ; i++) {
        int out_offset = reverse_bits(i,stride) * DFT_SIZE;
        f_dft(r_i+i, i_i+i,  r_o+out_offset, i_o+out_offset, DFT_SIZE, stride);
    }

    stride = DFT_SIZE*2;
    step = FFT_SIZE/stride;
    while(stride <= FFT_SIZE) {
       for(i = 0; i < FFT_SIZE; i+= stride) {
          float *real = r_o+i;
          float *imag = i_o+i;
          size = stride/2;
          for(int j = 0; j < size; j++) {
             float c = d_cosines[j*step], s = d_sines[j*step];
             float rotated_r =  real[size] * c - imag[size] * s;
             float rotated_i =  real[size] * s + imag[size] * c;
             real[size] = real[0] - rotated_r;
             imag[size] = imag[0] - rotated_i;
             real[0]    = real[0] + rotated_r;
             imag[0]    = imag[0] + rotated_i;
             real++;
             imag++;
          }
       }
       stride *= 2;
       step   /= 2;
    }
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
static void i_dft(long long *r_i, long long *i_i, long long *r_o, long long *i_o, int size, int stride) {
   for(int bin = 0; bin < size; bin++) {
      int i = 0;
      long long total_i = 0.0;
      long long total_q = 0.0;
      for(int s = 0; s < size; s++) {
         total_i +=  r_i[s*stride] * i_cosines[i] - i_i[s*stride] *   i_sines[i];
         total_q +=  r_i[s*stride] *   i_sines[i] + i_i[s*stride] * i_cosines[i];
         i += bin*(FFT_SIZE/size);
         if(i >= FFT_SIZE) i -= FFT_SIZE;
      }
      r_o[bin] = total_i >> I_TRIG_SCALE;
      i_o[bin] = total_q >> I_TRIG_SCALE;
   }
}

static void i_fft(long long *r_i, long long *i_i, long long *r_o, long long *i_o, int size) {
    int i, stride, step;

    stride = size/DFT_SIZE;
    for(i = 0; i < stride ; i++) {
        int out_offset = reverse_bits(i,stride) * DFT_SIZE;
        i_dft(r_i+i, i_i+i,  r_o+out_offset, i_o+out_offset, DFT_SIZE, stride);
    }

    stride = DFT_SIZE*2;
    step = FFT_SIZE/stride;
    while(stride <= FFT_SIZE) {
       for(i = 0; i < FFT_SIZE; i+= stride) {
          long long *real = r_o+i;
          long long *imag = i_o+i;
          size = stride/2;
          for(int j = 0; j < size; j++) {
             long long c = i_cosines[j*step], s = i_sines[j*step];
             long long rotated_r =  (real[size] * c - imag[size] * s)>>I_TRIG_SCALE;
             long long rotated_i =  (real[size] * s + imag[size] * c)>>I_TRIG_SCALE;
             real[size] = real[0] - rotated_r;
             imag[size] = imag[0] - rotated_i;
             real[0]    = real[0] + rotated_r;
             imag[0]    = imag[0] + rotated_i;
             real++;
             imag++;
          }
       }
       stride *= 2;
       step   /= 2;
    }
    for(i = 0; i < FFT_SIZE ; i++) {
       r_o[i] /= FFT_SIZE;
       i_o[i] /= FFT_SIZE;
    }
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

void ts_sub(struct timespec *r, struct timespec *a, struct timespec *b) {
   r->tv_sec = a->tv_sec - b->tv_sec;
   if(a->tv_nsec > b->tv_nsec) {
      r->tv_nsec = a->tv_nsec - b->tv_nsec;
   } else {
      r->tv_sec--;
      r->tv_nsec = a->tv_nsec - b->tv_nsec+1000000000;
   }
}

int main(int argc, char *argv[]) {
   // Setup
   int range = (1<<24);
   int i_min, i_max;
   for(int i = 0; i < FFT_SIZE; i++) {
      d_din_r[i] =(rand()%(range)-range/2);
      d_din_i[i] =(rand()%(range)-range/2);
      f_din_r[i] = d_din_r[i];
      f_din_i[i] = d_din_i[i];
      i_din_r[i] = d_din_r[i]*(1<<I_SCALE);
      i_din_i[i] = d_din_i[i]*(1<<I_SCALE);

      if(i == 0)             i_min = i_din_r[i];
      if(i_din_r[i] < i_min) i_min = i_din_r[i];
      if(i_din_i[i] < i_min) i_min = i_din_r[i];

      if(i == 0)             i_max = i_din_r[i];
      if(i_din_r[i] > i_max) i_max = i_din_r[i];
      if(i_din_i[i] > i_max) i_max = i_din_r[i];
   }
   init_tables();
   printf("Input integer data ranges from %i to %i\n", i_min>>I_SCALE, i_max>>I_SCALE);
   printf("Transform of %5i random complex numbers\n", FFT_SIZE);
   printf("=========================================\n");

   struct timespec tv_d_dft_start, tv_d_dft_end, tv_d_dft;
   struct timespec tv_d_fft_start, tv_d_fft_end, tv_d_fft;
   struct timespec tv_f_fft_start, tv_f_fft_end, tv_f_fft;
   struct timespec tv_i_fft_start, tv_i_fft_end, tv_i_fft;

   clock_gettime(CLOCK_MONOTONIC, &tv_d_dft_start);  
   d_dft(d_din_r,  d_din_i, d_dout_r_ref, d_dout_i_ref, FFT_SIZE,1);
   clock_gettime(CLOCK_MONOTONIC, &tv_d_dft_end);  

   clock_gettime(CLOCK_MONOTONIC, &tv_d_fft_start);  
   d_fft(d_din_r,  d_din_i, d_dout_r,     d_dout_i,     FFT_SIZE);
   clock_gettime(CLOCK_MONOTONIC, &tv_d_fft_end);  

   clock_gettime(CLOCK_MONOTONIC, &tv_f_fft_start);  
   f_fft(f_din_r,  f_din_i, f_dout_r,     f_dout_i,     FFT_SIZE);
   clock_gettime(CLOCK_MONOTONIC, &tv_f_fft_end);  

   clock_gettime(CLOCK_MONOTONIC, &tv_i_fft_start);  
   i_fft(i_din_r,  i_din_i, i_dout_r,     i_dout_i,     FFT_SIZE);
   clock_gettime(CLOCK_MONOTONIC, &tv_i_fft_end);  


   print_out( d_dout_r, d_dout_i, f_dout_r, f_dout_i, i_dout_r, i_dout_i, FFT_SIZE);

   printf("Time taken\n==========\n");
   ts_sub(&tv_d_dft, &tv_d_dft_end, &tv_d_dft_start);
   printf("Double DFT            %u.%09u secs\n",(unsigned)tv_d_dft.tv_sec, (unsigned)tv_d_dft.tv_nsec);

   ts_sub(&tv_d_fft, &tv_d_fft_end, &tv_d_fft_start);
   printf("Double FFT            %u.%09u secs\n",(unsigned)tv_d_fft.tv_sec, (unsigned)tv_d_fft.tv_nsec);

   ts_sub(&tv_f_fft, &tv_f_fft_end, &tv_f_fft_start);
   printf("Float FFT             %u.%09u secs\n",(unsigned)tv_f_fft.tv_sec, (unsigned)tv_f_fft.tv_nsec);

   ts_sub(&tv_i_fft, &tv_i_fft_end, &tv_i_fft_start);
   printf("Integer FFT           %u.%09u secs\n",(unsigned)tv_i_fft.tv_sec, (unsigned)tv_i_fft.tv_nsec);
}

