#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <cmath>
#include <iostream>
#include <ctime>


constexpr int N = 64;
constexpr double REYNOLDS_NUMBER = 1E6;
constexpr int GRAPHICS_UPDATE = 100;

constexpr int Q = 9;
constexpr double DENSITY = 2.7;
constexpr double LID_VELOCITY = 0.05;

constexpr int WIDTH  = 800;
constexpr int HEIGHT = 800;

float *scalar;

void showGraphics(double xmin, double xmax, double ymin, double ymax, const double *ux, const double *uy) {

  glViewport(0, 0, WIDTH, HEIGHT);
  glClearColor (1.0, 1.0, 1.0, 0.0);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glOrtho(xmin, xmax, ymin, ymax, -1.0, 1.0);

  glClear(GL_COLOR_BUFFER_BIT);

  const int NX = WIDTH;
  const int NY = HEIGHT;

  float dx = (xmax - xmin) / NX;
  float dy = (ymax - ymin) / NY;

  float min_curl = -0.02;
  float max_curl =  0.02;

  // Assign color to each pixel (i, j)
  for(int i = 0; i < NX-1; i++) {
	for(int j = 0; j < NY-1; j++) {

	  // map pixel coordinate (i,j) to LBM lattice coordinates (x,y)
	  int xin = i * N / NX;
	  int yin = j * N / NY;

	  // get locations of 4 data points inside which this pixel lies
	  int idx00 = (xin) * N + (yin);
	  int idx01 = (xin) * N + (yin + 1);
	  int idx10 = (xin+1) * N + (yin);
	  int idx11 = (xin+1) * N + (yin + 1);

	  //			   Neighbors
	  //
	  //               0p      1p 
	  //               |       |
	  //               |       |
	  //        m1-----01------11----p1
	  //               |       |
	  //               | Pixel |
	  //               |       |
	  //        m0-----00------10----p0
	  //               |       |
	  //               |       |
	  //               0m      1m

	  int idxm0 = (xin > 0)   ? (xin-1) * N + (yin  ) : idx00;
	  int idx0m = (yin > 0)   ? (xin  ) * N + (yin-1) : idx00;
	  int idx1m = (yin > 0)   ? (xin+1) * N + (yin-1) : idx10;
	  int idxp0 = (xin < N-1) ? (xin+2) * N + (yin  ) : idx10;
	  int idxp1 = (xin < N-1) ? (xin+2) * N + (yin+1) : idx11;
	  int idx1p = (yin < N-1) ? (xin+1) * N + (yin+2) : idx11;
	  int idx0p = (yin < N-1) ? (xin  ) * N + (yin+2) : idx01;
	  int idxm1 = (xin > 0)   ? (xin-1) * N + (yin+1) : idx01;

	  // calculate the normalized coordinates of the pixel
	  float xfl = (float)i * (float)N / (float) NX;
	  float yfl = (float)j * (float)N / (float) NY;
	  float x = xfl - (float)xin;
	  float y = yfl - (float)yin;

	  // calculate curl of the velocity field at the 4 data points
	  float dVdx_00 = uy[idx10] - uy[idxm0];
	  float dVdx_10 = uy[idxp0] - uy[idx00];
	  float dVdx_01 = uy[idx11] - uy[idxm1];
	  float dVdx_11 = uy[idxp1] - uy[idx01];

	  float dUdy_00 = ux[idx01] - ux[idx0m];
	  float dUdy_10 = ux[idx11] - ux[idx1m];
	  float dUdy_01 = ux[idx0p] - ux[idx00];
	  float dUdy_11 = ux[idx1p] - ux[idx10];

	  float curl_z_00 = dVdx_00 - dUdy_00;
	  float curl_z_10 = dVdx_10 - dUdy_10;
	  float curl_z_01 = dVdx_01 - dUdy_01;
	  float curl_z_11 = dVdx_11 - dUdy_11;

	  // bilinear interpolation
	  // float ux_interp = ux[idx00] * (1.0 - x) * (1.0 - y) + ux[idx10] * x * (1.0 - y) + ux[idx01] * (1.0 - x) * y + ux[idx11] * x * y;
	  // float uy_interp = uy[idx00] * (1.0 - x) * (1.0 - y) + uy[idx10] * x * (1.0 - y) + uy[idx01] * (1.0 - x) * y + uy[idx11] * x * y;
	  float curl_z_in = curl_z_00 * (1.0 - x) * (1.0 - y) + curl_z_10 * x * (1.0 - y) + curl_z_01 * (1.0 - x) * y + curl_z_11 * x * y;

	  // scalar[i*WIDTH + j] = pow((ux_interp*ux_interp + uy_interp*uy_interp), 0.5) / LID_VELOCITY;   // normalized velocity magnitude
	  scalar[i * WIDTH + j] = (max_curl - curl_z_in) / (max_curl - min_curl);                         // normalized vorticity

	  float x_actual = xmin + i * dx;
	  float y_actual = ymin + j * dy;
	  float VAL = scalar[i * WIDTH + j];

	  float R, G, B;

	  if(VAL <= 0.5) {
		// yellow to blue transition
		R = 2 * VAL;
		G = 2 * VAL;
		B = 1 - 2 * VAL;
	  } else {
		// red to yellow transition
		R = 1;
		G = 2 - 2 * VAL;
		B = 0;
	  }

	  glColor3f(std::abs(R), std::abs(G), std::abs(B));
	  glRectf(x_actual, y_actual, x_actual + dx, y_actual + dy);
	}
  }

}

void D3Q9(double *ex, double *ey, int *oppos, double *wt) {

  ex[0] =  0.0;   ey[0] =  0.0;   wt[0] = 4.0 /  9.0;
  ex[1] =  1.0;   ey[1] =  0.0;   wt[1] = 1.0 /  9.0;
  ex[2] =  0.0;   ey[2] =  1.0;   wt[2] = 1.0 /  9.0;
  ex[3] = -1.0;   ey[3] =  0.0;   wt[3] = 1.0 /  9.0;
  ex[4] =  0.0;   ey[4] = -1.0;   wt[4] = 1.0 /  9.0;
  ex[5] =  1.0;   ey[5] =  1.0;   wt[5] = 1.0 / 36.0;
  ex[6] = -1.0;   ey[6] =  1.0;   wt[6] = 1.0 / 36.0;
  ex[7] = -1.0;   ey[7] = -1.0;   wt[7] = 1.0 / 36.0;
  ex[8] =  1.0;   ey[8] = -1.0;   wt[8] = 1.0 / 36.0;

  oppos[0] = 0;      //      6        2        5
  oppos[1] = 3;      //               ^
  oppos[2] = 4;      //               |
  oppos[3] = 1;      //               |
  oppos[4] = 2;      //      3 <----- 0 -----> 1
  oppos[5] = 7;      //               |
  oppos[6] = 8;      //               |
  oppos[7] = 5;      //               v
  oppos[8] = 6;      //      7        4        8
}

void initialize(const int N, const int Q, const double DENSITY, const double LID_VELOCITY, 
	double *ex, double *ey, int *oppos, double *wt,
	double *rho, double *ux, double *uy, double* sigma, 
	double *f, double *feq, double *f_new) {

  for(int i = 0; i < N; i++) {
	for(int j = 0; j < N; j++) {

	  int index = i * N + j;

	  rho[index] = DENSITY;
	  ux[index] = 0.0;       // x-component of velocity
	  uy[index] = 0.0;       // x-component of velocity
	  sigma[index] = 0.0;    // rate-of-strain field

	  if(j == N - 1)
		ux[index] = LID_VELOCITY;

	  for(int a = 0;a < Q; a++) {

		int index_f = a + index * Q;

		double edotu = ex[a] * ux[index] + ey[a] * uy[index];
		double udotu = ux[index] * ux[index] + uy[index] * uy[index];

		feq[index_f]   = rho[index] * wt[a] * (1.0 + 3.0 * edotu + 4.5 * edotu * edotu - 1.5 * udotu);
		f[index_f]     = feq[index_f];
		f_new[index_f] = feq[index_f];

	  }

	}
  }
}

// this function updates the values of the distribution functions at all points along all directions
// carries out one lattice time-step (streaming + collision) in the algorithm
void collideAndStream(
	const int N, const int Q, const double DENSITY, const double LID_VELOCITY, const double REYNOLDS_NUMBER,
	const double *ex, const double *ey, const int *oppos, const double *wt,
	double *rho,         // density
	double *ux,         // X-velocity
	double *uy,         // Y-velocity
	double *sigma,      // rate-of-strain
	double *f,          // distribution function
	double *feq,        // equilibrium distribution function
	double *f_new)      // new distribution function
{

  for(int i = 1; i < N-1; i++) {
	for(int j = 1; j < N-1; j++) {

	  // natural index
	  int index = i * N + j;  // column-major ordering

	  // calculate fluid viscosity based on the Reynolds number
	  double kinematicViscosity = LID_VELOCITY * (double) N / REYNOLDS_NUMBER;

	  // calculate relaxation time tau
	  double tau =  0.5 + 3.0 * kinematicViscosity;

	  // collision
	  for(int a = 0; a < Q; a++) {
		int index_f = a + index*Q;
		double edotu = ex[a]*ux[index] + ey[a]*uy[index];
		double udotu = ux[index]*ux[index] + uy[index]*uy[index];
		feq[index_f] = rho[index] * wt[a] * (1 + 3*edotu + 4.5*edotu*edotu - 1.5*udotu);
	  }

	  // streaming from interior node points
	  for(int a = 0; a < Q; a++) {

		int index_f = a + index*Q;
		int index_nbr = (i+ex[a])*N + (j+ey[a]);
		int index_nbr_f = a + index_nbr * Q;
		int indexoppos = oppos[a] + index*Q;

		double tau_eff, tau_t, C_Smagorinsky;  // turbulence model parameters

		C_Smagorinsky = 0.16;

		tau_t = 0.5*(pow(pow(tau,2) + 18.0*pow(C_Smagorinsky,2)*sigma[index],0.5) - tau);
		tau_eff = tau + tau_t;

		// post-collision distribution at (i,j) along "a"
		double f_plus = f[index_f] - (f[index_f] - feq[index_f])/tau_eff;

		int iS = i + ex[a]; int jS = j + ey[a];

		if((iS == 0) || (iS == N-1) || (jS == 0) || (jS == N-1)) {
		  double ubdote = ux[index_nbr]*ex[a] + uy[index_nbr]*ey[a];
		  f_new[indexoppos] = f_plus - 6.0 * DENSITY * wt[a] * ubdote;
		} else {
		  f_new[index_nbr_f] = f_plus;
		}
	  }

	}
  }
}

void macroVar(
	const int N, const int Q, const double DENSITY, const double LID_VELOCITY, const double REYNOLDS_NUMBER,
	const double *ex, const double *ey, const int *oppos, const double *wt,
	double *rho,         // density
	double *ux,         // X-velocity
	double *uy,         // Y-velocity
	double *sigma,      // rate-of-strain
	double *f,          // distribution function
	double *feq,        // equilibrium distribution function
	double *f_new)      // new distribution function
{

  for(int i = 1; i < N-1; i++) {
	for(int j = 1; j < N-1; j++) {

	  // natural index
	  int index = i * N + j;  // column-major ordering

	  // push f_new into f
	  for(int a = 0; a < Q; a++) {
		int index_f = a + index*Q;
		f[index_f] = f_new[index_f];
	  }

	  // update density at interior nodes
	  rho[index] = 0.0;
	  for(int a = 0; a < Q; a++) {
		int index_f = a + index*Q;
		rho[index] += f_new[index_f];
	  }

	  // update velocity at interior nodes
	  double velx = 0.0;
	  double vely = 0.0;
	  for(int a = 0; a < Q; a++) {
		int index_f = a + index*Q;
		velx += f_new[index_f]*ex[a];
		vely += f_new[index_f]*ey[a];
	  }
	  ux[index] = velx/rho[index];
	  uy[index] = vely/rho[index];

	  // update the rate-of-strain field
	  double sum_xx = 0.0, sum_xy = 0.0, sum_xz = 0.0;
	  double sum_yx = 0.0, sum_yy = 0.0, sum_yz = 0.0;
	  double sum_zx = 0.0, sum_zy = 0.0, sum_zz = 0.0;
	  for(int a = 1; a < Q; a++) {
		int index_f = a + index*Q;

		sum_xx = sum_xx + (f_new[index_f] - feq[index_f])*ex[a]*ex[a];
		sum_xy = sum_xy + (f_new[index_f] - feq[index_f])*ex[a]*ey[a];
		sum_xz = 0.0;
		sum_yx = sum_xy;
		sum_yy = sum_yy + (f_new[index_f] - feq[index_f])*ey[a]*ey[a];
		sum_yz = 0.0;
		sum_zx = 0.0;
		sum_zy = 0.0;
		sum_zz = 0.0;
	  }

	  // evaluate |S| (magnitude of the strain-rate)
	  sigma[index] = pow(sum_xx,2) + pow(sum_xy,2) + pow(sum_xz,2) + pow(sum_yx,2) + pow(sum_yy,2) + pow(sum_yz,2) + pow(sum_zx,2) + pow(sum_zy,2) + pow(sum_zz,2);
	  sigma[index] = pow(sigma[index],0.5);

	}
  }
}

int main(int argc, char* argv[]) {

  GLFWwindow *window;

  if(!glfwInit())
	return -1;

  window = glfwCreateWindow(WIDTH, HEIGHT, "Lattice Boltzman", NULL, NULL);
  if(!window) {
	glfwTerminate();
	return -1;
  }

  glewInit();

  glfwMakeContextCurrent(window);

  // distribution functions
  double *f = new double[N*N*Q];
  double *feq = new double[N*N*Q];
  double *f_new = new double[N*N*Q];

  // density and velocity
  double *rho = new double[N*N];
  double *ux = new double[N*N];
  double *uy = new double[N*N];

  // rate-of-strain
  double *sigma = new double[N*N];

  // D3Q9 parameters
  double *ex = new double[Q];
  double *ey = new double[Q];
  int *oppos = new int[Q];
  double *wt = new double[Q];

  scalar = new float[WIDTH * HEIGHT];

  D3Q9(ex, ey, oppos, wt);
  initialize(N, Q, DENSITY, LID_VELOCITY, ex, ey, oppos, wt, rho, ux, uy, sigma, f, feq, f_new);

  int time = 0;
  double xmin = 0, xmax = N, ymin = 0, ymax = N;

  while(!glfwWindowShouldClose(window)) {
	time++;     

	collideAndStream(N, Q, DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, ex, ey, oppos, wt, rho, ux, uy, sigma, f, feq, f_new);
	macroVar(N, Q, DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, ex, ey, oppos, wt, rho, ux, uy, sigma, f, feq, f_new);

	if(time % GRAPHICS_UPDATE == 0) {
	  showGraphics(xmin, xmax, ymin, ymax, ux, uy);
	  glfwSwapBuffers(window);
	  glfwPollEvents();
	}
  }

  delete [] f;
  delete [] feq;
  delete [] f_new;
  delete [] rho;
  delete [] ux;
  delete [] uy;
  delete [] sigma;
  delete [] ex;
  delete [] ey;
  delete [] oppos;
  delete [] wt;

  delete[] scalar;

  glfwDestroyWindow(window);
  glfwTerminate();

  return 0;
}
