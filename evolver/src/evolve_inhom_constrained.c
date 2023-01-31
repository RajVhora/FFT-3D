#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include <fftw3.h>

#include "../headers/functions.h"
#include "../headers/nrutil.h"

void evolve_inhom_constrained(int n_x, int n_y, int n_z, double delta_x,
							  double delta_y, double delta_z, double kappa, double A, double delta_t, int time_steps,
							  double t0, fftw_complex *comp, double Ceff[4][4][4][4],
							  double DeltaC[4][4][4][4], double **sigma_T, double **epsilon_T,
							  double ave_comp, double **sig_app, int STEPS, int MAXITR, double MAXERR,
							  int n_alpha, int n_beta, double *alpha, double *beta, double *alpha_prime,
							  double *beta_prime, int plan_indicator, int SNA, int NOISEADDER,
							  double noise_str)
{

	FILE *fpd;
	FILE *fptmp;
	FILE *fpstr;

	int INDEX = 0;
	int i1, i2, i3, i4, i5, i6, i7;
	int J;

	int half_nx, half_ny, half_nz;
	double *kx, *ky, *kz;
	double delta_kx, delta_ky, delta_kz;
	double k2, k4;
	double inv_denom;

	int nxnynz;
	double inv_nxnynz;

	char NAME[50];

	double epsstr11, epsstr22, epsstr33, epsstr12, epsstr13, epsstr23;

	double *c;
	double **Del_sigma_T;
	double **E;
	double **strain;
	double *Omega11;
	double *Omega22;
	double *Omega33;
	double *Omega12;
	double *Omega21;
	double *Omega13;
	double *Omega31;
	double *Omega23;
	double *Omega32;
	double S[4][4][4][4];

	double temp1, temp2, temp3, temp4, temp5, temp6;

	double a, b;
	double ap, bp;
	double realc;
	size_t tmp;

	unsigned FLAG;

	fftw_complex *g;
	fftw_complex *u1_old, *u2_old, *u3_old;
	fftw_complex *u1_new, *u2_new, *u3_new;
	fftw_complex *eps_star11, *eps_star22, *eps_star33, *eps_star12, *eps_star13, *eps_star23;
	fftw_complex *eta;
	fftw_complex *mu_el;

	fftw_plan planF, planB;

	nxnynz = n_x * n_y * n_z;
	inv_nxnynz = 1.0 / nxnynz;

	g = fftw_malloc(nxnynz * sizeof(fftw_complex));
	u1_old = fftw_malloc(nxnynz * sizeof(fftw_complex));
	u2_old = fftw_malloc(nxnynz * sizeof(fftw_complex));
	u3_old = fftw_malloc(nxnynz * sizeof(fftw_complex));
	u1_new = fftw_malloc(nxnynz * sizeof(fftw_complex));
	u2_new = fftw_malloc(nxnynz * sizeof(fftw_complex));
	u3_new = fftw_malloc(nxnynz * sizeof(fftw_complex));
	eps_star11 = fftw_malloc(nxnynz * sizeof(fftw_complex));
	eps_star22 = fftw_malloc(nxnynz * sizeof(fftw_complex));
	eps_star33 = fftw_malloc(nxnynz * sizeof(fftw_complex));
	eps_star12 = fftw_malloc(nxnynz * sizeof(fftw_complex));
	eps_star13 = fftw_malloc(nxnynz * sizeof(fftw_complex));
	eps_star23 = fftw_malloc(nxnynz * sizeof(fftw_complex));
	eta = fftw_malloc(2 * 2 * sizeof(fftw_complex));
	mu_el = fftw_malloc(nxnynz * sizeof(fftw_complex));

	if (plan_indicator == 0)
	{
		FLAG = FFTW_ESTIMATE;
	}
	else if (plan_indicator == 1)
	{
		FLAG = FFTW_MEASURE;
	}
	else if (plan_indicator == 2)
	{
		FLAG = FFTW_PATIENT;
	}
	else if (plan_indicator == 3)
	{
		FLAG = FFTW_EXHAUSTIVE;
	}

	// fftw_plan_with_nthreads(8);

	planF =
		fftw_plan_dft_3d(n_x, n_y, n_z, mu_el, mu_el, FFTW_FORWARD, FLAG);
	planB =
		fftw_plan_dft_3d(n_x, n_y, n_z, eps_star11, eps_star11, FFTW_BACKWARD, FLAG);

	/* Declare and initialise the acoustic tensor */

	half_nx = (int)n_x / 2;
	half_ny = (int)n_y / 2;
	half_nz = (int)n_z / 2;

	delta_kx = (2.0 * M_PI) / (n_x * delta_x);
	delta_ky = (2.0 * M_PI) / (n_y * delta_y);
	delta_kz = (2.0 * M_PI) / (n_z * delta_z);

	kx = (double *)malloc((size_t)n_x * sizeof(double));
	ky = (double *)malloc((size_t)n_y * sizeof(double));
	kz = (double *)malloc((size_t)n_z * sizeof(double));

	Omega11 = dvector(0, nxnynz);
	Omega22 = dvector(0, nxnynz);
	Omega33 = dvector(0, nxnynz);
	Omega12 = dvector(0, nxnynz);
	Omega21 = dvector(0, nxnynz);
	Omega13 = dvector(0, nxnynz);
	Omega31 = dvector(0, nxnynz);
	Omega23 = dvector(0, nxnynz);
	Omega32 = dvector(0, nxnynz);

	for (i1 = 0; i1 < n_x; ++i1)
	{
		if (i1 < half_nx)
			kx[i1] = i1 * delta_kx;
		else
			kx[i1] = (i1 - n_x) * delta_kx;
	}
	for (i2 = 0; i2 < n_y; ++i2)
	{
		if (i2 < half_ny)
			ky[i2] = i2 * delta_ky;
		else
			ky[i2] = (i2 - n_y) * delta_ky;
	}
	for (i3 = 0; i3 < n_z; ++i3)
	{
		if (i3 < half_nz)
			kz[i3] = i3 * delta_kz;
		else
			kz[i3] = (i3 - n_z) * delta_kz;
	}

	for (i1 = 0; i1 < n_x; ++i1)
	{
		for (i2 = 0; i2 < n_y; ++i2)
		{
			for (i3 = 0; i3 < n_z; ++i3)
			{
				J = i3 + n_z * (i2 + n_y * i1);
				Omega11[J] = 0.0;
				Omega22[J] = 0.0;
				Omega33[J] = 0.0;
				Omega12[J] = 0.0;
				Omega21[J] = 0.0;
				Omega13[J] = 0.0;
				Omega31[J] = 0.0;
				Omega23[J] = 0.0;
				Omega32[J] = 0.0;
			}
		}
	}

	/* Write the initial configuration */

	c = (double *)malloc((size_t)nxnynz * sizeof(double));
	for (i1 = 0; i1 < n_x; ++i1)
	{
		for (i2 = 0; i2 < n_y; ++i2)
		{
			for (i3 = 0; i3 < n_z; ++i3)
			{
				J = i3 + n_z * (i2 + n_y * i1);
				c[J] = creal(comp[J]);
			}
		}
	}

	sprintf(NAME, "../output/data/time%d.dat",
			(int)(INDEX + t0));
	fpd = fopen(NAME, "w");
	tmp = fwrite(&c[0], sizeof(double), (size_t)nxnynz, fpd);
	fclose(fpd);

	/** Calculate the Del_sigma_T tensor **/

	Del_sigma_T = dmatrix(1, 3, 1, 3);
	calculate_Del_sigma_T(DeltaC, epsilon_T, Del_sigma_T);
	E = dmatrix(1, 3, 1, 3);
	strain = dmatrix(1, 3, 1, 3);

	/* Calculate the homogeneous strain (if there is applied stress). In
	 * the absence of applied strain, we assume the homogeneous strain to
	 * be zero */

	calculate_S_exp(n_x, n_y, n_z, delta_x, delta_y, delta_z, comp, ave_comp, Ceff, DeltaC, S, n_alpha,
					alpha);
	/*
	for (i1 = 1; i1 < 4; ++i1)
	{
		for (i2 = 1; i2 < 4; ++i2)
		{
			E[i1][i2] = 0.0;
		}
	}

	for (i1 = 1; i1 < 4; ++i1)
	{
		for (i2 = 1; i2 < 4; ++i2)
		{
			for (i3 = 1; i3 < 4; ++i3)
			{
				for (i4 = 1; i4 < 4; ++i4)
				{
					E[i1][i2] = E[i1][i2] + S[i1][i2][i3][i4] * sig_app[i3][i4];
				}
			}
		}
	}
	*/
	/*
	E[1][1] = E[2][2] = 0.01;
	E[1][2] = E[2][1] = 0.0;
	*/

	printf("%le\n", E[1][1]);
	printf("%le\n", E[2][2]);
	printf("%le\n", E[3][3]);
	printf("%le\n", E[1][2]);
	printf("%le\n", E[2][1]);
	printf("%le\n", E[1][3]);
	printf("%le\n", E[3][1]);
	printf("%le\n", E[2][3]);
	printf("%le\n", E[3][2]);

	/* Evolve */

	fptmp = fopen("../output/noise_info", "a");

	for (INDEX = 0; INDEX < time_steps + 1; ++INDEX)
	{
		printf("INDEX %d\n", INDEX);
		/* Calculate g and its Fourier transform */

		for (i1 = 0; i1 < n_x; ++i1)
		{
			for (i2 = 0; i2 < n_y; ++i2)
			{
				for (i3 = 0; i3 < n_z; ++i3)
				{
					J = i3 + n_z * (i2 + n_y * i1);
					realc = creal(comp[J]);
					g[J] = 2.0 * A * realc * (1.0 - realc) * (1.0 - 2.0 * realc) + 0.0 * _Complex_I;
				}
			}
		}

		fftw_execute_dft(planF, g, g);

		if (INDEX == 0)
		{

			/* The acoustic tensor Omega and the (homogeneous) base solution */

			calculate_Omega(n_x, n_y, n_z, half_nx, half_ny, half_nz, kx, ky, kz, Ceff,
							Omega11, Omega12, Omega21, Omega22, Omega13, Omega31, Omega33, Omega23, Omega32);

			calculate_uzero(n_x, n_y, n_z, half_nx, half_ny, half_nz, delta_kx, delta_ky, delta_kz, ave_comp,
							comp, Ceff, sigma_T, u1_old, u2_old, u3_old, n_beta, beta, Omega11, Omega12, Omega21,
							Omega22, Omega13, Omega31, Omega33, Omega23, Omega32, planF);
		}

		/* Refine the displacements */

		refine_u(n_x, n_y, n_z, half_nx, half_ny, half_nz, delta_x, delta_y, delta_z, epsilon_T,
				 delta_kx, delta_ky, delta_kz, comp, ave_comp, Ceff, DeltaC, sigma_T,
				 u1_old, u2_old, u3_old, u1_new, u2_new, u3_new, Del_sigma_T, sig_app, E, MAXITR, MAXERR,
				 n_alpha, n_beta, alpha, beta, Omega11, Omega12, Omega21, Omega22, Omega13, Omega31, Omega33, Omega23, Omega32, planF, planB);

		for (i1 = 0; i1 < n_x; ++i1)
		{
			for (i2 = 0; i2 < n_y; ++i2)
			{
				for (i3 = 0; i3 < n_z; ++i3)
				{
					J = i3 + n_z * (i2 + n_y * i1);
					eps_star11[J] = _Complex_I * kx[i1] * u1_new[J];
					eps_star22[J] = _Complex_I * ky[i2] * u2_new[J];
					eps_star33[J] = _Complex_I * kz[i3] * u3_new[J];
					eps_star12[J] = 0.5 * _Complex_I * (kx[i1] * u2_new[J] + ky[i2] * u1_new[J]);
					eps_star13[J] = 0.5 * _Complex_I * (kx[i1] * u3_new[J] + kz[i3] * u1_new[J]);
					eps_star23[J] = 0.5 * _Complex_I * (ky[i2] * u3_new[J] + kz[i3] * u2_new[J]);
				}
			}
		}

		/* Get the strain back to the real space */

		fftw_execute(planB);
		fftw_execute_dft(planB, eps_star22, eps_star22);
		fftw_execute_dft(planB, eps_star33, eps_star33);
		fftw_execute_dft(planB, eps_star12, eps_star12);
		fftw_execute_dft(planB, eps_star13, eps_star13);
		fftw_execute_dft(planB, eps_star23, eps_star23);

		/* Calculate mu_el and take it to the Fourier space */
		 //if(INDEX != 0 && INDEX%STEPS ==0){
		// Declare file pointer
		FILE *fpstr;
		if(INDEX%STEPS ==0){
		 sprintf(NAME,"../output/data/stress_at_time%d.dat",(int) (INDEX + t0));
		 fpstr = fopen(NAME,"w");
		 }

		for (i1 = 0; i1 < n_x; ++i1)
		{
			for (i2 = 0; i2 < n_y; ++i2)
			{
				for (i7 = 0; i7 < n_z; ++i7)
				{

					J = i7 + n_z * (i2 + n_y * i1);

					epsstr11 = creal(eps_star11[J] * inv_nxnynz);
					epsstr22 = creal(eps_star22[J] * inv_nxnynz);
					epsstr33 = creal(eps_star33[J] * inv_nxnynz);
					epsstr12 = creal(eps_star12[J] * inv_nxnynz);
					epsstr13 = creal(eps_star13[J] * inv_nxnynz);
					epsstr23 = creal(eps_star23[J] * inv_nxnynz);

					realc = creal(comp[J]);
					mu_el[J] = 0.0;

					a = alpha[(int)(n_alpha * realc)];
					b = beta[(int)(n_beta * realc)];

					ap = alpha_prime[(int)(n_alpha * realc)];
					bp = beta_prime[(int)(n_beta * realc)];

					strain[1][1] =
						E[1][1] + epsstr11 - epsilon_T[1][1] * b;
					strain[1][2] =
						E[1][2] + epsstr12 - epsilon_T[1][2] * b;
					strain[2][1] =
						E[2][1] + epsstr12 - epsilon_T[2][1] * b;
					strain[2][2] =
						E[2][2] + epsstr22 - epsilon_T[2][2] * b;
					strain[1][3] =
						E[1][3] + epsstr13 - epsilon_T[1][3] * b;
					strain[3][1] =
						E[3][1] + epsstr13 - epsilon_T[3][1] * b;
					strain[3][3] =
						E[3][3] + epsstr33 - epsilon_T[3][3] * b;
					strain[2][3] =
						E[2][3] + epsstr23 - epsilon_T[2][3] * b;
					strain[3][2] =
						E[3][2] + epsstr23 - epsilon_T[3][2] * b;

					temp1 = 0.;
					temp2 = 0.;
					temp3 = 0.;
					temp4 = 0.;
					temp5 = 0.;
					temp6 = 0.;
					for (i3 = 1; i3 < 4; ++i3)
					{
						for (i4 = 1; i4 < 4; ++i4)
						{
							for (i5 = 1; i5 < 4; ++i5)
							{
								for (i6 = 1; i6 < 4; ++i6)
								{
									if (i3 == 1 && i4 == 1)
										temp1 = temp1 + strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);
									if (i3 == 2 && i4 == 2)
										temp2 = temp2 + strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);
									if ((i3 == 1 && i4 == 2) || (i3 == 2 && i4 == 1))
										temp3 = temp3 + strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);
									if ((i3 == 1 && i4 == 3) || (i3 == 3 && i4 == 1))
										temp4 = temp4 + strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);
									if (i3 == 3 && i4 == 3)
										temp5 = temp5 + strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);
									if ((i3 == 2 && i4 == 3) || (i3 == 3 && i4 == 2))
										temp6 = temp6 + strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);

									mu_el[J] = mu_el[J] + 0.5 * ap * DeltaC[i3][i4][i5][i6] * strain[i5][i6] * strain[i3][i4] - bp * epsilon_T[i3][i4] * strain[i5][i6] * (Ceff[i3][i4][i5][i6] + DeltaC[i3][i4][i5][i6] * a);
								}
							}
						}
					}
					//if(INDEX != 0 && INDEX%STEPS ==0){
					if(INDEX%STEPS ==0){
					fprintf(fpstr,"%d %d %d %le %le %le %le %le %le\n",i1,i2,i7,temp1,temp2,temp3,temp4,temp5,temp6);}
				}
			}
		}
		if(INDEX != 0 && INDEX%STEPS ==0){
		fclose(fpstr);
		}
		fftw_execute(planF);
		fftw_execute_dft(planF, comp, comp);

		/* Evolve the composition profile */

		for (i1 = 0; i1 < n_x; ++i1)
		{
			for (i2 = 0; i2 < n_y; ++i2)
			{
				for (i3 = 0; i3 < n_z; ++i3)
				{

					J = i3 + n_z * (i2 + n_y * i1);

					k2 = kx[i1] * kx[i1] + ky[i2] * ky[i2] + kz[i3] * kz[i3];
					k4 = k2 * k2;

					inv_denom = 1.0 + 2.0 * kappa * k4 * delta_t;

					inv_denom = 1.0 / inv_denom;

					comp[J] = inv_denom * (comp[J] - k2 * delta_t * (g[J] + mu_el[J]));
				}
			}
		}

		/* Get the composition back to real space */

		fftw_execute_dft(planB, comp, comp);
		for (i1 = 0; i1 < n_x; ++i1)
		{
			for (i2 = 0; i2 < n_y; ++i2)
			{
				for (i3 = 0; i3 < n_z; ++i3)
				{
					J = i3 + n_z * (i2 + n_y * i1);
					comp[J] = comp[J] * inv_nxnynz;
				}
			}
		}

		/* Once in a while check that the composition is within bounds,
		 * and write the output **/

		if (INDEX != 0 && INDEX % STEPS == 0)
		{
			for (i1 = 0; i1 < n_x; ++i1)
			{
				for (i2 = 0; i2 < n_y; ++i2)
				{
					for (i3 = 0; i3 < n_z; ++i3)
					{
						J = i3 + n_z * (i2 + n_y * i1);
						c[J] = creal(comp[J]);
						if (c[J] <= -1.0 || c[J] >= 2.0)
						{
							printf("The composition goes out of bounds\n");
							printf("Exiting\n");
							exit(0);
						}
					}
				}
			}

			sprintf(NAME, "../output/data/time%dc.vtk",
					(int)(INDEX + t0));
			fpd = fopen(NAME, "w");

			fprintf(fpd, "# vtk DataFile Version 4.2.0\n");
			fprintf(fpd, "Order Parameter data\n");
			fprintf(fpd, "ASCII\n");
			fprintf(fpd, "DATASET STRUCTURED_POINTS\n");
			fprintf(fpd, "DIMENSIONS %d %d %d\n", n_x, n_y, n_z);
			fprintf(fpd, "ORIGIN 0 0 0\n");
			fprintf(fpd, "SPACING 1 1 1\n");
			fprintf(fpd, "POINT_DATA %d\n", n_x * n_y * n_z);
			fprintf(fpd, "SCALARS order_parameter double\n");
			fprintf(fpd, "LOOKUP_TABLE default\n");
			for (i3 = 0; i3 < n_z; ++i3)
			{
				for (i2 = 0; i2 < n_y; ++i2)
				{
					for (i1 = 0; i1 < n_x; ++i1)
					{
						J = i3 + n_z * (i2 + n_y * i1);
						fprintf(fpd, "%le\n", comp[J]);
					}
				}
			}
			fclose(fpd);
		}

		if (INDEX != 0 && INDEX % NOISEADDER == 0 && SNA == 1)
		{
			fprintf(fptmp, "Step %d: Noise added\n", INDEX);
			for (i1 = 0; i1 < n_x; ++i1)
			{
				for (i2 = 0; i2 < n_y; ++i2)
				{
					for (i3 = 0; i3 < n_z; ++i3)
					{
						J = i3 + n_z * (i2 + n_y * i1);
						c[J] = creal(comp[J]);
					}
				}
			}
			add_noise(n_x, n_y, n_z, c, noise_str);
			for (i1 = 0; i1 < n_x; ++i1)
			{
				for (i2 = 0; i2 < n_y; ++i2)
				{
					for (i3 = 0; i3 < n_z; ++i3)
					{
						J = i3 + n_z * (i2 + n_y * i1);
						comp[J] = c[J] + _Complex_I * 0.0;
					}
				}
			}
		}
	}
	fclose(fptmp);

	/* Free the variables and destroy the plans */

	fftw_free(g);
	fftw_free(eps_star11);
	fftw_free(eps_star12);
	fftw_free(eps_star22);
	fftw_free(eps_star13);
	fftw_free(eps_star33);
	fftw_free(eps_star23);
	fftw_free(u1_old);
	fftw_free(u2_old);
	fftw_free(u3_old);
	fftw_free(u1_new);
	fftw_free(u2_new);
	fftw_free(u3_new);
	fftw_free(mu_el);
	fftw_free(eta);

	fftw_destroy_plan(planF);
	fftw_destroy_plan(planB);

	free(c);
	free(kx);
	free(ky);
	free(kz);
	free_dmatrix(Del_sigma_T, 1, 3, 1, 3);
	free_dmatrix(E, 1, 3, 1, 3);
	free_dmatrix(strain, 1, 3, 1, 3);
	free_dvector(Omega11, 0, nxnynz);
	free_dvector(Omega12, 0, nxnynz);
	free_dvector(Omega21, 0, nxnynz);
	free_dvector(Omega22, 0, nxnynz);
	free_dvector(Omega13, 0, nxnynz);
	free_dvector(Omega31, 0, nxnynz);
	free_dvector(Omega33, 0, nxnynz);
	free_dvector(Omega23, 0, nxnynz);
	free_dvector(Omega32, 0, nxnynz);
}
