/** My header file for function declarations **/

extern void evolve_diffusional(int n_x, int n_y, int n_z, double delta_x, double delta_y, double delta_z, 
double kappa, double A, double delta_t, int time_steps, double t0,
fftw_complex *comp, double ave_comp, int STEPS, int plan_indicator, 
int SNA, int NOISEADDER, double noise_str);

extern void generate_elasticity_tensor(double c11m, double c12m, 
double c44m, double c11p,double c12p, double c44p, 
double Ceff[4][4][4][4], double DeltaC[4][4][4][4],
int alpha_interptype, int prob_type, double ave_comp);

extern void calculate_sigma_T(double Ceff[4][4][4][4], 
double eps0_11, double eps0_22, double eps0_33, double eps0_12,
double eps0_13, double eps0_23, double **sigma_T, 
double **epsilon_T);

extern void calculate_Omega(int n_x, int n_y, int n_z, int half_nx, int half_ny, int half_nz,
double *kx, double *ky, double *kz, double Ceff[4][4][4][4], double *Omega11, 
double *Omega12, double*Omega21, double *Omega22,
double *Omega13, double*Omega31, double *Omega33,
double *Omega23, double*Omega32);



extern void evolve_inhom_constrained(int n_x, int n_y, int n_z, double delta_x, 
double delta_y, double delta_z, double kappa, double A, double delta_t, int
time_steps, double t0, fftw_complex *comp, double Ceff[4][4][4][4],
double DeltaC[4][4][4][4], double **sigma_T, double **epsilon_T, 
double ave_comp, double **sig_app, int STEPS, int MAXITR, 
double MAXERR, int n_alpha, int n_beta, double *alpha, double *beta,
double *alpha_prime, double *beta_prime, int plan_indicator,
int SNA, int NOISEADDER, double noise_str);

extern void calculate_Del_sigma_T(double DeltaC[4][4][4][4], 
double **epsilon_T, double **Del_sigma_T);

extern void calculate_uzero(int n_x, int n_y, int n_z, int half_nx, int half_ny, int half_nz,
double delta_kx, double delta_ky, double delta_kz, double ave_comp, fftw_complex *comp, 
double Ceff[4][4][4][4], double **sigma_zero, fftw_complex *u1_old, 
fftw_complex *u2_old, fftw_complex *u3_old, int n_beta, double *beta, double *Omega11,
double *Omega12, double *Omega21, double *Omega22, double *Omega13, double *Omega31, double *Omega33,
 double *Omega23, double *Omega32, fftw_plan planF);

extern void refine_u(int n_x, int n_y, int n_z, int half_nx, int half_ny, int half_nz,
double delta_x, double delta_y, double delta_z, double **epsilon_T,
double delta_kx, double delta_ky, double delta_kz, fftw_complex *comp,
double ave_comp, double Ceff[4][4][4][4], double DeltaC[4][4][4][4], 
double **sigma_zero, fftw_complex *u1_old, fftw_complex *u2_old, fftw_complex *u3_old, 
fftw_complex *u1_new, fftw_complex *u2_new, fftw_complex *u3_new, double **Del_sigma_T, 
double **sig_app, double **E, int MAXITS, double MAXERR,
int n_alpha, int n_beta, double *alpha, double *beta,
double *Omega11, double *Omega12, double *Omega21,
 double *Omega22, double *Omega13, double *Omega31, double *Omega33, double *Omega23, double *Omega32,
fftw_plan planF, fftw_plan planB);

extern void calculate_hom_strain(int n_x, int n_y, int n_z, 
double delta_x, double delta_y, double delta_z, double Ceff[4][4][4][4], 
double DeltaC[4][4][4][4], double S[4][4][4][4], 
double **sig_zero_exp, fftw_complex *comp, double ave_comp,
fftw_complex *eps_star11, fftw_complex *eps_star12,
fftw_complex *eps_star22, fftw_complex *eps_star13, fftw_complex *eps_star33,
fftw_complex *eps_star23, double **sig_app, double **E, 
int n_alpha, double *alpha);

extern void calculate_alpha(int n_alpha, double *alpha, 
int alpha_interptype, double steepness_factor, double ave_comp);

extern void calculate_beta(int n_beta, double *beta,
int beta_interptype, double steepness_factor, double ave_comp);

extern void calculate_alpha_prime(int n_alpha, double *alpha_prime,
int alpha_interptype, double steepness_factor);

extern void calculate_beta_prime(int n_beta, double *beta_prime,
int beta_interptype, double steepness_factor);

extern void calculate_S_exp(int n_x, int n_y, int n_z, double delta_x, 
double delta_y, double delta_z, fftw_complex *comp, double ave_comp, 
double Ceff[4][4][4][4], double DeltaC[4][4][4][4], 
double S[4][4][4][4], int n_alpha, double *alpha);

extern void calculate_sig_zero_exp(int n_x, int n_y, int n_z, double delta_x,
double delta_y, double delta_z, fftw_complex *comp, double ave_comp, 
double Ceff[4][4][4][4], double DeltaC[4][4][4][4], 
double **epsilon_T, double **sig_zero_exp, int n_alpha,
int n_beta, double *alpha, double *beta);

extern void write_ppfiles(int n_x, int n_y, int n_z, double delta_t1, 
int time_steps1, double t0, double delta_t2, int time_steps2,
double t1, double delta_t3, int time_steps3, double t2, int STEPS);

extern void add_noise(int n_x, int n_y, int n_z, double *comp, 
double noise_str);

extern void calculate_S(double Ceff[4][4][4][4], double S[4][4][4][4]);

extern void calculate_aleph(double c11m, double c12m, double c44m,
double c11p, double c12p, double c44p, double **sig_app,
double **epsilon_T);
