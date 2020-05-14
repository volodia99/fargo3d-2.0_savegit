#include "fargo3d.h"

void _CondInit() {
  int i,j,k;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *rho;
  real h;
  
  real omega;
  real r, r3;

  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;

  real rhog, rhod;	

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      h = ASPECTRATIO*Ymed(j);
      r = Ymed(j);
      r3 = r*r*r;
      omega = sqrt(G*MSTAR/(r3));
      for (i=NGHX; i<Nx+NGHX; i++) {
	v2[l] = v3[l] = 0.0;
	v1[l] = omega*r;
#ifdef CYLINDRICAL
	rhog = SIGMA0*pow(r/R0,-SIGMASLOPE)*exp(-pow(Zmed(k)/h,2.0)/2.0)/(ZMAX-ZMIN);
#else
	real xi = SIGMASLOPE+1.+FLARINGINDEX;
	real beta = 1.-2*FLARINGINDEX;
	real h = ASPECTRATIO*pow(r/R0,FLARINGINDEX);
	if (FLARINGINDEX == 0.0) {
	  rhog = SIGMA0/sqrt(2.0*M_PI)/(R0*ASPECTRATIO)*pow(r/R0,-xi)* \
	    pow(sin(Zmed(k)),-beta-xi+1./(h*h));
	} else {
	  rhog = SIGMA0/sqrt(2.0*M_PI)/(R0*ASPECTRATIO)*pow(r/R0,-xi)* \
	    pow(sin(Zmed(k)),-xi-beta)*					\
	    exp((1.-pow(sin(Zmed(k)),-2.*FLARINGINDEX))/2./FLARINGINDEX/(h*h));
	}
	rhod = rhog*EPSILON;
	
	if (Fluidtype == GAS) {
	  rho[l] = rhog;
	  v1[l] *= sqrt(pow(sin(Zmed(k)),-2.*FLARINGINDEX)-(beta+xi)*h*h);
	  v1[l] -= OMEGAFRAME*r*sin(Zmed(k));
	}
	if (Fluidtype == DUST) {
	  rho[l] = rhod;
	  v1[l] -= OMEGAFRAME*r*sin(Zmed(k));
	}
#endif
	if (Fluidtype == GAS) {
#ifdef ISOTHERMAL
	  e[l] = h*sqrt(G*MSTAR/r);
#else
	  e[l] = rho[l]*h*h*G*MSTAR/r/(GAMMA-1.0);
#endif
	}
	if (Fluidtype == DUST) {
	  e[l] = 0.0;
	}
      }
    }
  }
}

void CondInit() {

  int id_gas = 0;
  int feedback = YES;
  //We first create the gaseous fluid and store it in the array Fluids[]
  Fluids[id_gas] = CreateFluid("gas",GAS);

  //We now select the fluid
  SelectFluid(id_gas);

  //and fill its fields
  _CondInit();

  //We repeat the process for the dust fluids
  char dust_name[MAXNAMELENGTH];
  int id_dust;

  for(id_dust = 1; id_dust<NFLUIDS; id_dust++) {
    sprintf(dust_name,"dust%d",id_dust); //We assign different names to the dust fluids

    Fluids[id_dust]  = CreateFluid(dust_name, DUST);
    SelectFluid(id_dust);
    _CondInit();

  }

  /*We now fill the collision matrix (Feedback from dust included)
   Note: ColRate() moves the collision matrix to the device.
   If feedback=NO, gas does not feel the drag force.*/
  
  ColRate(INVSTOKES1, id_gas, 1, feedback);
  ColRate(INVSTOKES2, id_gas, 2, feedback);
  ColRate(INVSTOKES3, id_gas, 3, feedback);

}
