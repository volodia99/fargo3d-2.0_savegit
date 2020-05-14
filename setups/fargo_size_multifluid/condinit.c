#include "fargo3d.h"

void _CondInit() {

  int i,j,k;
  real r, omega;

  real *rho  = Density->field_cpu;
  real *cs   = Energy->field_cpu;
  real *vphi = Vx->field_cpu;
  real *vr   = Vy->field_cpu;

  real rhog, rhod;
  real vk;

  i = j = k = 0;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {

	r     = Ymed(j);
	omega = sqrt(G*MSTAR/r/r/r);                       //Keplerian frequency
	rhog  = SIGMA0*pow(r/R0,-SIGMASLOPE);              //Gas surface density
  rhod  = rhog*EPSILON;                              //Dust surface density

	if (Fluidtype == GAS) {
	  rho[l]   = rhog;
	  vphi[l]  = omega*r*sqrt(1.0 + pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
				  (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));
	  vr[l]    = 0.0;
	  cs[l]    = ASPECTRATIO*pow(r/R0,FLARINGINDEX)*omega*r;
    // printf("%f,", rho[l]);
}

	if (Fluidtype == DUST) {
	  rho[l]  = rhod;
	  vphi[l] = omega*r;
	  vr[l]   = 0.0;
	  cs[l]   = 0.0;
	}

	vphi[l] -= OMEGAFRAME*r;

      }
    }
  }
}

void CondInit() {
  int i,j,k;
  i = j = k = 0;
  real sumsize = 0.0;
  int id_gas = 0;
  int feedback = YES;
  // int feedback = NO;

  real *rho[NFLUIDS];
  real size[NFLUIDS];

  //We first create the gaseous fluid and store it in the array Fluids[]
  Fluids[id_gas] = CreateFluid("gas",GAS);

  //We now select the fluid
  SelectFluid(id_gas);

  //and fill its fields
  _CondInit();

  rho[id_gas]  = Fluids[id_gas]->Density->field_cpu;

  //We repeat the process for the dust fluids
  char dust_name[MAXNAMELENGTH];
  int id_dust;

  #ifdef EPSTEINDRAG
  for(id_dust = 1; id_dust<NFLUIDS; id_dust++) {
    if (NFLUIDS==2){
      size[id_dust] = DUSTSIZEMAX;
    }
    else {
      size[id_dust] = DUSTSIZEMIN*exp(((float)(id_dust-1)/(NFLUIDS-2))*log(DUSTSIZEMAX/DUSTSIZEMIN));
    }
    sumsize+=pow(size[id_dust],4.0-DUSTSIZESLOPE);
  }
  #endif

  for(id_dust = 1; id_dust<NFLUIDS; id_dust++) {
    sprintf(dust_name,"dust%d",id_dust); //We assign different names to the dust fluids

    Fluids[id_dust]  = CreateFluid(dust_name, DUST);
    SelectFluid(id_dust);
    _CondInit();

    rho[id_dust]  = Fluids[id_dust]->Density->field_cpu;

    #ifdef EPSTEINDRAG
    for (k=0; k<Nz+2*NGHZ; k++) {
      for (j=0; j<Ny+2*NGHY; j++) {
        for (i=0; i<Nx+2*NGHX; i++) {
          rho[id_dust][l]*=pow(size[id_dust],4.0-DUSTSIZESLOPE)/sumsize;
        }
      }
    }

    ColRate(2.*FACTORUNITMASS*1.989e33/M_PI/INTERNALDENSITY/size[id_dust]/100/(pow(FACTORUNITLENGTH*1.49597870e13,2.0)), id_gas, id_dust, feedback);
    #endif

  }

  /*We now fill the collision matrix (Feedback from dust included)
   Note: ColRate() moves the collision matrix to the device.
   If feedback=NO, gas does not feel the drag force.*/

#ifndef EPSTEINDRAG
  ColRate(INVSTOKES1, id_gas, 1, feedback);
  ColRate(INVSTOKES2, id_gas, 2, feedback);
  ColRate(INVSTOKES3, id_gas, 3, feedback);
  ColRate(INVSTOKES4, id_gas, 4, feedback);
  ColRate(INVSTOKES5, id_gas, 5, feedback);
  ColRate(INVSTOKES6, id_gas, 6, feedback);
  ColRate(INVSTOKES7, id_gas, 7, feedback);
#endif
}
