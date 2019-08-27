//in general one does the following:
//1. declare grid as global: Grid* gridptr = NULL;
//2. define grid: MLS_DEFGRID(m, n, x_lo, x_hi, y_lo, y_hi) 
//3. initialize grid level set: MLS_INITLS(LS,gridptr);
//4. start computational loop: at any time step do
//MLS_READVEL(U,ngrow,gridptr);
//MLS_ADVECT(gridptr,dt);
//MLS_REINIT(gridpts,dt,totaltime,dx);
//MLS_WRITELS(LS,gridptr);

void MLS_DEFGRID(int m, int n, float x_lo, float x_hi, float y_lo, float y_hi){ 
	gridptr = new Grid(m,n,x_lo,x_hi,y_lo,y_hi);//build grid
}
void MLS_INITLS(FArrayBox& LS, Grid& thegrid, float xloM, float xhiM, float yloM, float yhiM, int lo, int hi, float dxM, float dyM){
	//lo = lo(1) and hi = hi(1) 
	int m = thegrid->m, n = thegrid->n, i, j;
	float dxV = thegrid->dx, dyV = thegrid->dy, xloV = thegrid->x_lo, yloV = thegrid->y_lo;
	int ilo = int((xloM-xloV)/dxV),jlo = int((yloM-yloV)/dyV),ihi = int((xhiM-xloV)/dxV),jhi = int((yhiM-yloV)/dyV);
	//now we have to initialize MLS data on [ilo ihi]x[jlo jhi]
	int Nx0 = int((xloV+dxV*ilo-xloM)/dxM), Ny0 = int((yloV+dyV*jlo-yloM)/dyM), Nx = int(dxV/dxM), Ny = int(dyV/dyM), iM, jM;
	float NE=0,NW=0,SW=0,SE=0;
	int ngrow = 1, n_hor = hi-lo+1+2*ngrow;
	//write level set values
	for (j = 0; j < (jhi-jlo+1); j++)
		for (i = 0; i < (ihi-ilo+1); i++){
			iM = Nx0+i*Nx; jM = Ny0+j*Ny;
			SW = LS->dataPtr(0)[iM + n_hor*jM];
			SE = LS->dataPtr(0)[iM+1 + n_hor*jM];
			NW = LS->dataPtr(0)[iM + n_hor*(jM+1)];
			NE = LS->dataPtr(0)[iM+1 + n_hor*(jM+1)];
			thegrid->phi[i+ilo][j+jlo] = (SW+SE+NW+NE)/4.0;
		}
}
void MLS_READVEL(FArrayBox& U, Grid& thegrid, float xloM, float xhiM, float yloM, float yhiM, int lo, int hi, float dxM, float dyM){
	//lo = lo(1) and hi = hi(1) 
	int m = thegrid->m, n = thegrid->n, i, j;
	float dxV = thegrid->dx, dyV = thegrid->dy, xloV = thegrid->x_lo, yloV = thegrid->y_lo;
	int ilo = int((xloM-xloV)/dxV),jlo = int((yloM-yloV)/dyV),ihi = int((xhiM-xloV)/dxV),jhi = int((yhiM-yloV)/dyV);
	//now we have to write velocity data on the [ilo ihi]x[jlo jhi] subset of the MLS grid
	int Nx0 = int((xloV+dxV*ilo-xloM)/dxM), Ny0 = int((yloV+dyV*jlo-yloM)/dyM), Nx = int(dxV/dxM), Ny = int(dyV/dyM), iM, jM;
	float NE=0,NW=0,SW=0,SE=0;
	int ngrow = 1, n_hor = hi-lo+1+2*ngrow;
	//write u velocity values
	for (j = 0; j < (jhi-jlo+1); j++)
		for (i = 0; i < (ihi-ilo+1); i++){
			iM = Nx0+i*Nx; jM = Ny0+j*Ny;//the 'hot' index inside U(0)
			SW = U->dataPtr(0)[iM + n_hor*jM];
			SE = U->dataPtr(0)[iM+1 + n_hor*jM];
			NW = U->dataPtr(0)[iM + n_hor*(jM+1)];
			NE = U->dataPtr(0)[iM+1 + n_hor*(jM+1)];
			thegrid->u[i+ilo][j+jlo] = (SW+SE+NW+NE)/4.0;
			SW = U->dataPtr(1)[iM + n_hor*jM];
			SE = U->dataPtr(1)[iM+1 + n_hor*jM];
			NW = U->dataPtr(1)[iM + n_hor*(jM+1)];
			NE = U->dataPtr(1)[iM+1 + n_hor*(jM+1)];
			thegrid->v[i+ilo][j+jlo] = (SW+SE+NW+NE)/4.0;
		}
}
void MLS_ADVECT(Grid& thegrid, float dt){
	thegrid->update_levelsetSL(dt);
	thegrid->AdvectParticles(dt, 0, 1);
}
void MLS_REINIT(Grid& thegrid, float dt, float totaltime){
	int reinitflag = fmod(totaltime/dt,20);
	float dxv = thegrid->dx; 
	thegrid->Redistance1();
	if (reinitflag==0){
		thegrid->CheckReinit(6,dxv/2.0);
		thegrid->FindZeroLevelCurve(1,1, -dxv, 0);
		thegrid->FindZeroLevelCurve(1,1, dxv, 0);
		thegrid->RedistributeParticlesGrid(0.0f);
	}
}
void MLS_WRITELS(FArrayBox& LS, Grid& thegrid, float xloM, float yloM, int lo1, int lo2, int hi1, int hi2, float dxM, float dyM){

	int m = thegrid->m, n = thegrid->n, i, j, ilo ,jlo, ngrow = 1, n_hor = hi1-lo1+1+2*ngrow;
	float dxV = thegrid->dx, dyV = thegrid->dy, xloV = thegrid->x_lo, yloV = thegrid->y_lo,ctrx,ctry,NE,NW,SW,SE,val,s,t;

	for (j = 1; j <= hi2-lo2+1; j++)//we skip the 0th row, as it's a ghost
		for (i = 1; i <= hi1-lo1+1; i++){//same with columns
			ctrx = xloM+(i-0.5)*dxM;
			ctry = yloM+(j-0.5)*dyM;
			ilo = int((ctrx-xloV)/dxV);
			jlo = int((ctry-yloV)/dyV);
			SW=thegrid->phi[ilo][jlo];
			SE=thegrid->phi[ilo+1][jlo];
			NW=thegrid->phi[ilo][jlo+1];
			NE=thegrid->phi[ilo+1][jlo+1];
			s = (ctrx-xloV)/dxV-ilo;
			t = (ctry-yloV)/dyV-jlo;
			val = (1-t)*((1-s)*SW+s*SE)+t*((1-s)*NW+s*NE);
			LS->dataPtr(0)[i+n_hor*j] = val;
		}
}