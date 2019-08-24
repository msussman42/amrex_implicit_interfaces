//   Stage| Description
//   -----+--------------------------------------
//    -1  +-initialization (0)
//     0  | Liquid cooling
//     1  +-Nucleation
//     2  | Recalescence in progress
//     3  +-Recalescence finished, frost starting (1)
//     4  | Frost  
//     5  +-Regular freezing starts (2) 
//     6  | Regular freezing in progress
//     6  | Solid cooling in progress
//

// (0): Change the ice surface tensions to values similar to
//      substarte. It's done in the calling function.
//      Noting to be done here.
// (1): Solid sruface tensions are set to be the same as the substrate
//      Surface tension to emulate ice growing out of the substarte,
//      T_w->T_f, and Z_w->Z_mush, Z={Cp, mu, k}
// (2): Ice surface tensions are set to be the original value

// Right now stage 1, 3, and 5 are internal only and cannot seen by
// the calling function.

#define INVALID_TIME -1000

void multistage_freezing(
	
 int nmat, // number of materials
 int nten, // number of surface tension coef
 Array<int>& recal_material,  // 0=no 1=static 2=dynamic 
   // treatment, size=nmat
 Array<Real>& freezing_state_old, 
  // nmat*(0-stage,1-touch_time,2-nucleation_time,3-frost_time, 
  // 4-reg_freezing_time, 5-mushy_frac). size=6*nmat
 Array<Real>& freezing_state_new, 
  // nmat*(0-stage,1-touch_time,2-nucleation_time, 3-frost_time,
  // 4-reg_freezing_time, 5-mushy_frac). size=6*nmat
 Array<Real>& material_state, 
  // nmat*(0-T_average, 1-T_top, 2-dist_top, 3-dist_bottom, 
  // 4-touched, 5-centroid, 6-mass). size=nmat* (6+BL_SPACEDIM)
 Array<Real>& exper_var, 
 // nmat*(0-trigger_temperature, 1-recal_speed, 2-frost_period). size=3*nmat
	Array<Real>& density;   // size=nmat
 Array<Real>& latent_heat, // size=2*nten
 Array<Real>& heat_conduct, // size=nmat
 Array<Real>& spec_heat, // size=nmat
 Array<Real>& visc_coef, // size=nmat
 Array<Real>& surf_ten, // size=nten
 Real time)
{
	Real drop_height, recal_duration
		// Sanity checks
  for (int im=0; im<nmat; im++) {
			// touched==1 -> touch_time <= time
			if((material_state[im*(6+BLSpaceDIM)+4]==1.0 &&
							(freezing_state_old[im*nmat+1]>time))) {
				BoxLib::Error("multistage_freezing: invalid touch time!");
			}
		}

	// Copy old stage data to new stage data...
	for (int i=0; i<freezing_stage_old.size(); i++) {
		freezing_stage_new[i]=freezing_stage_old[i];
	}
	// and make nessecary changes for each material
	for (int im=0; im<nmat; im++) {
		for (int im_opp=im+1;im_opp<=nmat;im_opp++) {
			for (int ireverse=0;ireverse<=1;ireverse++) {
				int iten;
				get_iten_cpp(im,im_opp,iten,nmat);
				int indexEXP=iten+ireverse*nten-1;
				Real LL=latent_heat[indexEXP];
				int im_source;
				int im_dest;
				if (LL!=0.0) {

					if (ireverse==0) {
						im_source=im;
						im_dest=im_opp;
					} else if (ireverse==1) {
						im_source=im_opp;
						im_dest=im;
					}

					if ((recal_material==1) || (recal_material==2)) {

						////// stage -1: initialization
						if (freezing_stage_new[im_source*6+0]==-1.0) {
							freezing_stage_new[im_source*6+0]=0.0;
						} // end of stage -1

						////// stage 0: liquid cooling before impact
						if (freezing_stage_new[im_source*6+0]==0.0) {
							// if no contact
							if (material_state[im_source*(6+BL_SPACEDIM)+4]==0){
								// do nothing
								freezing_stage_new[im_source*6+1]=INVALID_TIME;
							}
							else { // if has contact
								// record the contact time if not recorded yet
								if (freezing_stage_old[im_source*6+1]==INVALID_TIME){
									freezing_stage_new[im_source*6+1]=time;
								}	
								// check for the nucleation critrria. Necleation starts at
								// specific top temperature that comes from experiment for
								// static droplet, or an deduced average droplet temperature for
								// the impact case.
								if (((recal_material==1) && 
													(freezing_stage_new[im_source*6+1]!=INVALID_TIME) &&
													(material_state[im_source*(6+BL_SPACEDIM)+1]<=exper_var[im_source*2+0])) ||
												((recal_material==2) && 
													(freezing_stage_new[im_source*6+1]!=INVALID_TIME) &&
													(material_state[im_source*(6+BL_SPACEDIM)+0]<=exper_var[im_source*2+0])))	{
									// tansit to nucleation stage
									freezing_stage_new[im_source*6+0]=1.0;
									// set nuclation_time
									freezing_stage_new[im_source*6+2]=time;
								}
							}
						} // end of stage 0

						////// stage 1: nucleation
						if (freezing_stage_new[im_source*6+0]==1.0) {
							// tansit to recalescence stage
							freezing_stage_new[im_source*6+0]=2.0;
						} // end of stage 1

						////// stage 2: recalescence
						if (freezing_stage_new[im_source*6+0]==2.0) {
							drop_height = material_state[im_source*(6+BL_SPACEDIM)+2] -
								material_state[im_source*(6+BL_SPACEDIM)+3];
							// recal_duration = drop_height / recal_speed
							recal_duration = drop_height / exper_var[im_source*2+1];
							// update end of recal stage
							freezing_stage_new[im_source*6+3]=freezing_stage_new[im_source*6+2] 
								+ recal_duration;
							// check if passed the end of recal stage
							if(time>=freezing_stage_new[im_source*6+3]){
								freezing_stage_new[im_source*6+0]=3.0;
							}
						} // end of stage 2
				
						////// stage 3: end of recal, initiating frost
						if (freezing_stage_new[im_source*6+0]==3.0) {
							// tansit to frost stage
							freezing_stage_new[im_source*6+0]=4.0;
							freezing_stage_new[im_source*6+4]=freezing_stage_new[im_source*6+3] 
								+ exper_var[im_source*2+2] ;
							// droplet temperature is changed to T_m (in the calling
							// function), and water properties are changed to the mushy
							// material properties (in the calling function)
						} // end of stage 3

						////// stage 4: frost
						if (freezing_stage_new[im_source*6+0]==4.0) {
							// wait to pass the frost
							if (time>=freezing_stage_new[im_source*6+4]){
								freezing_stage_new[im_source*6+0]=5.0;
							}
						} // end of stage 4

						////// stage 5: initiation of regular freezing
						if (freezing_stage_new[im_source*6+0]==5.0) {
							// The surface tensions are changed for solid material to
							// the physical values (in the calling function). 
							freezing_stage_new[im_source*6+0]=6.0;
						} // end of stage 5

						////// stage 6: regular freezung or solid cooling
						if (freezing_stage_new[im_source*6+0]==6.0) {
							// do nothing
						} // end of stage 6
					} // recal_material 1 or 2

					else if (recal_material==0){
						// do nothing
					}
					else {
						BoxLib::Error("multistage_freezing: invalid freezing model!");
					} // recal_material
				} // LL!=0
			} // ireverse
		} // im_opp
	} // im
} // multistage_freezing
