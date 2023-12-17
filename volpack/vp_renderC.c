/*
 * vp_renderC.c
 *
 * Function to render classified volumes.
 *
 * Copyright (c) 1994 The Board of Trustees of The Leland Stanford
 * Junior University.  All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided
 * that the above copyright notice and this permission notice appear in
 * all copies of this software and that you do not sell the software.
 * Commercial licensing is available by contacting the author.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
 * WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author:
 *    Phil Lacroute
 *    Computer Systems Laboratory
 *    Electrical Engineering Dept.
 *    Stanford University
 *
 * This version of the VolPack library is a modified version of the
 * original available from Stanford University.  Modifications were
 * made by the Center for Computational Sciences and Engineering,
 * Lawrence Berkeley National Laboratory.  Modifications to VolPack
 * (c) 2000 The Regents of the University of California (through
 * E.O. Lawrence Berkeley National Laboratory), subject to approval by
 * the U.S. Department of Energy.  Your use of this software is under
 * license from Stanford University with respect to the underlying
 * original VolPack code (see copyright notice and permission above)
 * and The Regents of the University of California with respect to
 * modifications thereto (see AmrVis license.txt file for applicable
 * license terms).  Contact Berkeley Lab's Center for Computational
 * Sciences and Engineering at webmaster@mothra.lbl.gov or Berkeley
 * Lab's Technology Transfer Department at TTD@lbl.gov for questions
 * or to receive more information.
 *
 */

/*
 * $Date: 2000-10-02 18:12:03 $
 * $Revision: 1.3 $
 */

#include "vp_global.h"

#define COMP_AC1PB_FUNC		VPCompAC1PB
extern void VPCompAC1PB(
			vpContext *vpc,			/* context */
			int icount,			/* slice size */
			int jcount,
			int k,				/* slice number */
			double slice_depth_cueing_dbl,	/* depth cueing factor for slice */
			GrayIntPixel *intimage,		/* intermediate image pixels */
			double weightTLdbl,		/* resampling weights */
			double weightBLdbl,
			double weightTRdbl,
			double weightBRdbl,
			unsigned char *run_lengths,	/* run lengths for slice */
			void *voxel_data		/* voxel data for slice */
			);

#define COMP_AC3PB_FUNC		VPCompAC3PB
extern void VPCompAC3PB(
			vpContext *vpc,			/* context */
			int icount,			/* slice size */
			int jcount,
			int k,				/* slice number */
			double slice_depth_cueing_dbl,	/* depth cueing factor for slice */
			RGBIntPixel *intimage,		/* intermediate image pixels */
			double weightTLdbl,		/* resampling weights */
			double weightBLdbl,
			double weightTRdbl,
			double weightBRdbl,
			unsigned char *run_lengths,	/* run lengths for slice */
			void *voxel_data		/* voxel data for slice */
			);

#ifdef COMP_AC11B
#define COMP_AC11B_FUNC		VPCompAC11B
extern void VPCompAC11B(
			vpContext *vpc,			/* context */
			int icount,			/* slice size */
			int jcount,
			int k,				/* slice number */
			double slice_depth_cueing_dbl,	/* depth cueing factor for slice */
			GrayIntPixel *intimage,		/* intermediate image pixels */
			double weightTLdbl,		/* resampling weights */
			double weightBLdbl,
			double weightTRdbl,
			double weightBRdbl,
			unsigned char *run_lengths,	/* run lengths for slice */
			void *voxel_data		/* voxel data for slice */
			);
#else
#define COMP_AC11B_FUNC		VPCompAC1NB
#endif

#ifdef COMP_AC31B
#define COMP_AC31B_FUNC		VPCompAC31B
extern void VPCompAC31B(
			vpContext *vpc,			/* context */
			int icount,			/* slice size */
			int jcount,
			int k,				/* slice number */
			double slice_depth_cueing_dbl,	/* depth cueing factor for slice */
			RGBIntPixel *intimage,		/* intermediate image pixels */
			double weightTLdbl,		/* resampling weights */
			double weightBLdbl,
			double weightTRdbl,
			double weightBRdbl,
			unsigned char *run_lengths,	/* run lengths for slice */
			void *voxel_data		/* voxel data for slice */
			);
#else
#define COMP_AC31B_FUNC		VPCompAC3NB
#endif

#ifdef COMP_AC12B
#define COMP_AC12B_FUNC		VPCompAC12B
extern void VPCompAC12B();
#else
#define COMP_AC12B_FUNC		VPCompAC1NB
#endif

#ifdef COMP_AC32B
#define COMP_AC32B_FUNC		VPCompAC32B
extern void VPCompAC32B(
			vpContext *vpc,			/* context */
			int icount,			/* slice size */
			int jcount,
			int k,				/* slice number */
			double slice_depth_cueing_dbl,	/* depth cueing factor for slice */
			RGBIntPixel *intimage,		/* intermediate image pixels */
			double weightTLdbl,		/* resampling weights */
			double weightBLdbl,
			double weightTRdbl,
			double weightBRdbl,
			unsigned char *run_lengths,	/* run lengths for slice */
			void *voxel_data		/* voxel data for slice */
			);
#else
#define COMP_AC32B_FUNC		VPCompAC3NB
#endif

#define COMP_AC1NB_FUNC		VPCompAC1NB
extern void VPCompAC1NB(
			vpContext *vpc,			/* context */
			int icount,			/* slice size */
			int jcount,
			int k,				/* slice number */
			double slice_depth_cueing_dbl,	/* depth cueing factor for slice */
			GrayIntPixel *intimage,		/* intermediate image pixels */
			double weightTLdbl,		/* resampling weights */
			double weightBLdbl,
			double weightTRdbl,
			double weightBRdbl,
			unsigned char *run_lengths,	/* run lengths for slice */
			void *voxel_data		/* voxel data for slice */
			);

#define COMP_AC3NB_FUNC		VPCompAC3NB
extern void VPCompAC3NB(
			vpContext *vpc,			/* context */
			int icount,			/* slice size */
			int jcount,
			int k,				/* slice number */
			double slice_depth_cueing_dbl,	/* depth cueing factor for slice */
			RGBIntPixel *intimage,		/* intermediate image pixels */
			double weightTLdbl,		/* resampling weights */
			double weightBLdbl,
			double weightTRdbl,
			double weightBRdbl,
			unsigned char *run_lengths,	/* run lengths for slice */
			void *voxel_data		/* voxel data for slice */
			);

#define COMP_AC1PS_FUNC		VPCompAC1PB

#define COMP_AC3PS_FUNC		VPCompAC3PB

#ifdef COMP_AC11S
#define COMP_AC11S_FUNC		VPCompAC11S
extern void VPCompAC11S();
#else
#define COMP_AC11S_FUNC		VPCompAC1NS
#endif

#ifdef COMP_AC31S
#define COMP_AC31S_FUNC		VPCompAC31S
extern void VPCompAC31S();
#else
#define COMP_AC31S_FUNC		VPCompAC3NS
#endif

#ifdef COMP_AC12S
#define COMP_AC12S_FUNC		VPCompAC12S
extern void VPCompAC12S();
#else
#define COMP_AC12S_FUNC		VPCompAC1NS
#endif

#ifdef COMP_AC32S
#define COMP_AC32S_FUNC		VPCompAC32S
extern void VPCompAC32S();
#else
#define COMP_AC32S_FUNC		VPCompAC3NS
#endif

#define COMP_AC1NS_FUNC		VPCompAC1NS
extern void VPCompAC1NS(
			vpContext *vpc,			/* context */
			int icount,			/* slice size */
			int jcount,
			int k,				/* slice number */
			double slice_depth_cueing_dbl,	/* depth cueing factor for slice */
			GrayIntPixel *intimage,		/* intermediate image pixels */
			double weightTLdbl,		/* resampling weights */
			double weightBLdbl,
			double weightTRdbl,
			double weightBRdbl,
			unsigned char *run_lengths,	/* run lengths for slice */
			void *voxel_data,		/* voxel data for slice */
			GrayIntPixel *shadow_buffer
			);

#define COMP_AC3NS_FUNC		VPCompAC3NS
extern void VPCompAC3NS(
			vpContext *vpc,			/* context */
			int icount,			/* slice size */
			int jcount,
			int k,				/* slice number */
			double slice_depth_cueing_dbl,	/* depth cueing factor for slice */
			RGBIntPixel *intimage,		/* intermediate image pixels */
			double weightTLdbl,		/* resampling weights */
			double weightBLdbl,
			double weightTRdbl,
			double weightBRdbl,
			unsigned char *run_lengths,	/* run lengths for slice */
			void *voxel_data,		/* voxel data for slice */
			GrayIntPixel *shadow_buffer
			);

#ifdef INDEX_VOLUME
extern void VPCompAI11B();
#endif

#define SHADOWS_OFF		0
#define SHADOWS_ON		1
#define SHADOW_OPTS		2

#define MATERIAL_CALLBACK	0
#define MATERIAL_ONE		1
#define MATERIAL_TWO		2
#define MATERIAL_MORE		3
#define MATERIAL_OPTS		4

#define COLOR_GRAY		0
#define COLOR_RGB		1
#define COLOR_OPTS		2

static void (*AffineProcTable[SHADOW_OPTS][MATERIAL_OPTS][COLOR_OPTS])(vpContext*,...) = {
    {
	{ (void(*)(vpContext*,...))COMP_AC1PB_FUNC, (void(*)(vpContext*,...))COMP_AC3PB_FUNC },
	{ (void(*)(vpContext*,...))COMP_AC11B_FUNC, (void(*)(vpContext*,...))COMP_AC31B_FUNC },
	{ (void(*)(vpContext*,...))COMP_AC12B_FUNC, (void(*)(vpContext*,...))COMP_AC32B_FUNC },
	{ (void(*)(vpContext*,...))COMP_AC1NB_FUNC, (void(*)(vpContext*,...))COMP_AC3NB_FUNC }
    },
    {
	{ (void(*)(vpContext*,...))COMP_AC1PS_FUNC, (void(*)(vpContext*,...))COMP_AC3PS_FUNC },
	{ (void(*)(vpContext*,...))COMP_AC11S_FUNC, (void(*)(vpContext*,...))COMP_AC31S_FUNC },
	{ (void(*)(vpContext*,...))COMP_AC12S_FUNC, (void(*)(vpContext*,...))COMP_AC32S_FUNC },
	{ (void(*)(vpContext*,...))COMP_AC1NS_FUNC, (void(*)(vpContext*,...))COMP_AC3NS_FUNC }
    }
};

/*
 * vpRenderClassifiedVolume
 *
 * Render a classified, run-length encoded volume using the shear-warp
 * algorithm.
 */

vpResult
vpRenderClassifiedVolume(vpc)
vpContext *vpc;
{
    int retcode;
    void (*composite_func)(vpContext*,...);
    int shadow_option, material_option, color_option;
#ifdef INDEX_VOLUME
    int have_index;
#endif

    /* check for errors and initialize */
    if ((retcode = VPCheckShader(vpc)) != VP_OK)
	return(retcode);
    if ((retcode = VPFactorView(vpc)) != VP_OK)
	return(retcode);
    if ((retcode = VPCheckClassifiedVolume(vpc, vpc->best_view_axis)) != VP_OK)
	return(retcode);
    if ((retcode = VPCheckImage(vpc)) != VP_OK)
	return(retcode);
    if ((retcode = VPCheckShadows(vpc)) != VP_OK)
	return(retcode);

#ifdef INDEX_VOLUME
    switch (vpc->best_view_axis) {
    case VP_X_AXIS:
	if (vpc->rle_x->voxel_index != NULL)
	    have_index = 1;
	else
	    have_index = 0;
	break;
    case VP_Y_AXIS:
	if (vpc->rle_y->voxel_index != NULL)
	    have_index = 1;
	else
	    have_index = 0;
	break;
    case VP_Z_AXIS:
	if (vpc->rle_z->voxel_index != NULL)
	    have_index = 1;
	else
	    have_index = 0;
	break;
    default:
	VPBug("invalid viewing axis in vpRenderClassifiedVolume");
    }
#endif /* INDEX_VOLUME */

    Debug((vpc, VPDEBUG_RENDER, "Algorithm: affine RLEvolume\n"));

    /* determine which options are enabled */
    if (vpc->enable_shadows)
	shadow_option = SHADOWS_ON;
    else
	shadow_option = SHADOWS_OFF;
    if (vpc->shading_mode == CALLBACK_SHADER)
	material_option = MATERIAL_CALLBACK;
    else if (vpc->num_materials == 1)
	material_option = MATERIAL_ONE;
    else if (vpc->num_materials == 2)
	material_option = MATERIAL_TWO;
    else
	material_option = MATERIAL_MORE;
    if (vpc->color_channels == 1)
	color_option = COLOR_GRAY;
    else
	color_option = COLOR_RGB;

    /* render */
    if (vpc->affine_view) {
	/* choose a compositing function */
	composite_func = AffineProcTable[shadow_option][material_option]
					[color_option];
#ifdef INDEX_VOLUME
	if (have_index && shadow_option == SHADOWS_OFF &&
	    material_option == MATERIAL_ONE && color_option == COLOR_GRAY) {
	    composite_func = VPCompAI11B;
	}
#endif
	VPRenderAffine(vpc, USE_RLEVOLUME, composite_func);
    } else {
	/* XXX perspective rendering not available yet */
	return(VPSetError(vpc, VPERROR_BAD_OPTION));
    }

    return(VP_OK);
}
