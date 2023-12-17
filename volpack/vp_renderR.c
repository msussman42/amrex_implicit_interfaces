/*
 * vp_renderR.c
 *
 * Function to render raw (unclassified) volumes.
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
 * $Date: 2000-10-02 18:12:04 $
 * $Revision: 1.3 $
 */

#include "vp_global.h"

#define COMP_AR1PB_FUNC		VPCompAR1PB
extern void VPCompAR1PB(
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
			void *voxel_data,		/* voxel data for slice */
			int voxel_istride,		/* strides for voxel data */
			int voxel_jstride
			);

#define COMP_AR3PB_FUNC		VPCompAR3PB
extern void VPCompAR3PB(
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
			void *voxel_data,		/* voxel data for slice */
			int voxel_istride,		/* strides for voxel data */
			int voxel_jstride
			);

#ifdef COMP_AR11B
#define COMP_AR11B_FUNC		VPCompAR11B
extern void VPCompAR11B(
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
			void *voxel_data,		/* voxel data for slice */
			int voxel_istride,		/* strides for voxel data */
			int voxel_jstride
			);
#else
#define COMP_AR11B_FUNC		VPCompAR1NB
#endif

#ifdef COMP_AR31B
#define COMP_AR31B_FUNC		VPCompAR31B
extern void VPCompAR31B(
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
			void *voxel_data,		/* voxel data for slice */
			int voxel_istride,		/* strides for voxel data */
			int voxel_jstride
			);
#else
#define COMP_AR31B_FUNC		VPCompAR3NB
#endif

#ifdef COMP_AR12B
#define COMP_AR12B_FUNC		VPCompAR12B
extern void VPCompAR12B();
#else
#define COMP_AR12B_FUNC		VPCompAR1NB
#endif

#ifdef COMP_AR32B
#define COMP_AR32B_FUNC		VPCompAR32B
extern void VPCompAR32B(
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
			void *voxel_data,		/* voxel data for slice */
			int voxel_istride,		/* strides for voxel data */
			int voxel_jstride
			);
#else
#define COMP_AR32B_FUNC		VPCompAR3NB
#endif

#define COMP_AR1NB_FUNC		VPCompAR1NB
extern void VPCompAR1NB(
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
			void *voxel_data,		/* voxel data for slice */
			int voxel_istride,		/* strides for voxel data */
			int voxel_jstride
			);

#define COMP_AR3NB_FUNC		VPCompAR3NB
extern void VPCompAR3NB(
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
			void *voxel_data,		/* voxel data for slice */
			int voxel_istride,		/* strides for voxel data */
			int voxel_jstride
			);


#define COMP_AR1PS_FUNC		VPCompAR1PB

#define COMP_AR3PS_FUNC		VPCompAR3PB

#ifdef COMP_AR11S
#define COMP_AR11S_FUNC		VPCompAR11S
extern void VPCompAR11S();
#else
#define COMP_AR11S_FUNC		VPCompAR1NS
#endif

#ifdef COMP_AR31S
#define COMP_AR31S_FUNC		VPCompAR31S
extern void VPCompAR31S();
#else
#define COMP_AR31S_FUNC		VPCompAR3NS
#endif

#ifdef COMP_AR12S
#define COMP_AR12S_FUNC		VPCompAR12S
extern void VPCompAR12S();
#else
#define COMP_AR12S_FUNC		VPCompAR1NS
#endif

#ifdef COMP_AR32S
#define COMP_AR32S_FUNC		VPCompAR32S
extern void VPCompAR32S();
#else
#define COMP_AR32S_FUNC		VPCompAR3NS
#endif

#define COMP_AR1NS_FUNC		VPCompAR1NS
extern void VPCompAR1NS(
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
			void *voxel_data,		/* voxel data for slice */
			int voxel_istride,		/* strides for voxel data */
			int voxel_jstride,
			GrayIntPixel *shadow_buffer
			);

#define COMP_AR3NS_FUNC		VPCompAR3NS
extern void VPCompAR3NS(
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
			void *voxel_data,		/* voxel data for slice */
			int voxel_istride,		/* strides for voxel data */
			int voxel_jstride,
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

static void (*AffineProcTable[SHADOW_OPTS][MATERIAL_OPTS][COLOR_OPTS])() = {
    {
	{ (void(*)())COMP_AR1PB_FUNC, (void(*)())COMP_AR3PB_FUNC },
	{ (void(*)())COMP_AR11B_FUNC, (void(*)())COMP_AR31B_FUNC },
	{ (void(*)())COMP_AR12B_FUNC, (void(*)())COMP_AR32B_FUNC },
	{ (void(*)())COMP_AR1NB_FUNC, (void(*)())COMP_AR3NB_FUNC }
    },
    {
	{ (void(*)())COMP_AR1PS_FUNC, (void(*)())COMP_AR3PS_FUNC },
	{ (void(*)())COMP_AR11S_FUNC, (void(*)())COMP_AR31S_FUNC },
	{ (void(*)())COMP_AR12S_FUNC, (void(*)())COMP_AR32S_FUNC },
	{ (void(*)())COMP_AR1NS_FUNC, (void(*)())COMP_AR3NS_FUNC }
    }
};

/*
 * vpRenderRawVolume
 *
 * Render an uclassified volume using the shear-warp algorithm.
 */

vpResult
vpRenderRawVolume(vpc)
vpContext *vpc;
{
    int retcode;
    void (*composite_func)(vpContext*,...);
    int shadow_option, material_option, color_option;

    /* check for errors and initialize */
    if ((retcode = VPCheckRawVolume(vpc)) != VP_OK)
	return(retcode);
    if ((retcode = VPCheckClassifier(vpc)) != VP_OK)
	return(retcode);
    if ((retcode = VPCheckShader(vpc)) != VP_OK)
	return(retcode);
    if ((retcode = VPCheckImage(vpc)) != VP_OK)
	return(retcode);
    if (vpc->num_clsfy_params > 2)
	return(VPSetError(vpc, VPERROR_LIMIT_EXCEEDED));
    if ((retcode = VPFactorView(vpc)) != VP_OK)
	return(retcode);

    Debug((vpc, VPDEBUG_RENDER, "Algorithm: affine RAWvolume (%s)\n",
	   vpc->mm_octree == NULL ? "no octree" : "with octree"));

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
	composite_func = (void(*)(vpContext*,...))AffineProcTable[shadow_option][material_option]
					[color_option];
	VPRenderAffine(vpc, USE_RAWVOLUME, composite_func);
    } else {
	/* XXX perspective rendering not available yet */
	return(VPSetError(vpc, VPERROR_BAD_OPTION));
    }

    return(VP_OK);
}
