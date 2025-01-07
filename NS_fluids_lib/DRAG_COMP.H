#ifndef DRAG_COMP_H_
#define DRAG_COMP_H_

#if (0)
This enumerated type indexes into DRAG_MF which has
ngrow_make_distance=3 ghost cells.  Information on the
stress applied to a given material "im_test" is stored.
#endif

#define ENUM_T11 (0)
#define ENUM_T12 (1)
#define ENUM_T22 (2)
#define ENUM_T33 (3)
#define ENUM_T13 (4)
#define ENUM_T23 (5)
#define ENUM_NUM_TENSOR_TYPE (2*AMREX_SPACEDIM)
#define NUM_CELL_ELASTIC (num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE)

#define DRAG_TYPE_T11 (0)
#define DRAG_TYPE_T12 (1)
#define DRAG_TYPE_T22 (2)
#define DRAG_TYPE_T33 (3)
#define DRAG_TYPE_T13 (4)
#define DRAG_TYPE_T23 (5)
#define DRAG_TYPE_UVEC (6)
#define DRAG_TYPE_VVEC (7)
#define DRAG_TYPE_WVEC (8)
#define DRAG_TYPE_SCALAR (9)
#define DRAG_TYPE_FLAG (10)
#define DRAG_TYPE_NEXT (11)

#define DRAG_TYPE_IQ_UFORCE (0)
#define DRAG_TYPE_IQ_VFORCE (1)
#define DRAG_TYPE_IQ_WFORCE (2)
#define DRAG_TYPE_IQ_SCALAR (3)
#define DRAG_TYPE_IQ_NEXT (4)

#define DRAGCOMP_IQ_BODYFORCE (0)
#define DRAGCOMP_IQ_FORCE (DRAGCOMP_IQ_BODYFORCE+3*num_materials)
#define DRAGCOMP_IQ_PFORCE (DRAGCOMP_IQ_FORCE+3*num_materials)
#define DRAGCOMP_IQ_VISCOUSFORCE (DRAGCOMP_IQ_PFORCE+3*num_materials)
#define DRAGCOMP_IQ_VISCOUS0FORCE (DRAGCOMP_IQ_VISCOUSFORCE+3*num_materials)
#define DRAGCOMP_IQ_VISCOFORCE (DRAGCOMP_IQ_VISCOUS0FORCE+3*num_materials)

#define DRAGCOMP_IQ_BODYTORQUE (DRAGCOMP_IQ_VISCOFORCE+3*num_materials)
#define DRAGCOMP_IQ_TORQUE (DRAGCOMP_IQ_BODYTORQUE+3*num_materials)
#define DRAGCOMP_IQ_PTORQUE (DRAGCOMP_IQ_TORQUE+3*num_materials)
#define DRAGCOMP_IQ_VISCOUSTORQUE (DRAGCOMP_IQ_PTORQUE+3*num_materials)
#define DRAGCOMP_IQ_VISCOUS0TORQUE (DRAGCOMP_IQ_VISCOUSTORQUE+3*num_materials)
#define DRAGCOMP_IQ_VISCOTORQUE (DRAGCOMP_IQ_VISCOUS0TORQUE+3*num_materials)

#define DRAGCOMP_IQ_COM (DRAGCOMP_IQ_VISCOTORQUE+3*num_materials)
#define DRAGCOMP_IQ_MOMINERTIA (DRAGCOMP_IQ_COM+3*num_materials)
#define DRAGCOMP_IQ_MASS (DRAGCOMP_IQ_MOMINERTIA+3*num_materials)
#define DRAGCOMP_IQ_PERIM (DRAGCOMP_IQ_MASS+num_materials)
#define N_DRAG_IQ (DRAGCOMP_IQ_PERIM+num_materials)

#define DRAGCOMP_STRESS (0)
#define DRAGCOMP_PSTRESS (DRAGCOMP_STRESS+6*num_materials)
#define DRAGCOMP_VISCOUSSTRESS (DRAGCOMP_PSTRESS+6*num_materials)
#define DRAGCOMP_VISCOUS0STRESS (DRAGCOMP_VISCOUSSTRESS+6*num_materials)
#define DRAGCOMP_VISCOSTRESS (DRAGCOMP_VISCOUS0STRESS+6*num_materials)

#define DRAGCOMP_TORQUE_ARM (DRAGCOMP_VISCOSTRESS+6*num_materials)
#if (0)
DRAGCOMP_FLAG=0 => drag uninitialized
DRAGCOMP_FLAG=1 => drag initialized
DRAGCOMP_FLAG=2 => drag extended
#endif
#define DRAGCOMP_FLAG (DRAGCOMP_TORQUE_ARM+3*num_materials)
#define DRAGCOMP_VISCOUS0STRESSMAG (DRAGCOMP_FLAG+num_materials)
#define DRAGCOMP_VISCOUS0STRESSTAN (DRAGCOMP_VISCOUS0STRESSMAG+num_materials)

#define N_DRAG (DRAGCOMP_VISCOUS0STRESSTAN+3*num_materials)

#endif
