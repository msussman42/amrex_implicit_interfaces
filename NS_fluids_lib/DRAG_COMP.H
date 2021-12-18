#ifndef DRAG_COMP_H_
#define DRAG_COMP_H_

#if (0)
This enumerated type indexes into DRAG_MF which has
ngrow_make_distance=3 ghost cells and gives the local force/torque
applied to a given material "im_test"
#endif

#define DRAG_TYPE_UFORCE (0)
#define DRAG_TYPE_VFORCE (1)
#define DRAG_TYPE_WFORCE (2)
#define DRAG_TYPE_SCALAR (3)
#define DRAG_TYPE_FLAG (4)
#define DRAG_TYPE_NEXT (5)

#define DRAGCOMP_BODYFORCE (0)
#define DRAGCOMP_FORCE (DRAGCOMP_BODYFORCE+3*num_materials)
#define DRAGCOMP_PFORCE (DRAGCOMP_FORCE+3*num_materials)
#define DRAGCOMP_VISCOUSFORCE (DRAGCOMP_PFORCE+3*num_materials)
#define DRAGCOMP_VISCOUS0FORCE (DRAGCOMP_VISCOUSFORCE+3*num_materials)
#define DRAGCOMP_VISCOFORCE (DRAGCOMP_VISCOUS0FORCE+3*num_materials)

#define DRAGCOMP_BODYTORQUE (DRAGCOMP_VISCOFORCE+3*num_materials)
#define DRAGCOMP_TORQUE (DRAGCOMP_BODYTORQUE+3*num_materials)
#define DRAGCOMP_PTORQUE (DRAGCOMP_TORQUE+3*num_materials)
#define DRAGCOMP_VISCOUSTORQUE (DRAGCOMP_PTORQUE+3*num_materials)
#define DRAGCOMP_VISCOUS0TORQUE (DRAGCOMP_VISCOUSTORQUE+3*num_materials)
#define DRAGCOMP_VISCOTORQUE (DRAGCOMP_VISCOUS0TORQUE+3*num_materials)

#define DRAGCOMP_COM (DRAGCOMP_VISCOTORQUE+3*num_materials)
#define DRAGCOMP_MOMINERTIA (DRAGCOMP_COM+3*num_materials)
#define DRAGCOMP_MASS (DRAGCOMP_MOMINERTIA+3*num_materials)
#define DRAGCOMP_PERIM_VECTOR (DRAGCOMP_MASS+num_materials)
#define DRAGCOMP_PERIM (DRAGCOMP_PERIM_VECTOR+3*num_materials)
#if (0)
DRAGCOMP_FLAG=0 => drag uninitialized
DRAGCOMP_FLAG=1 => drag initialized
DRAGCOMP_FLAG=2 => drag extended
#endif
#define DRAGCOMP_FLAG (DRAGCOMP_PERIM+num_materials)
#define N_DRAG (DRAGCOMP_FLAG+num_materials)

#endif
