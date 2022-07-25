#ifndef INDEX_TYPE_MACROS_H_
#define INDEX_TYPE_MACROS_H_

#define TheUMACType IndexType(AMREX_D_DECL(IndexType::NODE, IndexType::CELL, IndexType::CELL))
#define TheXUMACType IndexType(AMREX_D_DECL(IndexType::CELL, IndexType::CELL, IndexType::CELL))
#define TheYUMACType IndexType(AMREX_D_DECL(IndexType::NODE, IndexType::NODE, IndexType::CELL))
#define TheZUMACType IndexType(AMREX_D_DECL(IndexType::NODE, IndexType::CELL, IndexType::NODE))
#define TheVMACType IndexType(AMREX_D_DECL(IndexType::CELL, IndexType::NODE, IndexType::CELL))
#define TheXVMACType IndexType(AMREX_D_DECL(IndexType::NODE, IndexType::NODE, IndexType::CELL))
#define TheYVMACType IndexType(AMREX_D_DECL(IndexType::CELL, IndexType::CELL, IndexType::CELL))
#define TheZVMACType IndexType(AMREX_D_DECL(IndexType::CELL, IndexType::NODE, IndexType::NODE))
#define TheWMACType IndexType(AMREX_D_DECL(IndexType::CELL, IndexType::CELL, IndexType::NODE))
#define TheXWMACType IndexType(AMREX_D_DECL(IndexType::NODE, IndexType::CELL, IndexType::NODE))
#define TheYWMACType IndexType(AMREX_D_DECL(IndexType::CELL, IndexType::NODE, IndexType::NODE))
#define TheZWMACType IndexType(AMREX_D_DECL(IndexType::CELL, IndexType::CELL, IndexType::CELL))

#endif
