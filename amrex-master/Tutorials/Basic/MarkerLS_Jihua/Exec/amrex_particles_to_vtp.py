import numpy as np
from sys import argv
from write_pview_file import write_paraview_file_particles

def read_amrex_ascii_particle_file(fname):

    infile=open(fname,'r')
    line=infile.readline()#number of particles
    nparticles=int(line.split()[0])

    line_ignore1=infile.readline()#num SoA  real
    line_ignore2=infile.readline()#num SoA  int
    line_ignore3=infile.readline()#num AoS real
    line_ignore4=infile.readline()#num AoS  int

    particle_pos=np.zeros((nparticles,3))
    particle_scalar=np.zeros((1,nparticles))

    for i in range(nparticles):
        spltline=infile.readline().split()
        particle_pos[i][0]=float(spltline[0])
        particle_pos[i][1]=float(spltline[1])
        particle_pos[i][2]=float(spltline[2])
        particle_scalar[0][i]=float(spltline[3]) 

    return(particle_pos,particle_scalar,nparticles)


if __name__ == "__main__":

    filenum_min=int(argv[1])
    filenum_max=int(argv[2])
    file_prefix=argv[3]
    
    skip = int(argv[4])

    #nfiles=filenum_max-filenum_min+1
    ncdata=np.array([])
    for i in range(filenum_min, filenum_max+skip, skip):
        filename=file_prefix+"%5.5d"%(filenum_min+i)
        (ppos,ncdata,npart)=read_amrex_ascii_particle_file(filename)
        write_paraview_file_particles(filename+".vtp",ppos,ncdata)


#particle data:
#x,y,z,id,cpu,mask,i_x,i_y,i_z