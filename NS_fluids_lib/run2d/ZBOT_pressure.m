probeF=load('pressure_T1probe_uniform_64x256');
probeC=load('pressure_T1probe_32x128');
avgF=load('pressure_average_uniform_64x256');
avgC=load('pressure_average_32x128');
plot(probeC(:,1),probeC(:,2))
hold on
plot(avgC(:,1),avgC(:,2))
plot(probeF(:,1),probeF(:,2))
plot(avgF(:,1),avgF(:,2))
avgF_AMR=load('pressure_average_64x256');
probeF_AMR=load('pressure_T1probe_64x256');
plot(probeF_AMR(:,1),probeF_AMR(:,2))
plot(avgF_AMR(:,1),avgF_AMR(:,2))
legend("T1 probe pressure 32x128","average ullage pressure 32x128","T1 probe pressure (uniform) 64x256","average ullage pressure (uniform) 64x256","T1 probe pressure (AMR) 64x256","average ullage pressure (AMR) 64x256")
