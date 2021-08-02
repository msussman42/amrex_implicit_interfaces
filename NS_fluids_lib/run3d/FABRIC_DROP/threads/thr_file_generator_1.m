clear
close all
clc

sub_dir="./";
write_files=0;

if(write_files)
    system(strcat("rm -f ",sub_dir,"thread_*"));
end


% Periodic domain in X and Y direction 
domain_lo=[-1.5 -1.5 -1.4];
domain_hi=[ 1.5  1.5  1.6];
thread_z=0;

drop_location=[0.0 0.0 1.05];
drop_radius=0.5;

% numbers should be even for the sake of symmetry in the periodic domain
n_1 = 4; % number of threads prallel to Y axis (flat)
n_2 = 4; % number of threads prallel to X axis (wavy)
% number of thread files will be (n_1 + 2) + (n_2 + 2)

dx_th=0.01; % separation of sample points in X direction
grid_cell=[192 192 192];  % fluid grid cells of the finest level

r_1 = 0.14; % thread radius (flat threads)
r_2 = r_1; % thread radius (wavy threads)
a_wavy=0.4; % amplitude of the wavy threads


figure('units','normalized','outerposition',[-1 0 1 1])
hold on

%% flat threads

d_1=(domain_hi(1)-domain_lo(1))/n_1;
X_1=domain_lo(1)-d_1/2 : d_1 : domain_hi(1)+d_1/2;
if(write_files)
    for i=1:n_1+2
        filename=sprintf("%sthread_%03d.txt",sub_dir,i);
        fid=fopen(filename,'w');
        fprintf(fid,"%d %16.15f\n",2,r_1);
        fprintf(fid,"%16.15f %16.15f %16.15f\n",...
            X_1(i),domain_lo(2)-(domain_hi(2)-domain_lo(2))/2,thread_z);
        fprintf(fid,"%16.15f %16.15f %16.15f\n",...
            X_1(i),domain_hi(2)+(domain_hi(2)-domain_lo(2))/2,thread_z);
        fclose(fid);
    end
end
%% wavy threads
omega = pi/d_1;

n_dx= ceil((domain_hi(1)-domain_lo(1))/dx_th);
X_2= linspace(domain_lo(1)-(domain_hi(1)-domain_lo(1))/2,...
              domain_hi(1)+(domain_hi(1)-domain_lo(1))/2,n_dx+1);

          z_wavy_p= a_wavy*sin(omega*(X_2-domain_lo(1)))+thread_z;
z_wavy_n=-a_wavy*sin(omega*(X_2-domain_lo(1)))+thread_z;

%% Plot wavy threads and r_2 level set 
plot(X_2,z_wavy_p,'k');
wavy_p_p=zeros(2,length(X_2)-1);
wavy_p_n=zeros(2,length(X_2)-1);
for i=1+1:length(X_2)-1
    p_1=[X_2(i-1); z_wavy_p(i-1)];
    p_2=[X_2(i)  ; z_wavy_p(i)  ];
    p_3=[X_2(i+1); z_wavy_p(i+1)];
    vec=p_3-p_1;
    vec=vec/norm(vec,2);
    vec=[0 -1;1 0]*vec;
    wavy_p_p(:,i)=p_2+r_2*vec;
    wavy_p_n(:,i)=p_2-r_2*vec;
end
plot(wavy_p_p(1,:),wavy_p_p(2,:),'r.');
plot(wavy_p_n(1,:),wavy_p_n(2,:),'r.');

plot(X_2,z_wavy_n,'k');
wavy_n_p=zeros(2,length(X_2)-1);
wavy_n_n=zeros(2,length(X_2)-1);
for i=1+1:length(X_2)-1
    p_1=[X_2(i-1); z_wavy_n(i-1)];
    p_2=[X_2(i)  ; z_wavy_n(i)  ];
    p_3=[X_2(i+1); z_wavy_n(i+1)];
    vec=p_3-p_1;
    vec=vec/norm(vec,2);
    vec=[0 -1;1 0]*vec;
    wavy_n_p(:,i)=p_2+r_2*vec;
    wavy_n_n(:,i)=p_2-r_2*vec;
end
plot(wavy_n_p(1,:),wavy_n_p(2,:),'r.');
plot(wavy_n_n(1,:),wavy_n_n(2,:),'r.');

%% plot thread circular cross sections
theta=linspace(0,2*pi,100);
circle_X=r_1*cos(theta);
circle_Z=r_1*sin(theta);

for i = 1:n_1+2
    plot(X_1(i)+circle_X,thread_z+circle_Z,'b');
end

%% plot fluid grid
xgrid=linspace(domain_lo(1),domain_hi(1),grid_cell(1)+1);
for i=1:length(xgrid)
    plot([xgrid(i) xgrid(i)],[domain_lo(3) domain_hi(3)],'g');
end
zgrid=linspace(domain_lo(3),domain_hi(3),grid_cell(3)+1);
for i=1:length(zgrid)
    plot([domain_lo(1) domain_hi(1)],[zgrid(i) zgrid(i)],'g');
end

%% write wavy tread files
d_2=(domain_hi(2)-domain_lo(2))/n_2;
Y_2=domain_lo(2)-d_2/2 : d_2 : domain_hi(2)+d_2/2;
if(write_files)
    for i=n_1+3:n_1+3+n_2+1
        if(mod(i,2)==0)
            z_wavy=z_wavy_p;
        else
            z_wavy=z_wavy_n;
        end
        filename=sprintf("%sthread_%03d.txt",sub_dir, i);
        fid=fopen(filename,'w');
        fprintf(fid,"%d %16.15f\n",length(X_2),r_2);
        for j=1:length(X_2)
            fprintf(fid,"%16.15f %16.15f %16.15f\n",...
                X_2(j),Y_2(i-(n_1+2)),z_wavy(j));
        end
        fclose(fid);
    end
end

%% plot drop
theta=linspace(0,2*pi,100);
circle_X=drop_radius*cos(theta);
circle_Z=drop_radius*sin(theta);

for i = 1:n_1+2
    plot(drop_location(1)+circle_X,drop_location(3)+circle_Z,'b');
end

xlabel('X');
ylabel('Z');
axis equal



