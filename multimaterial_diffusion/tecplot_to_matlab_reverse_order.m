hold off
clc
clear all
close all

box_old=zeros(1,3);
box_new=zeros(1,3);
num_plotted=0;
A3D_new=zeros(2,2,2);
A3D_old=zeros(2,2,2);
xaxis_new=zeros(1,2);
xaxis_old=zeros(1,2);
yaxis_new=zeros(1,2);
yaxis_old=zeros(1,2);

ok_to_continue=1

while (ok_to_continue==1)

 num_plotted
 [filename,pathname]=uigetfile('*.tec');
 f_full=fullfile(pathname,filename)

 delimiterIn=' ' %separates based on spaces
 % line 1: VARIABLES= ....
 % line 2: zone i=      ??  j=    ??   k=   ??   f=point
 % line 3: SOLUTIONTIME= ????   STRANDID=   ??
 Headerlines=3 %number of non-numerical rows
 A=[]; %empty array
 A=importdata(f_full,delimiterIn,Headerlines); 
 A_data_size=size(A.data) %finds dimensions of specific tab in array A

 varnames=A.textdata(1,1)
 zonedim=A.textdata(2,1)
 time_strand=A.textdata(3,1)

 varnames=strrep(varnames,'VARIABLES=',' ')
 varnames=strrep(varnames,'"',' ')
 varnames_char=char(varnames)
 menu_list=strsplit(varnames_char,',')
 menu_list_size=size(menu_list)

 answer=0;
 while (answer==0)

  answer=input('enter menu list number from 1 to menu_list_size : ')
  if ((answer>=1)&&(answer<=menu_list_size(2)))
   % do nothing
  else
   disp('answer out of range, try again')
   answer=0
  end       

 end

 v=input('Enter contour values in the form [x,x,...,x]:');

 zonedim=strrep(zonedim,'f=point',' ')
 zonedim=strrep(zonedim,'zone i=',' ')
 zonedim=strrep(zonedim,'j=',' ')
 zonedim_char=char(zonedim)
 ij_dimen=str2num(zonedim_char)

 time_strand=strrep(time_strand,'SOLUTIONTIME=',' ')
 time_strand=strrep(time_strand,'STRANDID=',' ')
 time_strand_char=char(time_strand)
 time_strand_vals=str2num(time_strand_char)
 time_real=time_strand_vals(1,1)
 strand_id_int=round(time_strand_vals(1,2))

 imax=ij_dimen(1,1)
 jmax=ij_dimen(1,2)
 imax*jmax
 A_data_size()
 if (ij_dimen(1,1)*ij_dimen(1,2)==A_data_size(1))
  % do nothing
 else
  exit()
 end

 ival=1
 jval=1
 xaxis=zeros(1,imax);
 yaxis=zeros(1,jmax);
 A3D=zeros(jmax,imax,A_data_size(2));

 for arow=1:A_data_size(1)
  if (ival==1)
   yaxis(1,jval)=A.data(arow,2);
  end
  if (jval==1)
   xaxis(1,ival)=A.data(arow,1);
  end
  for idata=1:A_data_size(2)
   A3D(jval,ival,idata)=A.data(arow,idata);
  end
  jval=jval+1;
  if (jval>jmax)
   ival=ival+1;
   jval=1;
  end
 end

 if (mod(num_plotted,4)==0)
  LineSpec='red';
 elseif (mod(num_plotted,4)==1)
  LineSpec='green';
 elseif (mod(num_plotted,4)==2)
  LineSpec='blue';
 elseif (mod(num_plotted,4)==3)
  LineSpec='black';
 end
 [M,c]=contour(yaxis,xaxis,A3D(:,:,answer),v,LineSpec);
 c.LineWidth=2;
 hold on

 num_plotted=num_plotted+1
 if (num_plotted>1)
  for dir=1:3
   box_old(1,dir)=box_new(1,dir);
  end
  A3D_old=zeros(box_old(1,2),box_old(1,1),box_old(1,3));
  xaxis_old=zeros(1,box_old(1,1)); 
  yaxis_old=zeros(1,box_old(1,2));
  for ival=1:box_old(1,1)
  for jval=1:box_old(1,2)
  for idata=1:box_old(1,3)
   A3D_old(jval,ival,idata)=A3D_new(jval,ival,idata);
  end
  end
  end
  for ival=1:box_old(1,1)
   xaxis_old(1,ival)=xaxis_new(1,ival);
  end
  for jval=1:box_old(1,2)
   yaxis_old(1,jval)=yaxis_new(1,jval);
  end

  box_new(1,1)=imax;
  box_new(1,2)=jmax;
  box_new(1,3)=A_data_size(2);

  A3D_new=zeros(box_new(1,2),box_new(1,1),box_new(1,3));
  xaxis_new=zeros(1,box_new(1,1));
  yaxis_new=zeros(1,box_new(1,2));
  for ival=1:box_new(1,1)
  for jval=1:box_new(1,2)
  for idata=1:box_new(1,3)
   A3D_new(jval,ival,idata)=A3D(jval,ival,idata);
  end
  end
  end
  for ival=1:box_new(1,1)
   xaxis_new(1,ival)=xaxis(1,ival);
  end
  for jval=1:box_new(1,2)
   yaxis_new(1,jval)=yaxis(1,jval);
  end
  error_counter=0;
  L1norm=0.0;
  dx=xaxis_old(1,2)-xaxis_old(1,1)

  lq_array=zeros(box_old(1,2),box_old(1,1));

  for ival=1:box_old(1,1)
  for jval=1:box_old(1,2)
   data_old=A3D_old(jval,ival,answer);
   xq=xaxis_old(1,ival);
   yq=yaxis_old(1,jval);
   % interpolate from the new (fine) grid to the old (coarse) grid.
   data_new=interp2(yaxis_new,xaxis_new,A3D_new(:,:,answer),xq,yq);
   lq_array(jval,ival)=data_new;
   if ((abs(data_old)<=dx)||(abs(data_new)<=dx))
    error_counter=error_counter+1;
    L1norm=L1norm+abs(data_old-data_new);
   end
  end
  end
%  [M,c]=contour(yaxis,xaxis,lq_array,v);
%  c.LineWidth=3;

  if (error_counter>0) 
   L1norm=L1norm/error_counter;
  end
  error_counter
  L1norm
 elseif (num_plotted==1)
  box_new(1,1)=imax;
  box_new(1,2)=jmax;
  box_new(1,3)=A_data_size(2);
  A3D_new=zeros(box_new(1,2),box_new(1,1),box_new(1,3));
  xaxis_new=zeros(1,box_new(1,1));
  yaxis_new=zeros(1,box_new(1,2));
  for ival=1:box_new(1,1)
  for jval=1:box_new(1,2)
  for idata=1:box_new(1,3)
   A3D_new(jval,ival,idata)=A3D(jval,ival,idata);
  end
  end
  end
  for ival=1:box_new(1,1)
   xaxis_new(1,ival)=xaxis(1,ival);
  end
  for jval=1:box_new(1,2)
   yaxis_new(1,jval)=yaxis(1,jval);
  end
 end

 ok_to_continue=input('enter "1" to compare another data set, "0" to end')
end
