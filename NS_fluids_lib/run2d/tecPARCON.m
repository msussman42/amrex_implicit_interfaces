%% PARCON Data

hold off
clc
clear
close all

ok_to_continue=1

while (ok_to_continue==1)

 [filename,pathname]=uigetfile('*.tec')
 f_full=fullfile(pathname,filename)

 delimiterIn=' '; %separates based on spaces
 % line 1: TITLE= ....
 % line 2: VARIABLES= ....
 % line 3: ZONE F="POINT", I=      ??  J= 1   K= 1, SOLUTIONTIME= ????   
 %         STRANDID= ??
 Headerlines=3; %number of non-numerical rows
 A=[]; %empty array
 A=importdata(f_full,delimiterIn,Headerlines); 
 A_data_size=size(A.data); %finds dimensions of specific tab in array A
 
 ok_to_continue=0;
 
 varnames=A.textdata(2,1)
 zonedim=A.textdata(3,1) 
 time_strand=A.textdata(3,1) 

 varnames=strrep(varnames,'"',' ');
 varnames_char=char(varnames);
 varnames_char=strrep(varnames_char,'VARIABLES = ',' ');
 menu_list=strsplit(varnames_char,',')
 menu_list_size=size(menu_list)
 
 zonedim_char=char(zonedim);
 zonedim_char=strrep(zonedim_char,'f="point"',' ');
 zonedim_char=strrep(zonedim_char,'F="POINT"',' ');
 zonedim_char=strrep(zonedim_char,'i=',' ');
 zonedim_char=strrep(zonedim_char,'I=',' ');
 zonedim_char=strrep(zonedim_char,'J=',' ');
 zonedim_char=strrep(zonedim_char,'K=',' ');
 zonedim_char=strrep(zonedim_char,'ZONE',' ');
 zonedim_char=strrep(zonedim_char,'SOLUTIONTIME=',' ');
 zonedim_char=strrep(zonedim_char,'STRANDID=',' ');
 zonedim_char=strrep(zonedim_char,',',' ');
 zonedata=str2num(zonedim_char)

 number_particles=round(zonedata(1,1))
 solutiontime=zonedata(1,4)
 strand_id_int=round(zonedata(1,5))

 answer=0
 while (answer==0)

  answer=input('enter menu list number from 1 to menu_list_size : ')
  if ((answer>=1)&&(answer<=menu_list_size(2)))
   % do nothing
  else
   disp('answer out of range, try again')
   answer=0
  end  

  arow=1
  for np=1:number_particles
   xpart=A.data(arow,1)
   ypart=A.data(arow,2)
   arow=arow+1;
   % plot the particle here ...
  end
  
 end
end
 
 
