%Copyright (c) <2017>, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory
%Written by Tanya Kostova Vassilevska, kostova@llnl.gov, tan.v.kos@gmail.com
%LLNL Release number LLNL-CODE-735916
%All rights reserved.

%This file is part of <ERROM>, the main routine. For details, see the comments. 

%Licensed under the Apache License, Version 2.0 (the “Licensee”); you may not use this file except in compliance with the License.  You may obtain a copy of the License at:  http://www.apache.org/licenses/LICENSE-2.0

%Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the license.


%This is a code to calculate the approximation error from Proper Orthogonal Decomposition (POD) Model Reduction methods and to compare and analyze the errors for two variants of the method. 
%The code can be used to develop further variants of the methods 

%%%%%%%%%%%%%%%%%%%
% The code solves an ODE system, defined in the subroutine ODEdef.m and keeps the solution in array q
% The solutions is calculated using Matlab’s solver ode15s
% Then the code creates two TYPES of POD Reduced Order Models (ROMs):
% M1. Using snapshots of the solution at l equidistant time points kept in array y. The lifted ROM solution is kept in array b0. Two different ROMs are created - of dimension l and 2l. 
% M2. Using snapshots of the solution and the derivative at the same time points (equally spaced). The derivative snapshots are kept in array z. Two different ROMs are created - of dimension l and 2l. 
%
%The L2 norm of the point-wise in time errors between the solution and the ROM solutions is calculated and its decimal log is plotted. 
%
%Other calculated quantities include the singular values and the eigenvalues of a projected Jacobian. These are saved in files eefile.eps and Eigenvalues.txt
%
%The code is easily adaptable to calculate ROMs of different - user defined dimensions, to compare the errors for changed parameters of the system, etc.
%
%The output is plots of the error from the 3 methods, as well as a text
%file (Results.txt) outputting numeric outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear all;
colormap;


%%%%%%%%%%%%%% Define parameters of the FOM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Here we define parameters for the FitzHugh - Nagumo model; 
%%%%%%%%%%%%%  if another FOM is solved, user has to define the model-specific parameters as global  and assign values to them here

global delta1 delta2 a lambda gamma r0 r1 xr xl P L mu dx2 dx offset1 offset2 offset3

s_error=1.0e-15;% the cutoff for the singular values 


xl=0.0; 
xr=100.0; %left and right ends of x-interval

w0=0; %boundary  conditions
w1=0;

offset1=0.0; 
offset2=0.0; 
offset3=0.0;%offset parameters are meant to control the stability of the linear model


delta10=2; 
delta20=0.1;
a0=0.1; 
lambda0=2; 
gamma0=0.1; 
r00=10;
r10=30;
mu0=50;

%%%%%%%%%%%% This version allows for changing the parameters of the model to analyze the error of approximation of the ROMs for different parameters
%%%%%%%%%%%% First setup the parameters already defined to be the default
delta1=delta10; 
delta2=delta20; 
a=a0; 
lambda=lambda0; 
gamma=gamma0; 
r0=r00;
r1=r10;
mu=mu0;

%%%%%%%%%%% Then define the quantities to change the parameters; here they have been set to 0 but can be defined by the user
dsum=0;
dlambda0=0; 
ddelta10=0.; 
ddelta20=0.; 
da0=0.; 
dgamma0=0.; 
dr00=0.;
dr10=0.;
dmu0=0.;
dp0=0.;
dp1=0.;
df0=0.;
df1=0.;

dsum=dsum+dlambda0^2;
dsum=dsum+ddelta10^2;
dsum=dsum+ddelta20^2;
dsum=dsum+da0^2;
dsum=dsum+dgamma0^2;
dsum=dsum+dr00^2;
dsum=dsum+dr10^2;
dsum=dsum+dmu0^2;
dsum=dsum+dp0^2;
dsum=dsum+dp1^2;
dsum=dsum+df0^2;
dsum=dsum+df1^2;



fid=fopen('Results.txt', 'w'); % Dump all results in this file
feig=fopen('Eigenvalues.txt','w'); % the eigenvalues of the Jacobian of the FOM linearization

n=802; %n is the number of ODEs resulting from method of lines discretization of the FH-N equations; This is the dimension of the full order model

%%%%%%%%%% The following definitions (L, dx, dx2) are specific for the FH-N model and just define dx and dx2; they are not needed for other models

L= n/2-1; %CHOOSE n so that L/4 is integer

% L=n/2-1 determines step in space (xfinal-x0)/L is the space step; 
% n=2(L+1) determines the dimension of the ODE resulting from method of lines
%since the PDE has two variables - V and w - the first half of ODEs
%correspond to the discretization of V and the second to discretization of
% w; i.e. the first L+1 coordinates belong to the V variable, 
% the second L+1 belong to w
% 

dx=((xr-xl)/L);
dx2=dx^2;

%dx=1.0;
%dx2=1.0;


lambda0=2; 

%%%%%%%%%%% Define simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%
t0=0.0;%initial time
tfinal=2; %snapshots are collected from the interval [t0, tfinal]
tend=2; %solutions are compared on the interval [t0,tend]- it can be of different length compared to [t0,tfinal], smaller or larger

N=300;  %number of snapshots
N0=N;

m=100; % m determines step in time for the output, used for plotting; calculations otherwise are carried out with step that Matlab provides

X=5;% number of different parameter sets used to produce examples
Step=2;% used to produce a change in the parameter that is varied when defining a new set 

P=eye(n); %Initialize the Projection matrix P as the Unit matrix


%%%%%%%%%%%%% Define the initial conditions %%%%%%%%%%%%%%%%%%%%
for i=1:n
    yinit(i)=0.;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%% Define files for graphic output %%%%%%%%%%%%%%%%%%%
%%%%%%%% NOTE: All plots and labels  are FH-N specific %%%%

ffnames={'logErrorSolDerComp1first_1', 'logErrorSolDerComp1first_2', 'logErrorSolDerComp1first_3'}; % plots of  the decimal logarithms of the L2 distances between b0 and q and b1 and q; i.e. errors from M1 in space with dim = 2*k0=2*rank(y) and M2  in space with dim k1=rank (y,z) FOR A GIVEN TIME STEP; in this version only one plot of this type is created (because k=1 in the first loop), there is a possibility for 3 plots; this is valid for all other
fqnames={'logErrorSolDerComp1second_1', 'logErrorSolDerComp1second_2', 'logErrorSolDerComp1second_3'}; % plots of  the decimal logarithms of the L2 distances between b0 and q and b1 and q; i.e. errors from M1 in space with dim = 2*k0=2*rank(y) and M2  in space with dim k1=rank (y,z) FOR 10 times larger  TIME STEP
ddnames={'logErrorSolDerComp1third_1', 'dlogErrorSolDerComp1third_2', 'logErrorSolDerComp1third_3'}; % plots of  the decimal logarithms of the L2 distances between b0 and q and b1 and q; i.e. errors from M1 in space with dim = 2*k0=2*rank(y) and M2  in space with dim k1=rank (y,z) FOR 100 times larger  TIME STEP


ffnames2={'logErrorSolDerComp1first_12', 'logErrorSolDerComp1first_22', 'logErrorSolDerComp1first_32'};%plots of  the decimal logarithms of the L2 distances between b0 and q and b2 and q; i.e. errors from M1 in space with dim = 2*k0=2*rank(y) and M2  in space with dim k2=2*k0=2*rank(y) FOR A GIVEN TIME STEP
fqnames2={'logErrorSolDerComp1second_12', 'logErrorSolDerComp1second_22', 'logErrorSolDerComp1second_32'};%plots of  the decimal logarithms of the L2 distances between b0 and q and b2 and q; i.e. errors from M1 in space with dim = 2*k0=2*rank(y) and M2  in space with dim k2=2*k0=2*rank(y) FOR 10 times larger  TIME STEP
ddnames2={'logErrorSolDerComp1third_12', 'dlogErrorSolDerComp1third_22', 'logErrorSolDerComp1third_e32'};%plots of  the decimal logarithms of the L2 distances between b0 and q and b2 and q; i.e. errors from M1 in space with dim = 2*k0=2*rank(y) and M2  in space with dim k2=2*k0=2*rank(y) FOR 100 times larger  TIME STEP

fnames={'logErrorSolDerComp3_1', 'logErrorSolDerComp3_2', 'logErrorSolDerComp3_3'}; %the plots from ffnames, fqnames and ddnames on the same plot
fnames2={'logErrorSolDerComp3_12', 'logErrorSolDerComp3_22', 'logErrorSolDerComp3_32'};%the plots from ffnames2, fqnames2 and ddnames2 on the same plot

Vplots={'Vplot1', 'Vplot2', 'Vplot3', 'Vplot4', 'Vplot5', 'Vplot6', 'Vplot7', 'Vplot8', 'Vplot9', 'Vplot10','Vplot11', 'Vplot12', 'Vplot13', 'Vplot14', 'Vplot15', 'Vplot16', 'Vplot17', 'Vplot18', 'Vplot19', 'Vplot20'};% these are waterfall (3D) plots of the spatiotemporal dynamics of the solution V of the FGH-N FOM (i.e. these are model specific)
Wplots={'Wplot1','Wplot2','Wplot3','Wplot4','Wplot5','Wplot6','Wplot7','Wplot8','Wplot9','Wplot10','Wplot11','Wplot12','Wplot13','Wplot14','Wplot15','Wplot16','Wplot17','Wplot18','Wplot19','Wplot20' };% these are waterfall (3D) plots of the spatiotemporal dynamics of the solution W of the FH-N FOM

ggnames={'logErrorSolComp3_1', 'logErrorSolComp3_2', 'logErrorSolComp3_3'}; % plots of the dec log of the error from M1 (ROM dim=2*k0) for the three time steps
gvnames={'logVErrorSolComp3_1', 'logVErrorSolComp3_2', 'logVErrorSolComp3_3'};% plots of the dec log of the error only for V from M1 (ROM dim=2*k0) for the three time steps
gwnames={'logWErrorSolComp3_1', 'logWErrorSolComp3_2', 'logWErrorSolComp3_3'};% plots of the dec log of the error only for W from M1 (ROM dim=2*k0) for the three time steps

fdnames={'logErrorDerComp3_1', 'logErrorDerComp3_2', 'logErrorDerComp3_3'};% plots of the dec log of the error from M1 (ROM dim=k1) for the three time steps
fvnames={'logVErrorDerComp3_1', 'logVErrorDerComp3_2', 'logVErrorDerComp3_3'};% plots of the dec log of the error only for V from M1 (ROM dim=k1) for the three time steps
fwnames={'logWErrorDerComp3_1', 'logWErrorDerComp3_2', 'logWErrorDerComp3_3'};% plots of the dec log of the error only for W from M1 (ROM dim=k1) for the three time steps

fdnames2={'logErrorDerHalfDimComp3_1', 'logErrorDerHalfDimComp3_2', 'logErrorDerHalfDimComp3_3'};% plots of the dec log of the error from M2 (ROM dim=2*k0) for the three time steps
fvnames2={'logVErrorDerHalfDimComp3_1', 'logVErrorDerHalfDimComp3_2', 'logVErrorDerHalfDimComp3_3'};% plots of the dec log of the error only for V from M2 (ROM dim=2*k0) for the three time steps
fwnames2={'logWErrorDerHalfDimComp3_1', 'logWErrorDerHalfDimComp3_2', 'logWErrorDerHalfDimComp3_3'};% plots of the dec log of the error only for W from M2 (ROM dim=2*k0) for the three time steps

eignames={'efile1', 'efile2', 'efile3'};%plot of the decimal log of the singular values in M1 for each of the three steps
eeignames={'eefile1'};%plot of the decimal log of the singular values in M1 and M2 (dim=k1) for all three steps


fV={'VdiscretePlot1','VdiscretePlot2', 'VdiscretePlot3','VdiscretePlot4','VdiscretePlot5', 'VdiscretePlot6','VdiscretePlot7','VdiscretePlot8', 'VdiscretePlot9'};%plots of the FOM and ROM (M1 with dim 2*k0) V solutions at 4 discrete values of X for a given time step 
fw={'WdiscretePlot1','WdiscretePlot2','WdiscretePlot3','WdiscretePlot4','WdiscretePlot5','WdiscretePlot6','WdiscretePlot7','WdiscretePlot8','WdiscretePlot9'};
%plots of the FOM and ROM (M1 with dim 2*k0) W solutions at 4 discrete values of X for a given time step 

  
%%%%%%%%%%%

fprintf(fid, '   n= %3d, L= %3d m= %3d, t0= %1.1f, tfinal= %1.1f tend= %1.1f, N= %3d \n', n, L, m, t0, tfinal, tend, N);
fprintf(fid, 'Data:  \n   s_error= %1.2g, delta1= %1.4f, delta2= %1.4f, a= %1.4f, lambda= %1.4f, gamma= %1.4f \n', s_error, delta1, delta2, a, lambda, gamma);
fprintf(fid, '   r0= %1.4f, r1= %1.4f, xl= %1.4f, xr= %1.4f, mu= %1.4f \n', r0, r1, xl, xr, mu);
fprintf(fid, ' BC y0t='); 

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Calculate the eigenvalues of the Jacobian matrix (T) for the linearized FH-N model %%%%%%%%%%%%%%%
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par1=delta1/dx2;
par2=delta2/dx2;

clear T;
T=zeros(2*(L+1));%initialize a zero matrix T

T(1,1)=-par1*(1+offset1); 
T(1,2)=par1;

for i=2:L
    T(i,i-1)=par1;
    T(i,i)=-2*par1*(1+offset2); %%%% Changing the values of offset 1, 2, 3 will change the stability properties of T 
    T(i,i+1)=par1;
end

T(L+1,L)=-par1*(1+offset3);
T(L+1, L-1)=par1;

for i=2:L
T(L+i+1, L+i)= par2;
T(L+i+1, L+i+1)= -2*par2-gamma;
T(L+i+1,L+i+2)=par2;
T(L+i+1,i+1)=mu;
end


T(L+2,L+2)=-1;
T(2*L+2,2*L+2)=-1;

clear e eV;
eV=eig(T);  %%%  Calculate the eigenvalues of T
e=sort(eV); %%%  Sort the eigenvalues of T

fprintf(feig,'\n');
fprintf(feig,'eig(T)=  \n'); 

for i=1:n
fprintf(feig,'%1.3f\n',e(i)); % output eigenvalues in a file 
end
%%%%%%%%%%%%%%%%%%%%%%%%%% Jacobian T eigevalues have been calculated and output in a file %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Calculate FOM solution and output values at m equidistant points %%%%%%%%%%%%%%%%%%%%%%%

reltol=eps; 
abstol=eps;
options=odeset('RelTol',reltol,'AbsTol',abstol); 

tout=linspace(t0,tend,m); % output m equally spaced values of solution - this output will be used for plotting and comparing with the ROM solution
disp('solving the full ode');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,q]=ode15s(@ODEdef,tout,yinit,options); % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('solved the full ode');

%%%%%%%%%%%%%%%%% LABELS %%%%%%%%%%%%%%%%%%%%%%%%%

se= 'eps= ';
sa=' a=';
sd1=' D1=';
sd2=' D2=';
sla=' lambda=';
sga=' gamma=';
sm='  mu=';
sN1='  N=';
sN2='  N1=';
sl0=' *';
sl1=' .';
sl2=' o';
sx=' x=';
sr0='  r0=';
sr1='  r1=';
sof1=' of1=';sof2=' of2=';sof3=' of3=';

s0 = sprintf('With derivatives!');
s1 = sprintf('y(x,t):  n= %3d eqs, compared to ROM sol. with N snapshots and l1 degrees', n);
s3 = sprintf('w(x,t):  n= %3d eqs, compared to ROM sol. with N snapshots and l1 degrees', n);
s4=sprintf('y0(t)=-r0; y1(t)=-r1, %s %1.2f %s %1.2f, %s %3.3f ', sr0, r0, sr1, r1, sof1, offset1,sof2, offset2, sof3, offset3); 
s12 = sprintf('y(x,t):  n= %3d eqs, compared to ROM sol. with N snapshots and l2 degrees', n);
s32 = sprintf('w(x,t):  n= %3d eqs, compared to ROM sol. with N snapshots and l2 degrees', n);
s_vn = sprintf('v(x,t):  n= %3d eqs', n);
s_wn = sprintf('w(x,t):  n= %3d eqs', n);
s_param=sprintf('%s %1.2f %s %1.2f %s %1.2f %s %1.2f %s %1.2f %s %1.2f', sa, a, sd1, delta1, sd2, delta2, sla, lambda, sga, gamma, sm, mu);

s_init=sprintf('y0(t)=-r0; y1(t)=-r1, %s %1.2f %s %1.2f, %s %3.3f, %s %3.3f', sr0, r0, sr1, r1);
sxx=sprintf('%s %3.1f', sx, xr);

%%%%%%%%%%%% Plot FOM  V,w %%%%%%%%%%%%%%%%

clear v; clear time; clear dist;

for i=1:m
    time(i)=t0+(i-1)*(tend-t0)/m;
    for j=1:n/2
        dist(j)=xl+(j-1)*(xr-xl)/L; % calculated in order to plot waterfall below FOR V
        v(i,j)=q(i,j);
        %disp(i); disp(dist(j)); disp(v(i,j));
    end
end

%{
f1=figure;
mesh(time,dist,v');
waterfall(time,dist,v');
title([{s_vn},{s_init},{s_param}],'fontsize', 9);xlabel('Time');
print(f1,  Vplots{1},'-deps', '-r200');
%}

clear w; clear time; clear dist;

for i=1:m
    time(i)=t0+(i-1)*(tend-t0)/m;
    for j=n/2+1:n
        dist(j-n/2)=xl+(j-n/2-1)*(xr-xl)/L; % calculated in order to plot waterfall below
        w(i,j-n/2)=q(i,j);
        %disp(i); disp(dist(j)); disp(v(i,j));
    end    
end

%{
f2=figure;
mesh(time,dist,w');
waterfall(time,dist,w');
title([{s_wn},{s_init},{s_param}],'fontsize', 9);xlabel('Time');
print(f2,Wplots{1},'-deps', '-r200');
%}
%%%%%%%%%%%% END Plot FOM  V,w %%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Loop to calculate ROMs %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y=3;% number of different between snapshot time intervals to be experimented with 


for k=1:1 %This loop right now passes only once. It can be used to vary parameters of the simulation or of the model

for klm=1:Y %this loop explores the dependence on the size of the time step; we take snapshots at each t0+Deltat*j, j=0, 1...N timestep
    N=N0/(10^(klm-1)); NS(klm)=N; %will reduce number of snapshots with each new test 

    disp('N= '); 
    disp(N);


    Deltat=(tfinal-t0)/N; %Deltat is an interval defining how often we collect snapshots

% calculate solution and output values only at time points with step Deltat
P=eye(n); %initially P is the unit matrix

clear y;
 
tspan=t0:Deltat:tfinal; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[t,y]=ode15s(@ODEdef,tspan,yinit,options); % the output is a set of N+1 equistant snapshots y which will be used o form the ROM basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% now calculate the derivatives at the snapshot time points %%%%%

clear z; %Here will keep the time derivative snapshots

for i=1:N+1 %N is the number of snapshots, n is the dimension of the system
    for j=1:n
        z(i,j)=0;
    end
end

clear ytemp;clear yf;
    
for i=1:N+1
    time=t0+(i-1)*Deltat;
    for j=1:n
        ytemp(j)=y(i,j); % ytemp is a vector that contains snapshots
    end
    
    yf=ODEdef(time,ytemp); % yf is a vector that contains the derivative 
        
    for j=1:n
    z(i,j)=yf(j); % z is a set of N+1 vectors containing the derivatives at the time points at which snapshots were collected
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Calculate the SVD of y transposed - this is because of the way ode15s creates the solution; %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

clear U0; clear S0; clear V0;
[U0,S0,V0]=svd(y');

%[U0,S0,V0]=svds(y', N+1);

%display the left eigenvectors of the SVD, these are of dimension n
%disp(U0);
disp('size of y'); disp(size(y));
disp('size of U0'); disp(size(U0)); 
disp('size of S0'); disp(size(S0)); 

%{
%checking if the vectors are orthonormal
fprintf('U0 orthonormality: \n');
p=U0*U0'; disp('U0*U0'=');disp(p);

%end checking if the vectors are orthonormal
%}

% ---------------- find the dominant singular eigenvalues for a given s_error cutoff --------

%disp('S0:   '); disp(S0); 

clear rr cc;
[rr,cc]=size(S0); %how many columns and rows does the singular values matrix have
k0=min(cc,rr);

for i=1:min(cc,rr)
fprintf(fid,'S0(i i): %d %.3f \n', i,  S0(i,i)); 
end


disp('k0= '); disp k0;

for j=1:k0
   eigv0(klm, j)=log10(S0(j,j)); %since S0(i,i) are sorted, all k0 singular values are not 0
   axis0(klm, j)=j;
end


fprintf(feig,'k0 Singular values =  %3d \n', k0);

for i=1:k0
fprintf(feig,'%1.3f\n',eigv0(klm,i));
end

%{
plot(axis0,eigv0);
title('Method 1')
print('-deps', '-r200', eignames{klm});
%}

disp('k0='); disp(k0); disp('S0(k0,k0)='); disp(S0(k0,k0));
energy=0;

%%%%%%%%%%%%%%%%%% Option for defining the ROM basis dimension by a cutoff of the singular values
%{
while energy < s_error && k0>0
    energy=S0(k0,k0); 
    k0=k0-1;
end

k0=k0+1;
disp('S0(k0,k0)='); disp(S0(k0,k0));
% k0 now is the number of dominant eigenvectors to use

fprintf(fid,'k0= %3d \n',  k0);
%}
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Now forming the ROM basis and projection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear X0; clear P0;
X0=U0(:,1:2*k0); %the reduced basis dim

P0=X0*X0';  % now P0 is the projection matrix - it is k0xk0 


clear X00; clear P00;
X00=U0(:,1:k0); %the reduced basis dim

P00=X00*X00';

%%%%%%%%%%%%%%%%%%%%% Now find the eigenvalues of the projected jacobian matrix %%%%%%%%%%%%%%%%%%%%%%%
clear petV pet PtildeT;
PtildeT=P0*T;
petV=eig(PtildeT);
pet=sort(petV);

fprintf(feig,'eig(PtildeT)=  \n');
disp('eig(PtildeT)=  '); disp (pet);

for i=1:n
fprintf(feig,'%1.3f\n',pet(i));
end

clear petV0 pet0 PtildeT0;
PtildeT0=P00*T;
petV0=eig(PtildeT0);
pet0=sort(petV0);

fprintf(feig,'eig(PtildeT)=  \n');
disp('eig(PtildeT)=  '); disp (pet0);

for i=1:n
fprintf(feig,'%1.3f\n',pet0(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Now solve the lifted ROM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear b0;

clear binit0;
P=P0; % assigning P0 to P because P is global and used in ODEdef
binit0=P*yinit'; % initial conditions for lifetd ROM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,b0]=ode15s(@ODEdef,tout,binit0,options);  % This finds the lifted ROM solution where the ROM basis construction uses only solution snapshots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear b00;

clear binit00;
P=P00; % assigning P0 to P because P is global and used in ODEdef
binit00=P*yinit'; % initial conditions for lifetd ROM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,b00]=ode15s(@ODEdef,tout,binit00,options);  % This finds the lifted ROM solution where the ROM basis construction uses only solution snapshots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% visualize solution of lifted ROM 
NN=N; 
l0(klm)=2*k0; 

clear w0; clear time; clear dist;
for i=1:m
    time(i)=t0+(i-1)*(tend-t0)/m;
    for j=n/2+1:n
        dist(j-n/2)=xl+(j-n/2-1)*(xr-xl)/L; % calculated in order to plot waterfall below
        w0(i,j-n/2)=b0(i,j);
        %disp(i); disp(dist(j)); disp(w0(i,j));
    end    
end

%%%%%%%%%%%%%%%%%%%% Plot the error for a few discrete x %%%%%%%%%%%%%%%%%%%%%%%
%f3=figure;  
%{
plot(time,w(:,1),'k');  
hold on;
plot(time,w(:,L/4),'g');  
hold on;
plot(time,w(:,L/2),'b');  
hold on;
plot(time,w(:,3*L/4),'m');  
hold on;
plot(time,w(:,L),'r'); 
hold on;
plot(time,w0(:,1),'ko');  
hold on;
plot(time,w0(:,L/4),'go');  
hold on;
plot(time,w0(:,L/2),'bo');  
hold on;
plot(time,w0(:,3*L/4),'mo');  
hold on;
plot(time,w0(:,L),'ro'); 
titw1=sprintf('Number of snapshots = %3d, n= %3d, ROM dim.= %3d ', NS(klm)+1, n, 2*k0);  
%titw=sprintf('FOM (-)  and M1 ROM (o) solutions w(0,t)-black, w(X/4,t): g, w(X/2,t): b, w(3X/4,t)- pink, w(X,t)- red');    
titw=sprintf(' FOM (-)  and M1 ROM (o) solutions w(0,t)-black,  w(X/2,t)- pink,  w(X,t)- red');    
title([{titw},{titw1}, {sxx}], 'fontsize', 9);xlabel('Time'); ylabel('w(.,t)');
print (f3, fw{klm},'-depsc', '-r200');
clear f3;
%}

%{
mesh(time,dist,w');
waterfall(time,dist,w');
title([{s_wn},{s_init},{s_param}],'fontsize', 9);xlabel('Time');
print('-deps', '-r200', wplots{k});
%title([{s1}, {s2}, {s4}], 'fontsize', 9); xlabel('Time');  
%}

clear v0; clear time; clear dist;

for i=1:m
    time(i)=t0+(i-1)*(tend-t0)/m;
    for j=1:n/2
        dist(j)=xl+(j-1)*(xr-xl)/L; % calculated in order to plot waterfall below
        v0(i,j)=b0(i,j);
        %disp(i); disp(dist(j)); disp(v0(i,j));
    end
    
end

f4=figure;  
plot(time,v(:,1),'k');  
hold on;
plot(time,v(:,L/4),'g');  
hold on;
plot(time,v(:,L/2),'b');  
hold on;
plot(time,v(:,3*L/4),'m');  
hold on;
plot(time,v(:,L),'r'); 
hold on;
plot(time,v0(:,1),'ko');  
hold on;
plot(time,v0(:,L/4),'go');  
hold on;
plot(time,v0(:,L/2),'bo');  
hold on;
plot(time,v0(:,3*L/4),'mo');  
hold on;
plot(time,v0(:,L),'ro'); 
titv1=sprintf('Number of snapshots = %3d, n= %3d, ROM dim.= %3d ', NS(klm)+1, n, 2*k0);  
%titv=sprintf('FOM (-) and M1 ROM (o) solutions V(0,t): blk, V(X/4,t): g, V(X/2,t): b, V(3X/4,t): m, V(X,t): r');        
titv=sprintf('FOM (-)  and M1 ROM (o) solutions v(0,t)-black,  v(X/2,t)- pink,  v(X,t)- red');    
title([{titv},{titv1}, {sxx}], 'fontsize', 9);xlabel('Time'); ylabel('v(.,t)');
print (f4, fV{klm},'-depsc', '-r200');
clear f4;


%{
s_init=sprintf(', %s %1.2f %s %1.2f, %s %3.3f ', sr0, r0, sr1, r1);
mesh(time,dist,v');
waterfall(time,dist,v');
title([{s_vn},{s_init},{s_param}],'fontsize', 9);xlabel('Time');
print('-deps', '-r200', ffnames{1});

%title([{s1}, {s2}, {s4}], 'fontsize', 9);   
 
%}

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculating the error and the decimal log of the error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('difference b/n appr and exact solution for V \n'); 
fprintf(fid, 'difference b/n appr and exact solution for V \n'); 
 
errorV0=0;
errorW0=0;
error0=0;
error00=0;

%%%% Total error
for i=1:m %recall m = 100; so we calculate the error at each  point Tsart+i*(Tend-Tstart)/m
    
    diff=0; 
    diff0=0; 
    qnorm=0;
    
    for M=1:2*L+2      
    diff=diff+(b0(i,M)-q(i,M))^2; 
    diff0=diff0+(b00(i,M)-q(i,M))^2; 
    qnorm=qnorm+q(i,M)^2;        
    end
    
    Diff0(i)=sqrt(diff);%this is the total L2 error for the first test    
    Log2D0(i)=log2(Diff0(i));%Log2
    LD20(klm,i)= Log2D0(i);
    Log10D0(i)=log10(Diff0(i));%Decimal log of L2 error at time(i)
    LD100(klm,i)=Log10D0(i);

    Diff00(i)=sqrt(diff0);%this is the total L2 error for the first test    
    Log2D00(i)=log2(Diff00(i));%Log2
    LD200(klm,i)= Log2D00(i);
    Log10D00(i)=log10(Diff00(i));%Decimal log of L2 error at time(i)
    LD1000(klm,i)=Log10D00(i)
    
    qnorm=sqrt(qnorm); 
    sol0Norm(i)=qnorm; %L2 norm of exact solution (FOM)  NOT REALLY USED
 
    error0=max(error0, Diff0(i)); %This is the largest error for all i=1, . . , m
    error00=max(error00, Diff00(i));
end

for i=1:m
    fprintf(fid,'Error0: %1.3g %1.3g \n', time(i), Diff0(i));
    fprintf(fid,'Error00: %1.3g %1.3g \n', time(i), Diff00(i));
end

fprintf(fid,'Error0: %1.3g \n', error0); %This is the largest error for all i=1, . . , m;  NOT REALLY USED
fprintf(fid,'Error00: %1.3g \n', error00);
fprintf('difference b/n appr and exact solution for w \n'); 
fprintf(fid, 'difference b/n appr and exact solution for w \n'); 

%%%%%% Error in V
for i=1:m   
    diff=0;
    qnorm=0;
    
    for M=1:L+1
    diff=diff+(b0(i,M)-q(i,M))^2; 
    qnorm=qnorm+q(i,M)^2;
    end
    
    qnorm=sqrt(qnorm);
    DiffV0(i)=sqrt(diff);%this is the error for the first test
    Log2DV0(i)=log2(DiffV0(i));
    LD2V0(klm,i)=Log2DV0(i); 
    Log10DV0(i)=log10(DiffV0(i));%Decimal log of L2 error at time(i)
    LD10V0(klm,i)=Log10DV0(i);
   
    solV0Norm(i)=qnorm; 
    
    errorV0=max(errorV0, DiffV0(i));
        
end

for i=1:m
    fprintf(fid,'Error V0: %1.3g %1.3g \n', time(i), DiffV0(i));
end

fprintf(fid,'Error V0: %1.3g \n', errorV0);
fprintf('difference b/n appr and exact solution for w \n'); 
fprintf(fid, 'difference b/n appr and exact solution for w \n'); 


%%%%%% Error in w
for i=1:m
    diff=0;
    qnorm=0;
    
    for M=L+2:2*L+2 
    diff=diff+(b0(i,M)-q(i,M))^2;
    qnorm=qnorm+q(i,M)^2;
   
    end
    
    qnorm=sqrt(qnorm);
    Diffw0(i)=sqrt(diff);
    Log2DW0(i)=log2(Diffw0(i));
    LD2W0(klm,i)=Log2DW0(i);
    Log10DW0(i)=log10(Diffw0(i));
    LD10W0(klm,i)=Log10DW0(i);
    
    solW0Norm(i)=qnorm;
   
    errorW0=max(errorW0, Diffw0(i));
end


for i=1:m
    fprintf(fid,'Error W0: %1.3g %1.3g \n', time(i), Diffw0(i));
end

fprintf(fid,'Error W0: %1.3g \n', errorW0);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%%%%%%%%%%%%%% Now solve with derivatives and solution vectors as snapshots
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clear c1;clear U1; clear S1; clear V1;

c1=vertcat(y,z); %here derivative and ordinary snapshots are taken at N1<=N points

disp('Now solve with solution vectors and derivatives vectors  as snapshots');
disp('size of y ');disp(size(y));disp('size of z ');disp(size(z));disp('size of c1 ');disp(size(c1));
[U1,S1,V1]=svd(c1'); 

disp('size of U1'); disp(size(U1));
disp('size of S1'); disp(size(S1));

%fprintf('U1 orthonormality: \n');
%p=U01*U01'; disp('U1*U1*='); disp(p);
%display the left eigenvectors of the SVD, these are of dimension n

clear rr cc;
[row1,col1]=size(U1);
[rr,cc]=size(S1);

k1=min(rr,cc); 

for i=1:min(cc,rr)
fprintf(fid,'S1(i i): %d %.3f \n', i,  S1(i,i)); 
end


for j=1:k1
    eigv1(klm,j)=log10(S1(j,j));%Since S1 is sorted, all k1 singular values are non zero
    axis1(klm,j)=j;
end

fprintf(feig,'k1 Singular values =  %3d \n', k1);

for i=1:k1
fprintf(feig,'%1.3f\n',eigv1(klm,i));
end

%{
energy=0;

while energy< s_error && k1>0
    energy= S1(k1,k1); 
    k1=k1-1; 
end

k1=k1+1; 

disp('k1=');disp(k1); 
disp('S1(k1,k1)='); 
disp(S1(k1,k1));
% k1 now is the number of dominant eigenvectors to use

fprintf(fid,'k1= %3d \n',k1);
%}
%DISPLAYING partially THE FIRST 10 vectors of U1
%disp('U, U1 first 10 columns'); disp(U1(1:10,1:10));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear X1;
X1=U1(:,1:k0);  %the first k0 vectors as basis

clear P1;
P1=X1*X1';

%%%%%%%%%%%%%%Eigenvalues of projection
clear pet1 petV PtildeT1;
PtildeT1=P1*T; % Projection of Jacobian
petV=eig(PtildeT1);
pet1=sort(petV);

fprintf(feig,'eig(PtildeT1)=  \n');

for i=1:n
fprintf(feig,'%1.3f\n',pet1(i));
end

P=P1;

clear binit1; clear b1;
binit1=P*yinit';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,b1]=ode15s(@ODEdef,tout,binit1,options);   % Finds the lifted ROM solution where the ROM basis is constructed using solution and time derivative snapshots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% 
clear X2;
X2=U1(:,1:2*k0);  %the first k0 vectors as basis to compare error of 2 methods for same dimension

clear P2;
P2=X2*X2';

%%%%%%%%%%%%%%Eigenvalues of projection
clear petV pet2 PtildeT2;
PtildeT2=P2*T;
petV=eig(PtildeT2);
pet2=sort(petV);

fprintf(feig,'eig(PtildeT2)=  \n');

for i=1:n
fprintf(feig,'%1.3f\n',pet2(i));
end


P=P2;

clear binit2; clear b2;
binit2=P*yinit';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,b2]=ode15s(@ODEdef,tout,binit2,options); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualize solution of ROMs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% b1 %%%%%%%%%%%%%%%%%

clear w1; clear time; clear dist;
for i=1:m
    time(i)=t0+(i-1)*(tend-t0)/m;
    for j=n/2+1:n
        dist(j-n/2)=xl+(j-n/2-1)*(xr-xl)/L; % calculated in order to plot waterfall below
        w1(i,j-n/2)=b1(i,j);
        %disp(i); disp(dist(j)); disp(v(i,j));
    end    
end

%f5=figure;  
%{
plot(time,w(:,1),'k');  
hold on;
plot(time,w(:,L/4),'g');  
hold on;
plot(time,w(:,L/2),'b');  
hold on;
plot(time,w(:,3*L/4),'m');  
hold on;
plot(time,w(:,L),'r'); 
hold on;
plot(time,w1(:,1),'ko');  
hold on;
plot(time,w1(:,L/4),'go');  
hold on;
plot(time,w1(:,L/2),'bo');  
hold on;
plot(time,w1(:,3*L/4),'mo');  
hold on;
plot(time,w1(:,L),'ro'); 
titw21=sprintf('Number of snapshots = %3d, n= %3d, ROM dim.= %3d ', NS(klm)+1, n, k1);  
%titw20=sprintf('Method 2. FOM (-)  and ROM (*) solutions w(0,t): black, w(X/4,t): g, w(X/2,t): b, w(3X/4,t): m, w(X,t): r');    
%title([{titw20},{titw21}], 'fontsize', 9);xlabel('Time');
%titw20=sprintf('FOM (-)  and M2 ROM (o) solutions w(0,t)-black, w(X/4,t): g, w(X/2,t): b, w(3X/4,t)- pink, w(X,t)- red');    
titw20=sprintf('FOM (-)  and M2 ROM (*) solutions w(0,t)-black,  w(X/2,t)- pink,  w(X,t)- red');    
title([{titw20},{titw21}, {sxx}], 'fontsize', 9);xlabel('Time'); ylabel('w(.,t)');
print (f5, fw{3+klm},'-depsc', '-r200');
clear f5;
%}

clear v1; clear time; clear dist;
for i=1:m
    time(i)=t0+(i-1)*(tend-t0)/m;
    for j=1:n/2
        dist(j)=xl+(j-1)*(xr-xl)/L; %calculated in order to plot waterfall below
        v1(i,j)=b1(i,j);
    end
    
end

f6=figure;  
plot(time,v(:,1),'k');  
hold on;
plot(time,v(:,L/4),'g');  
hold on;
plot(time,v(:,L/2),'b');  
hold on;
plot(time,v(:,3*L/4),'m');  
hold on;
plot(time,v(:,L),'r'); 
hold on;
plot(time,v1(:,1),'ko');  
hold on;
plot(time,v1(:,L/4),'go');  
hold on;
plot(time,v1(:,L/2),'bo');  
hold on;
plot(time,v1(:,3*L/4),'mo');  
hold on;
plot(time,v1(:,L),'ro');
titv21=sprintf('Number of snapshots = %3d, n= %3d, ROM dim.= %3d ', NS(klm)+1, n, k1);  
%titv20=sprintf('Method 2 FOM (-)  and ROM (o) solutions V(0,t): black, V(X/4,t): g, V(X/2,t): b, V(3X/4,t): m, V(X,t) - r');    
%title([{titv20},{titv21}], 'fontsize', 9);xlabel('Time');
%titv20=sprintf(' FOM (-)  and M2 ROM (o) solutions v(0,t)-black, v(X/4,t): g, v(X/2,t): b, v(3X/4,t)- pink, v(X,t)- red');    
titv20=sprintf('FOM (-)  and M2 ROM (o) solutions. v(0,t)-black,  v(X/2,t)- pink,  v(X,t)- red');    
title([{titv20},{titv21}, {sxx}], 'fontsize', 9);xlabel('Time'); ylabel('v(.,t)');
print (f6, fV{3+klm},'-depsc', '-r200');
clear f6;

%%%%%%%%%%%%%%%%%%%% b2 %%%%%%%%%%%%%%%%
clear w2; clear time; clear dist;
for i=1:m
    time(i)=t0+(i-1)*(tend-t0)/m;
    for j=n/2+1:n
        dist(j-n/2)=xl+(j-n/2-1)*(xr-xl)/L; % calculated in order to plot waterfall below
        w2(i,j-n/2)=b2(i,j);
        %disp(i); disp(dist(j)); disp(v(i,j));
    end    
end

%f7=figure;  
%{
plot(time,w(:,1),'k');  
hold on;
plot(time,w(:,L/4),'g');  
hold on;
plot(time,w(:,L/2),'b');  
hold on;
plot(time,w(:,3*L/4),'m');  
hold on;
plot(time,w(:,L),'r');
hold on;
plot(time,w2(:,1),'ko');  
hold on;
plot(time,w2(:,L/4),'go');  
hold on;
plot(time,w2(:,L/2),'bo');  
hold on;
plot(time,w2(:,3*L/4),'mo');  
hold on;
plot(time,w2(:,L),'ro'); 
titw22=sprintf('Number of snapshots = %3d, n= %3d, ROM dim.= %3d ', NS(klm)+1, n, 2*k0);  
%titw20=sprintf('Method 2. FOM (-)  and ROM (*) solutions w(0,t): black, w(X/4,t): g, w(X/2,t): b, w(3X/4,t): m, w(X,t): r');    
%title([{titw20},{titw22}], 'fontsize', 9);xlabel('Time');
%titw20=sprintf('FOM (-)  and M2 ROM (o) solutions w(0,t)-black, w(X/4,t): g, w(X/2,t): b, w(3X/4,t)- pink, w(X,t)- red');    
titw20=sprintf('FOM (-)  and M2 ROM (o) solutions. w(0,t)-black,  w(X/2,t)- pink,  w(X,t)- red');    
title([{titw20},{titw22}, {sxx}], 'fontsize', 9);xlabel('Time'); ylabel('w(.,t)');
print (f7,fw{6+klm}, '-depsc', '-r200');
clear f7;
%}

clear v2; clear time; clear dist;
for i=1:m
    time(i)=t0+(i-1)*(tend-t0)/m;
    for j=1:n/2
        dist(j)=xl+(j-1)*(xr-xl)/L; %calculated in order to plot waterfall below
        v2(i,j)=b2(i,j);
    end
    
end

f8=figure;  
plot(time,v(:,1),'k');  
hold on;
plot(time,v(:,L/4),'g');  
hold on;
plot(time,v(:,L/2),'b');  
hold on;
plot(time,v(:,3*L/4),'m');  
hold on;
plot(time,v(:,L),'r');
hold on;
plot(time,v2(:,1),'ko');  
hold on;
plot(time,v2(:,L/4),'go');  
hold on;
plot(time,v2(:,L/2),'bo');  
hold on;
plot(time,v2(:,3*L/4),'mo');  
hold on;
plot(time,v2(:,L),'ro');
titv22=sprintf('Number of snapshots = %3d, n= %3d, ROM dim.= %3d ', NS(klm)+1, n, 2*k0);  
%titv20=sprintf('Method 2 FOM (-)  and ROM (o) solutions V(0,t): black, V(X/4,t): g, V(X/2,t): b, V(3X/4,t): m, V(X,t) - r');    
%title([{titv20},{titv22}], 'fontsize', 9);xlabel('Time');
%titv20=sprintf('FOM (-)  and M2 ROM (o) solutions v(0,t)-black, v(X/4,t): g, v(X/2,t): b, v(3X/4,t)- pink, v(X,t)- red');    
titv20=sprintf('FOM (-)  and M2 ROM (o) solutions. v(0,t)-black,  v(X/2,t)- pink,  v(X,t)- red');    
title([{titv20},{titv22}, {sxx}], 'fontsize', 9);xlabel('Time'); ylabel('v(.,t)');
print (f8, fV{6+klm},'-depsc', '-r200');
clear f8;

%{
clear w; clear time; clear dist;
for i=1:m
    time(i)=t0+(i-1)*(tend-t0)/m;
    for j=n/2+1:n
        dist(j-n/2)=xl+(j-n/2-1)*(xr-xl)/L; % calculated in order to plot waterfall below
        w(i,j-n/2)=b1(i,j);
        %disp(i); disp(dist(j)); disp(v(i,j));
    end

    
end

mesh(time,dist,w');
waterfall(time,dist,w');
title([{s_wn},{s_init},{s_param}],'fontsize', 9);xlabel('Time');
print('-deps', '-r200', wplots{k});
%}

l1(klm)= k0;
l2(klm)= 2*k0;


fprintf(fid, 'difference b/n appr and exact solution for V; k1 \n'); 

errorV1=0;
errorW1=0;  
error1=0;
errorV2=0;
errorW2=0;  
error2=0;

for i=1:m
    diff=0;
    qnorm=0; 
    
    for M=1:2*L+2
    diff=diff+(b1(i,M)-q(i,M))^2;
    qnorm=qnorm+q(i,M)^2;
    end
    
qnorm=sqrt(qnorm);%norm of the exact solution - we are calculating this many times!

Diff1(i)=sqrt(diff);
Log2D1(i)=log2(Diff1(i));
LD21(klm,i)=Log2D1(i);
Log10D1(i)=log10(Diff1(i));
LD101(klm,i)=Log10D1(i);

sol1Norm(i)=qnorm; 
error1=max(error1, Diff1(i));    

end

for i=1:m %here calculate only the error in V
    diff=0;
    qnorm=0;
    
    for M=1:L+1
    diff=diff+(b1(i,M)-q(i,M))^2;
    qnorm=qnorm+q(i,M)^2;   
    end
    
qnorm=sqrt(qnorm);
DiffV1(i)=sqrt(diff);
Log2DV1(i)=log2(DiffV1(i));
LD2V1(klm,i)=Log2DV1(i);
Log10DV1(i)=log10(DiffV1(i));
LD10V1(klm,i)=Log10DV1(i);

solV1Norm(i)=qnorm; 
errorV1=max(errorV1, DiffV1(i));    
end


for i=1:m %for w
    fprintf(fid,'Error V1: %1.3g %1.3g \n', time(i), DiffV1(i));
end
fprintf(fid,'Error V1: %1.3g \n', errorV1);


fprintf('difference b/n appr and exact solution for w \n'); 
fprintf(fid, 'difference b/n appr and exact solution for w \n'); 

for i=1:m
    diff=0;
    qnorm=0;
    
    for M=L+2:2*L+2
    diff=diff+(b1(i,M)-q(i,M))^2;
    qnorm=qnorm+q(i,M)^2;   
    end
    
    Diffw1(i)=diff;
    Log2DW1(i)=log2(Diffw1(i));
    LD2W1(klm,i)=Log2DW1(i);
    Log10DW1(i)=log10(Diffw1(i));
    LD10W1(klm,i)=Log10DW1(i);
    
    qnorm=sqrt(qnorm);
    solW1Norm(i)=qnorm; 
    errorW1=max(errorW1, Diffw1(i));   
end


for i=1:m
    fprintf(fid,'Error W1: %1.3g %1.3g \n', time(i), Diffw1(i));
end

fprintf(fid,'Error W1: %1.3g \n', errorW1);


%%%%%%

fprintf(fid, 'difference b/n appr and exact solution for V; 2*k0 \n'); 


for i=1:m
    diff=0;
    qnorm=0; 
    
    for M=1:2*L+2
    x=xl+(M-1)*(xr-xl)/L;   
    diff=diff+(b2(i,M)-q(i,M))^2;
    qnorm=qnorm+q(i,M)^2;
    end
    
qnorm=sqrt(qnorm);%norm of the exact solution 

Diff2(i)=sqrt(diff);
Log2D2(i)=log2(Diff2(i));
LD22(klm,i)=Log2D2(i);
Log10D2(i)=log10(Diff2(i));
LD102(klm,i)=Log10D2(i);

sol2Norm(i)=qnorm; 

error2=max(error2, Diff2(i));    

end

for i=1:m %here calculate only the error in V
    diff=0;
    qnorm=0;
    
    for M=1:L+1
    x=xl+(M-1)*(xr-xl)/L;
    diff=diff+(b2(i,M)-q(i,M))^2;
    qnorm=qnorm+q(i,M)^2;
    
    end
    
qnorm=sqrt(qnorm);
DiffV2(i)=sqrt(diff);
Log2DV2(i)=log2(DiffV2(i));
LD2V2(klm,i)=Log2DV2(i);
Log10DV2(i)=log10(DiffV2(i));
LD10V2(klm,i)=Log10DV2(i);

solV2Norm(i)=qnorm; 

errorV2=max(errorV2, DiffV2(i));    

end


for i=1:m %for w
    fprintf(fid,'Error V2: %1.3g %1.3g \n', time(i), DiffV2(i));
end
fprintf(fid,'Error V2: %1.3g \n', errorV2);

fprintf(fid, 'difference b/n appr and exact solution for w \n'); 

for i=1:m
    diff=0;
    qnorm=0;
    
    for M=L+2:2*L+2
    x=xl+(M-L-2)*(xr-xl)/L;         
    diff=diff+(b2(i,M)-q(i,M))^2;
    qnorm=qnorm+q(i,M)^2;
    
    end
    
    Diffw2(i)=diff;
    Log2DW2(i)=log2(Diffw2(i));
    LD2W2(klm,i)=Log2DW2(i);
    Log10DW2(i)=log10(Diffw2(i));
    LD10W2(klm,i)=Log10DW2(i);
    
    qnorm=sqrt(qnorm);
    solW2Norm(i)=qnorm; 
    fprintf(fid, 'x:  %2.2f   %1.3g %1.5g  \n', x, Diffw2(i), solW2Norm(i));  
    errorW2=max(errorW2, Diffw2(i));   
end


for i=1:m
    fprintf(fid,'Error W2: %1.3g %1.3g \n', time(i), Diffw2(i));
end

fprintf(fid,'Error W2: %1.3g \n', errorW2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'Differences in the error between methods \n');
fprintf(fid,'Error V0, V1 and V2: %1.3g %1.3g  \n', errorV0, errorV1);
fprintf(fid,'Error W0, W1 and W2: %1.3g %1.3g  \n', errorW0, errorW1);
fprintf(fid,'Error V0, V1 and V2: %1.3g %1.3g  \n', errorV0, errorV2);
fprintf(fid,'Error W0, W1 and W2: %1.3g %1.3g  \n', errorW0, errorW2);    

end

s4=sprintf('y0(t)=-r0; y1(t)=-r1, %s %1.2f %s %1.2f, %s %3.3f %s %3.3f %s %3.3f', sr0, r0, sr1, r1, sof1, offset1,sof2, offset2, sof3, offset3); 
s7=sprintf('%s %1.2f %s %1.2f %s %1.2f %s %1.2f %s %1.2f %s %1.2f', sa, a, sd1, delta1, sd2, delta2, sla, lambda, sga, gamma, sm, mu);
s710=sprintf('M1 (blue) ROM dimensions: %s %3d, %s %3d, %s %3d', sl0, l1(1), sl1, l1(2), sl2, l1(3));
s71=sprintf('M1 (blue) ROM dimensions: %s %3d, %s %3d, %s %3d', sl0, l0(1), sl1, l0(2), sl2, l0(3));
s72=sprintf('M2 (red) ROM dimensions: %s %3d, %s %3d, %s %3d', sl0, l1(1), sl1, l1(2), sl2, l1(3));
s73=sprintf('M2 (red) ROM dimensions: %s %3d, %s %3d, %s %3d', sl0, l2(1), sl1, l2(2), sl2, l2(3));

s8=sprintf('M2 errors for different values of N. N1= %3d, N2= %3d, N3= %3d, n= %3d ', NS(1), NS(2), NS(3), n);
s9=sprintf('M1 errors for different values of N. N1= %3d, N2= %3d, N3= %3d, n= %3d ', NS(1), NS(2), NS(3), n);
s9V=sprintf('v(t): M1 errors for different values of N. N1= %3d, N2= %3d, N3= %3d, n= %3d ', NS(1), NS(2), NS(3), n);
s9w=sprintf('w(t): M1 errors for different values of N. N1= %3d, N2= %3d, N3= %3d, n= %3d ', NS(1), NS(2), NS(3), n);
s990V=sprintf('v(t): M2 errors fordifferent values of N. N1= %3d, N2= %3d, N3= %3d, n= %3d ', NS(1), NS(2), NS(3), n);
s990w=sprintf('w(t): M2 errors for different 3 snapshot numbers: N1= %3d, N2= %3d, N3= %3d, n= %3d ', NS(1), NS(2), NS(3), n);

s99=sprintf('M1 and M2 errors for different values of N.  N1= %3d, N2= %3d, N3= %3d', NS(1), NS(2), NS(3));
s991=sprintf('Parameters: N1= %3d,n= %3d ', NS(1), n);
s992=sprintf('Parameters:  N2= %3d,  n= %3d ', NS(2), n);
s993=sprintf('Parameters:  N3= %3d, n= %3d ', NS(3), n);
s78=sprintf('%s %1.2f %s %1.2f %s %1.2f %s %1.2f %s %1.2f %s %1.2f %s %3d, %s %3d, %s %3d', sa, a, sd1, delta1, sd2, delta2, sla, lambda, sga, gamma, sm, mu);
ssv=sprintf('Distribution of the singular values for 3 snapshot numbers')


f9=figure;
plot(time,LD1000(1,:),'b*', time, LD1000(2,:), 'b.', time, LD1000(3,:), 'bo', time,LD101(1,:),'r*', time, LD101(2,:), 'r.', time, LD101(3,:), 'ro');
title([{s99},{s710},{s72}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f9,fnames{k},'-depsc', '-r200');

f10=figure;
plot(time,LD100(1,:),'b*', time, LD100(2,:), 'b.', time, LD100(3,:), 'bo', time,LD102(1,:),'r*', time, LD102(2,:), 'r.', time, LD102(3,:), 'ro');
title([{s99},{s71}, {s73}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f10, fnames2{k},'-depsc', '-r200');

%{
f11=figure;
plot(time,LD100(1,:),'b*', time, LD100(2,:), 'b.', time, LD100(3,:), 'bo');%Dec Log of Total L2 error Method 1
title([{s9},{s71}],'fontsize', 9);xlabel('Time'); ylabel('Log10 of the L2 error');
print(f11,ggnames{k},'-depsc', '-r200');

f12=figure;
plot(time,LD10V0(1,:),'b*', time, LD10V0(2,:), 'b.', time, LD10V0(3,:), 'bo');%Dec Log of V L2 error Method 1
title([{s9V},{s71}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f12,gvnames{k},'-depsc', '-r200');

f13=figure;
plot(time,LD10W0(1,:),'b*', time, LD10W0(2,:), 'b.', time, LD10W0(3,:), 'bo');%Dec Log of w L2 error Method 1
title([{s9w},{s71}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f13, gwnames{k}, '-depsc', '-r200');

f14=figure;
plot(time,LD101(1,:),'r*', time, LD101(2,:), 'r.', time, LD101(3,:), 'ro');%Dec Log of Total L2 error Method 2
title([{s8},{s72}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f14, fdnames{k}, '-depsc', '-r200');

f15=figure;
plot(time,LD10V1(1,:),'r*', time, LD10V1(2,:), 'r.', time, LD10V1(3,:), 'ro');%Dec Log of V L2 error Method 2
title([{s990V},{s72}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f15,fvnames{k}, '-depsc', '-r200');

f16=figure;
plot(time,LD10W1(1,:),'r*', time, LD10W1(2,:), 'r.', time, LD10W1(3,:), 'ro');%Dec Log of 2 L2 error Method 2
title([{s990w},{s72}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f16, fwnames{k},'-depsc', '-r200');

f17=figure;
plot(time,LD102(1,:),'r*', time, LD102(2,:), 'r.', time, LD102(3,:), 'ro');%Dec Log of Total L2 error Method 2, 1/2 dim
title([{s8},{s71}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f17, fdnames2{k},'-depsc', '-r200');

f18=figure;
plot(time,LD10V2(1,:),'r*', time, LD10V2(2,:), 'r.', time, LD10V2(3,:), 'ro');%Dec Log of V L2 error Method 2, 1/2 dim
title([{s990V},{s71}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f18, fvnames2{k}, '-depsc', '-r200');

f19=figure;
plot(time,LD10W2(1,:),'r*', time, LD10W2(2,:), 'r.', time, LD10W2(3,:), 'ro');%Dec Log of 2 L2 error Method 2, 1/2 dim
title([{s990w},{s71}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f19, fwnames2{k},'-depsc', '-r200');

f20=figure;
plot(time,LD100(1,:),'b*',  time,LD101(1,:),'r*');
title([{s991},{s71}, {s72}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f20, ffnames{k},'-depsc', '-r200');

f21=figure;
plot(time,LD100(2,:),'b.',  time,LD101(2,:),'r.');
title([{s992},{s71}, {s72}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f21, fqnames{k}, '-depsc', '-r200');

f22=figure;
plot(time,LD100(3,:),'bo',  time,LD101(3,:),'ro');
title([{s993},{s71}, {s72}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f22,ddnames{k}, '-depsc', '-r200');

f23=figure;
plot(time,LD100(1,:),'b*',  time,LD102(1,:),'r*');
title([{s991},{s71}, {s73}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f23, ffnames2{k}, '-depsc', '-r200');

f24=figure;
plot(time,LD100(2,:),'b.',  time,LD102(2,:),'r.');
title([{s992},{s71}, {s73}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f24,fqnames2{k}, '-depsc', '-r200');

f25=figure;
plot(time,LD100(3,:),'bo',  time,LD102(3,:),'ro');
title([{s993},,{s71}, {s73}],'fontsize', 9);xlabel('Time');ylabel('Log10 of the L2 error');
print(f25, ddnames2{k}, '-depsc', '-r200');

%}

f26=figure;
plot(axis0(1,:),eigv0(1,:),'r.',axis1(1,:),eigv1(1,:), 'rx', axis0(2,:),eigv0(2,:),'b.',axis1(2,:),eigv1(2,:), 'bx', axis0(3,:),eigv0(3,:),'g.',axis1(3,:),eigv1(3,:), 'gx');
title([{ssv}],'fontsize', 9); 
xlabel('Singular value number');ylabel('Log10 of singular value')
print(f26, eeignames{1},'-depsc', '-r200');

end

fclose('all');
