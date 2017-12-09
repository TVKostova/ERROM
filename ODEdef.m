%Copyright (c) <2017>, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory
%Written by Tanya Kostova Vassilevska, kostova@llnl.gov, tan.v.kos@gmail.com
%LLNL Release number LLNL-CODE-735916
%All rights reserved.

%This file is part of <ERROM>. For details, see the comments. 

%Licensed under the Apache License, Version 2.0 (the “Licensee”); you may not use this file except in compliance with the License.  You may obtain a copy of the License at:  http://www.apache.org/licenses/LICENSE-2.0

%Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the license.




%This is a routine accompanying the main ERROM Matlab code, which describes the discretization of the Fitzhugh-Nagumo system with 1D diffusion

function yt = FHNpdeA(t,y)
global delta1 delta2 a lambda gamma r0 r1 P L mu dx2 dx offset1 offset2 offset3

%n=length(y);
%L=n/2-1;

%y0t=-r0;
%y1t=-r1;

y0t=-r0; % this is the current at left
y1t=-r1; % this is the current at right

%offset patameters are meant to control the stability 
yt(1)=delta1/dx2*(-dx*y0t -y(1)*(1+offset1)+ y(2)) +lambda*(y(1)*(y(1)-a)*(1-y(1))-y(L+2));

for i=2:L
    yt(i)=delta1/dx2*(y(i-1)-2.0*y(i)*(1+offset2)+y(i+1))+lambda*(y(i)*(y(i)-a)*(1-y(i))-y(L+i+1));
end

yt(L+1)=delta1/dx2*(dx*y1t-y(L)*(1+offset3)+y(L-1))+lambda*(y(L+1)*(y(L+1)-a)*(1-y(L+1))-y(2*(L+1)));

yt(L+2)=0;% derivative is 0 because of constant recovery at left


for i=2:L
    yt(L+i+1)=delta2/dx2*(y(L+i+2)-2*y(L+i+1)+y(L+i))+ mu*y(i+1)-gamma*y(L+i+1);
end


yt(2*L+2)=0; % derivative is 0 because of constant recovery at right

yt=P*yt';

end


