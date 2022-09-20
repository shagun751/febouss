%% %%
% Reference
% Orszaghova, J., Borthwick, A.G.L., Taylor, P.H., 2012. 
% From the paddle to the beach - A Boussinesq shallow water 
% numerical wave tank based on Madsen and Sorensen equations. 
% Journal of Computational Physics 231, 328–344. 
% https://doi.org/10.1016/j.jcp.2011.08.028
%%%%%

%%
clc
clear all
close all

B = 1/15;
g = 9.81;
h = 1;
A = 0.2;

dx = 0.001;
xmax = 25.4;

C = getC(g,A,h);

qmax = C*A;

opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');
[x,q] = ode113(@(x,q) odefnc(x,q), [0:dx:xmax], [qmax;0],opts);

figure(1)
plot(x,q(:,1))

sx = [-xmax:dx:xmax]';
sp = flipud(q(:,1));
sz = max(size(sp));
sp(end:end+sz-1) = q(:,1);
se = sp/C;

figure(2)
plot(sx,sp);

figure(3)
plot(sx,se);

save('solitaryWave_dx0001.mat')

dlmwrite('solitaryWave_dx0001.dat',[sx,se,sp],'delimiter','\t','precision','%.6f')




%%
function C = getC(g,A,h)
    
    tmp1 = g*h*A*A*(A+3*h);
    tmp2 = 6*h*h*(A - h*log(1+A/h));
    
    C = sqrt(tmp1/tmp2);

end

%%
function dq = odefnc(x,q)

    B = 1/15;
    g = 9.81;
    h = 1;
    A = 0.2;
    C = 3.44066076038234;
    
    dq=zeros(2,1);

    tmp1 = (g*h-C*C)/C*q(1);
    tmp2 = g*q(1)*q(1)/2/C/C;
    tmp3 = C*q(1)*q(1)/(C*h + q(1));
    
    tmp5 = ( B*g*h*h*h/C - C*(B+1/3)*h*h );
    
    tmp6 = (tmp1 + tmp2 + tmp3) / tmp5;
   
    dq(1) = q(2);
    dq(2) = tmp6;

end