%finite difference method

clear
clc

%Declare Variables
    x = rand(12,1);
    xdot = rand(12,1);
    u = rand(6,1);
    Lm = rand(6,1);
%Compare Derivatives from muscle_dynamics.m
    [f,df_dx_m,df_dxdot_m,df_du_m,F]=muscle_dynamics(x,xdot,u,Lm);
%Preallocate
    df_dx = zeros(12,12);
    df_dxdot = zeros(12,12);
    df_du = zeros(12,6);
    hh=1e-7;
%Start Loop
    for i = 1:12
        xsave = x;
        x(i) = x(i) + hh;
        df_dx(:,i) = (muscle_dynamics(x,xdot,u,Lm) - f)/hh;
        x = xsave;
        xdotsave = xdot;
        xdot(i) = xdot(i)+hh;
        df_dxdot(:,i) = (muscle_dynamics(x,xdot,u,Lm) - f)/hh;
        xdot = xdotsave;
        if i < 7
            usave = u;
            u(i) = u(i)+hh;
            df_du(:,i)= (muscle_dynamics(x,xdot,u,Lm) - f)/hh;
            u = usave;
        end  
    end  