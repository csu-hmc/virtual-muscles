function [xnew,M,F]=timestep(x,u,h,q)

global muscle_data
%-------------------------------------------------------------
%Create zero initial conditions if called for the first time
%-------------------------------------------------------------
    global xdot_dyn u_dyn
    if ~(exist('xdot_dyn')==1) || isempty(xdot_dyn)
		xdot_dyn = zeros(size(x));
    end
    
    if ~(exist('u_dyn')==1) || isempty(u_dyn)
		u_dyn = zeros(size(u));
    end
%-----------------------------
%1st Order Rosenbrock Method
%-----------------------------
    %Obtain Muscle Lengths
        [Poly,Lm]=moment_arm(q);
    %Evaluate dynamics in current x and xdot
        [f,df_dx,df_dxdot,df_du,F]=muscle_dynamics(x,xdot_dyn,u_dyn,Lm);
        %Solve for the change in x
            du = u - u_dyn;
            dx = (df_dx + df_dxdot/h)\(df_dxdot*xdot_dyn - f - df_du*du);
            %keyboard
            xnew = x + dx;
        %Update Global Variables
            xdot_dyn = dx/h;
            u_dyn = u;
        %Calculate Torque at Joints   
            M=Poly'*F;
end