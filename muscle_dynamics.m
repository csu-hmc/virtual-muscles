function [ff ,df_dx, df_dxdot,df_du,force]=muscle_dynamics(x,xdot,u,Lm)

%Define Variables
    HillA=0.25;
    global muscle_data
    %Initialize Loop
    Nmus = size(Lm,1);
    ff=zeros(2*Nmus,1); 
    df_dx=zeros(2*Nmus,2*Nmus);
    df_dxdot=zeros(2*Nmus,2*Nmus);
    df_du=zeros(2*Nmus,Nmus);
    force=zeros(Nmus,1);
    for i=1:Nmus
        %Calculate LCE and cos(p) and derivatives (with respect to s)
            s=x(2*i-1); a=x(2*i);
            sdot=xdot(2*i-1); adot=xdot(2*i);
            b = sin(muscle_data.PennOpt(i));		
            Lce = sqrt(s^2 + b^2);
            cosp = s/Lce;
            dLce_ds = cosp;
            dcosp_ds = b^2/((s^2+b^2)^(3/2));
            
        %Calculate Lcedot and derivatives (with respect to s and sdot)
            Lcedot = sdot*cosp;
            dLcedot_dsdot = cosp;
            dLcedot_ds = sdot*dcosp_ds;
	
        %Normalized Force-Length Relationship at Maximum Activation (F1)
            f = (Lce - 1.0)/muscle_data.Width(i);
            F1 = exp(-f*f);
            dF1_dLce = (-2.0*f*F1)/muscle_data.Width(i);
            dF1_ds = dF1_dLce * dLce_ds;
        %Normalized Force-Velocity Relationship (F2) 
            if Lcedot < 0
                f = muscle_data.Vmax(i) - Lcedot/HillA;
                F2 = (muscle_data.Vmax(i) + Lcedot)/f;
                dF2_dLcedot = (1.0 + F2/HillA)/f;
            else
                f = Lcedot + muscle_data.c(i);
                F2 = (muscle_data.gmax(i)*Lcedot + muscle_data.c(i))/f;
                dF2_dLcedot = (muscle_data.gmax(i) - F2)/f;
            end
            dF2_dsdot =  dF2_dLcedot * dLcedot_dsdot;
            dF2_ds = dF2_dLcedot * dLcedot_ds;
            
        %PEE Force-Length Relationship (F3)
            k1 = 10.0/(muscle_data.Fmax(i)*muscle_data.Lceopt(i));		
            f = (Lce - muscle_data.PEEslack(i));		
            F3 = k1*f;						
            dF3_dLce = k1;
            if (f>0) 												
                F3 = F3 + muscle_data.kPEE(i)*f^2;
                dF3_dLce = dF3_dLce + 2*muscle_data.kPEE(i)*f;
            end
            dF3_ds = dF3_dLce * dLce_ds;
	
        %SEE Force-Length Relationship (F4) 
            k2 = 10.0/muscle_data.Fmax(i);				
            f = Lm(i) - s * muscle_data.Lceopt(i) - muscle_data.SEEslack(i);
            F4 = k2*f;								
            dF4_ds = -k2*muscle_data.Lceopt(i);
            if (f>0) 										
                F4 = F4 + muscle_data.kSEE(i)*f^2;
                dF4_ds = dF4_ds - 2 * muscle_data.kSEE(i) * muscle_data.Lceopt(i) * f;
            end
        
        %Viscous Damping (F5)
            F5 = 0.001*sdot;
            dF5_dsdot = 0.001;

        %Compute Force Imbalance in Muscle Contraction 
            ff(2*i-1,1) = F4 - (a*F1*F2 + F3)*cosp-F5;   
            ff(2*i,1) = adot - (u(i)-a)*((u(i)/muscle_data.Tact(i)+(1-u(i))/muscle_data.Tdeact(i)));
            df_du(2*i,i) = (-(2*u(i)-a)/muscle_data.Tact(i)) - ((a-2*u(i) + 1)/muscle_data.Tdeact(i));
            df_dx(2*i-1,2*i-1) = dF4_ds - (a*(dF1_ds*F2 + F1*dF2_ds) + dF3_ds)*cosp - (a*F1*F2 + F3)*dcosp_ds;
            df_dx(2*i-1,2*i) = -F1*F2*cosp;
            df_dx(2*i,2*i) = (u(i)/muscle_data.Tact(i))+((1-u(i))/muscle_data.Tdeact(i));
            df_dxdot(2*i-1,2*i-1) = -a*F1*dF2_dsdot*cosp - dF5_dsdot;
            df_dxdot(2*i,2*i)=1;
            
        %Total Muscle Force
            force(i,:) = muscle_data.Fmax(i)*F4;     
    end
    %keyboard
end
    