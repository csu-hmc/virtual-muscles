function [Lm_poly1,Lm]=moment_arm(q)

    global muscle_data
%----------------------
%Obtain Muscle Lengths
%----------------------
    %Polynomial Coefficients from Muscle Data
        Lm_poly=[muscle_data.L0 muscle_data.dRhip muscle_data.dRknee muscle_data.dRankle];
        Lm_poly1=Lm_poly(:,3:4);
    %Apply Walking Data to Polynomial
        Lm=Lm_poly(:,1)- Lm_poly(:,2)*q(1)-Lm_poly(:,3)*q(2)-Lm_poly(:,4)*q(3);
end