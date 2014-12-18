%-------------------------------------------------------------------------
%Simulator
%   simulates muscle dynamics
%-------------------------------------------------------------------------
    close all
    clear all
    clc
%------------------------------------------------
%Loading Muscle Properties and Parsing Variables
%------------------------------------------------
    %Import Data
        filename_muscle='muscle_data.txt';
        data = importdata(['Data' filesep filename_muscle]);
        muscle_names = {'R.Hamstrings','R.Rectus','R.Vasti','R.Gastroc',...
                        'R.Soleus','R.TibialisAnt'};
        Nmus = length(muscle_names);
    %Define Variables
        global muscle_data 
        HillA=0.25;
        muscle_data.PennOpt = zeros(6,1);
        muscle_data.gmax=repmat(1.5,6,1);
    %Obtain Remaining Variables
        muscle_data_new=zeros(Nmus,size(data.data,2));
        properties=data.textdata(2,2:end);
        for i=1:length(muscle_names)
            row=find(ismember(data.textdata(:,1),muscle_names{i}));
            muscle_data_new(i,:)=data.data(row-2,:);
            for j=1:length(properties);
                muscle_data.(char(properties(j)))=muscle_data_new(:,j);
            end
            muscle_data.c(i)=(muscle_data.Vmax(i)*HillA*(muscle_data.gmax(i)-1))/(HillA+1);
            muscle_data.kSEE(i) = 1/(muscle_data.umax(i)^2*muscle_data.SEEslack(i)^2);
        end
%---------------------
%Predict Joint Torque
%---------------------
    %Load Walking Data and Defining Joint Angles
        filename_walk='joint_angles.txt';
        data_walk = importdata(['Data' filesep filename_walk]);
        data_walk(:,2) = -data_walk(:,2);
        data_walk(:,3) = -(data_walk(:,3) + pi/2);
    %Filtering and Truncating Transient
        [num,den]=butter(2,6/(100/2));
        data_walk=filter(num,den,data_walk);
        data_walk=data_walk(100:end-100,:);
    %Defining Time
        frame=1:1:length(data_walk);                                        %Number of Frames in Data
        time_data=(frame-1)*0.01;                                           %Time of Data
        tfinal=time_data(end);                                              %Stop Time of Data
    %Vary Step Size 
        h = [0.00002 logspace(log10(0.0001), log10(0.016), 10)];
        RMS_knee = zeros(length(h)-1,1);
        RMS_ankle = zeros(length(h)-1,1);
        simulation_time = zeros(length(h),1);
        for i=1:length(h)
            fprintf('Simulating for %1.6f Step Size...',h(i))
            time = 0:h(i):tfinal;                                           %Time with New Step Size 
            data_walk_new = interp1(time_data,data_walk,time);              %Interpolate Data to New Step Size
            data_walk_new = data_walk_new(isfinite(data_walk_new(:,1)),:);  %Remove NANs from Data
            time = time(1:size(data_walk_new,1));                           %Clip Time Vector to Match Data without NANs
            tfinal = time(end);                                             %New End Time
            counter=1;
            tic
            %Start Simulation
                M=zeros(length(time),2);
                F=zeros(length(time),6);
                X=zeros(length(time),12);                                                   
                %Initial Conditions 
                    t=0;  
                    x1=repmat(2,Nmus,1);
                    x2=zeros(Nmus,1);             
                    x=[x1 x2];
                    x=reshape(x',Nmus*2,1);
                while t<tfinal + 1e-6
                    u=controller(t,Nmus);
                    q=data_walk_new(counter,:);
                    %q=zeros(1,3);
                    [x,Mnew,Fnew]=timestep(x,u,h(i),q);
                    t=t+h(i);
                    M(counter,:) = Mnew;
                    F(counter,:)=Fnew;
                    X(counter,:)=x';
                    counter=counter+1;
                end
            simulation_time(i,:)=toc;
            toc
        %RMS of Joint Torques
            if i <2
                %Store Values at h(1)
                    Knee_Moment0 = M(:,1);                              
                    Ankle_Moment0 = M(:,2);
                    tfinal0=tfinal;                                     
                    time0=time;  
                    peak_moment = max(abs(M));
            else 
                %Store Moments
                    Knee_Moment0_new = Knee_Moment0;
                    Ankle_Moment0_new = Ankle_Moment0;
                %Interpolate and Remove NANs
                    Knee_Moment = interp1(time,M(:,1),time0);           
                    Knee_Moment = Knee_Moment';
                    Ankle_Moment = interp1(time,M(:,2),time0);
                    Ankle_Moment = Ankle_Moment';
                    Knee_Nans = find(isnan(Knee_Moment));               
                    Ankle_Nans = find(isnan(Ankle_Moment));
                    Knee_Moment0_new(Knee_Nans)=[]; 
                    Ankle_Moment0_new(Ankle_Nans)=[];             
                    Knee_Moment(Knee_Nans)=[];
                    Ankle_Moment(Ankle_Nans)=[];
                    RMS_knee(i-1,:) = ((sqrt(mean((Knee_Moment0_new - Knee_Moment).^2)))/peak_moment(1))*100;
                    RMS_ankle(i-1,:) = ((sqrt(mean((Ankle_Moment0_new - Ankle_Moment).^2)))/peak_moment(2))*100;
            end
        end
%-------------------------------------------------------------------------
%Plot Results
%-------------------------------------------------------------------------
    figure(1)
        subplot(1,2,1)
        loglog(h(:,3:end),RMS_knee(2:end,:),'*-'); xlabel('Step Size','Fontweight','bold','Fontsize',16);
        ylabel('RMS Error (% Nm of Peak)','Fontweight','bold','Fontsize',16);
        title('Knee Torque','Fontweight','bold','Fontsize',18)
        subplot(1,2,2)
        loglog(h(:,3:end),RMS_ankle(2:end,:),'*-'); xlabel('Step Size','Fontweight','bold','Fontsize',16)
        ylabel('RMS Error (% Nm of Peak)','Fontweight','bold','Fontsize',16);
         title('Ankle Torque','Fontweight','bold','Fontsize',18)
        
    figure(2)
        subplot(2,1,1)
        plot(time, M,'Linewidth',2); xlabel('Time (s)','Fontweight','bold','Fontsize',20);
        ylabel('Torque (Nm)','Fontweight','bold','Fontsize',20);
        title('Predicted Torques','Fontweight','bold','Fontsize',22);
        legend('Knee','Ankle'); %xlim([15 18]); 
        set(gca,'FontSize',16)
        subplot(2,1,2)
        plot(time,F,'Linewidth',2); xlabel('Time (s)','Fontweight','bold','Fontsize',20)
        ylabel('Muscle Force (N)','Fontweight','bold','Fontsize',20)
        title('Predicted Muscle Forces','Fontweight','bold','Fontsize',22);
        legend(muscle_names); %xlim([15 18]); 
        set(gca,'FontSize',16)