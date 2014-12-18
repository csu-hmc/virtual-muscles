function u=controller(t,nMus)

%----------------------
%Sinusoid Input
%----------------------
      u=zeros(nMus,1);
      u([1 4 5]) = 0.5 + 0.5*sin(2*pi*t);
      u([2 3 6]) = 0.5 + 0.5*sin(2*pi*t+pi);

%----------------------
%Alternative Step Input
%----------------------
%     tfinal=28.9980;       
%     u=zeros(nMus,1);
%     increment=tfinal/(nMus+3);
%     for i=1:nMus
%         start_time=i*increment;
%             if t>start_time
%                 u(i)=0.25; 
%             end
%     end
end
