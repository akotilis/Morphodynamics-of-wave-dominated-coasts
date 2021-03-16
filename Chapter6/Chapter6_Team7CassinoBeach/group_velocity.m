function cg = group_velocity(c,n)

%Group velocity defined as Cg = c*n 
%Inputs: 
%c --> phase velocity [m/s]
%n --> propagation factor 
%Outputs:  
%Cg --> Group velocity [m/s] 

group_vel = n.*c;

cg = group_vel; 


end

