function phase_vel = phase_velocity(L, T) 
%Calculates the phase velocity as Vp = Lambda/T 
%Inputs: 
%L --> Wavelength [meters]
%T --> Period [seconds]
%Outputs: 
%phase_vel --> phase velocity [m/s]

Vp = L./T; 

phase = Vp; 

end

