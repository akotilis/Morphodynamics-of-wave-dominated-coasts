function rmse = rootmeansquare_error(x_obs,x_model)

%Inputs: 
%x_obs --> Wave heights of observations (meters) 
%x_model --> Wave heights of model (meters) 
%For the calculation of the rmse, both x_obs and x_model must be the same
%length. If not the same length, the function will send a message. 

%Output: 
%rmse --> Root mean square error 

%The if checks if both vectors are the same length 
if length(x_obs) == length(x_model)
   sum = 0;  
    for ii = 1:length(x_obs) 
        
        sum_modelobs = sum + power(x_model(ii) - x_obs(ii),2);
        
        sum = sum_modelobs;    
        
    end
    
coeff = 1/length(x_obs); 

rms_error = sqrt(coeff*sum);

rmse = rms_error; 
else
    disp('The observations and the model do not have the same length! To use this function please check they are the same length')
    
    rmse = NaN; 
           
end

end