classdef wenBayard_joint < handle
    %WENBAYARD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        prevError = [0;0]
        P = 10
        D = 5
    end
    
    methods
        function u = getU(obj, system)
            global timeStep stateVariables timeVector step
                  
            q1des = 0;
            q2des = -pi/6 + pi/2*sin(timeVector(step));
            dq1des = 0;
            dq2des = pi/2*cos(timeVector(step));
            ddq1des = 0;
            ddq2des = -pi/2*sin(timeVector(step));
  
            
            desiredState = vertcat(stateVariables(1,step), q1des, q2des, stateVariables(4,step), dq1des, dq2des);
            system.refreshMatricesForRK4(desiredState);
            
            ddq_des = [0; ddq1des; ddq2des;];
            usol = system.H \ (system.M * ddq_des + system.C);
            system.refreshMatrices();
            
            error = [q1des; q2des] - stateVariables(2:3, step);
            derror = (error-obj.prevError)/timeStep;
            
            u = usol + obj.P*error + obj.D*derror;
            obj.prevError = error;
            
        end
    end
    
end

