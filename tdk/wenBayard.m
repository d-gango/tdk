classdef wenBayard < handle
    %WENBAYARD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function u = getU(obj, system)
            global timeStep stateVariables timeVector maxStep step
            m0 = system.m0;
            m1 = system.m1;
            m2 = system.m2;
            l0 = system.l0;
            l1 = system.l1;
            l2 = system.l2;
            I0 = system.I0;
            I1 = system.I1;
            I2 = system.I2;
            fi0 = system.fi0;
            
           vx = -0.05*timeVector(step)*2*sin(2*timeVector(step));
           vy = 0.05*timeVector(step)*2*cos(2*timeVector(step));
            k = system.angularMomentum();
            
            if step == 1
                q_des = stateVariables(1:3, 1);
                dq_des = stateVariables(4:6, 1);
                ddq_des = [0;0;0];
            else
                q_des = stateVariables(1:3, step-1) + timeStep * stateVariables(4:6, step-1);
                q0des = q_des(1);
                q1des = q_des(2);
                q2des = q_des(3);
                dq_des = [-((m0 + m1 + m2)*(8*I1*l2*m0^2*vx*cos(q0des + q1des + q2des) - 4*I2*l1*m1^2*vy*sin(q0des + q1des) - 8*I2*l1*m0^2*vy*sin(q0des + q1des) + 8*I1*l2*m1^2*vx*cos(q0des + q1des + q2des) + 4*I1*l2*m2^2*vx*cos(q0des + q1des + q2des) + 8*I1*l2*m0^2*vy*sin(q0des + q1des + q2des) + 8*I1*l2*m1^2*vy*sin(q0des + q1des + q2des) + 4*I1*l2*m2^2*vy*sin(q0des + q1des + q2des) - 8*I2*l1*m0^2*vx*cos(q0des + q1des) - 4*I2*l1*m1^2*vx*cos(q0des + q1des) - 12*I2*l1*m0*m1*vx*cos(q0des + q1des) - 8*I2*l1*m0*m2*vx*cos(q0des + q1des) - 4*I2*l1*m1*m2*vx*cos(q0des + q1des) - 12*I2*l1*m0*m1*vy*sin(q0des + q1des) - 8*I2*l1*m0*m2*vy*sin(q0des + q1des) - 4*I2*l1*m1*m2*vy*sin(q0des + q1des) + 8*k*l1*l2*m0^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) - 8*k*l1*l2*m0^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 4*k*l1*l2*m1^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) - 4*k*l1*l2*m1^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 2*l1^2*l2*m0*m1^2*vx*cos(q0des + q1des + q2des) + 2*l1^2*l2*m0^2*m1*vx*cos(q0des + q1des + q2des) + 4*l1^2*l2*m0*m2^2*vx*cos(q0des + q1des + q2des) + 8*l1^2*l2*m0^2*m2*vx*cos(q0des + q1des + q2des) + l1^2*l2*m1*m2^2*vx*cos(q0des + q1des + q2des) + 2*l1^2*l2*m1^2*m2*vx*cos(q0des + q1des + q2des) + 2*l1^2*l2*m0*m1^2*vy*sin(q0des + q1des + q2des) + 2*l1^2*l2*m0^2*m1*vy*sin(q0des + q1des + q2des) + 4*l1^2*l2*m0*m2^2*vy*sin(q0des + q1des + q2des) + 8*l1^2*l2*m0^2*m2*vy*sin(q0des + q1des + q2des) + l1^2*l2*m1*m2^2*vy*sin(q0des + q1des + q2des) + 2*l1^2*l2*m1^2*m2*vy*sin(q0des + q1des + q2des) + 16*I1*l2*m0*m1*vx*cos(q0des + q1des + q2des) + 12*I1*l2*m0*m2*vx*cos(q0des + q1des + q2des) + 12*I1*l2*m1*m2*vx*cos(q0des + q1des + q2des) + 16*I1*l2*m0*m1*vy*sin(q0des + q1des + q2des) + 12*I1*l2*m0*m2*vy*sin(q0des + q1des + q2des) + 12*I1*l2*m1*m2*vy*sin(q0des + q1des + q2des) - 2*l1*l2^2*m0^2*m2*vx*cos(q0des + q1des) - l1*l2^2*m1^2*m2*vx*cos(q0des + q1des) - 2*l1*l2^2*m0^2*m2*vy*sin(q0des + q1des) - l1*l2^2*m1^2*m2*vy*sin(q0des + q1des) - 3*l1*l2^2*m0*m1*m2*vy*sin(q0des + q1des) + 2*l1*l2^2*m0*m2^2*vx*cos(q0des + q1des + q2des)*cos(q2des) + 4*l1*l2^2*m0^2*m2*vx*cos(q0des + q1des + q2des)*cos(q2des) + l1*l2^2*m1*m2^2*vx*cos(q0des + q1des + q2des)*cos(q2des) + 2*l1*l2^2*m1^2*m2*vx*cos(q0des + q1des + q2des)*cos(q2des) + 2*l1*l2^2*m0*m2^2*vy*sin(q0des + q1des + q2des)*cos(q2des) + 4*l1*l2^2*m0^2*m2*vy*sin(q0des + q1des + q2des)*cos(q2des) + l1*l2^2*m1*m2^2*vy*sin(q0des + q1des + q2des)*cos(q2des) + 2*l1*l2^2*m1^2*m2*vy*sin(q0des + q1des + q2des)*cos(q2des) + 12*k*l1*l2*m0*m1*cos(q0des + q1des + q2des)*sin(q0des + q1des) - 12*k*l1*l2*m0*m1*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 4*k*l1*l2*m0*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) - 4*k*l1*l2*m0*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 2*k*l1*l2*m1*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) - 2*k*l1*l2*m1*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 11*l1^2*l2*m0*m1*m2*vx*cos(q0des + q1des + q2des) + 11*l1^2*l2*m0*m1*m2*vy*sin(q0des + q1des + q2des) - 4*l1^2*l2*m0^2*m2*vx*cos(q0des + q1des)*cos(q2des) - l1^2*l2*m1^2*m2*vx*cos(q0des + q1des)*cos(q2des) - 4*l1^2*l2*m0^2*m2*vy*sin(q0des + q1des)*cos(q2des) - l1^2*l2*m1^2*m2*vy*sin(q0des + q1des)*cos(q2des) - 3*l1*l2^2*m0*m1*m2*vx*cos(q0des + q1des) + 6*l1*l2^2*m0*m1*m2*vx*cos(q0des + q1des + q2des)*cos(q2des) + 6*l1*l2^2*m0*m1*m2*vy*sin(q0des + q1des + q2des)*cos(q2des) + 4*l0*l1*l2*m0*m1^2*vx*cos(q0des + q1des + q2des)*cos(fi0 - q1des) + 4*l0*l1*l2*m0^2*m1*vx*cos(q0des + q1des + q2des)*cos(fi0 - q1des) + 4*l0*l1*l2*m0*m2^2*vx*cos(q0des + q1des + q2des)*cos(fi0 - q1des) + 8*l0*l1*l2*m0^2*m2*vx*cos(q0des + q1des + q2des)*cos(fi0 - q1des) - 4*l0*l1*l2*m0^2*m2*vx*cos(q1des - fi0 + q2des)*cos(q0des + q1des) + 4*l0*l1*l2*m0*m1^2*vy*sin(q0des + q1des + q2des)*cos(fi0 - q1des) + 4*l0*l1*l2*m0^2*m1*vy*sin(q0des + q1des + q2des)*cos(fi0 - q1des) + 4*l0*l1*l2*m0*m2^2*vy*sin(q0des + q1des + q2des)*cos(fi0 - q1des) + 8*l0*l1*l2*m0^2*m2*vy*sin(q0des + q1des + q2des)*cos(fi0 - q1des) - 4*l0*l1*l2*m0^2*m2*vy*cos(q1des - fi0 + q2des)*sin(q0des + q1des) - 4*l1^2*l2*m0*m1*m2*vx*cos(q0des + q1des)*cos(q2des) - 4*l1^2*l2*m0*m1*m2*vy*sin(q0des + q1des)*cos(q2des) + 10*l0*l1*l2*m0*m1*m2*vx*cos(q0des + q1des + q2des)*cos(fi0 - q1des) - 2*l0*l1*l2*m0*m1*m2*vx*cos(q1des - fi0 + q2des)*cos(q0des + q1des) + 10*l0*l1*l2*m0*m1*m2*vy*sin(q0des + q1des + q2des)*cos(fi0 - q1des) - 2*l0*l1*l2*m0*m1*m2*vy*cos(q1des - fi0 + q2des)*sin(q0des + q1des)))/(8*I1*l0*l2*m0^3*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 8*I1*l0*l2*m0^3*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 8*I0*l1*l2*m0^3*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 8*I0*l1*l2*m0^3*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 4*I0*l1*l2*m1^3*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 4*I0*l1*l2*m1^3*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 8*I2*l0*l1*m0^3*cos(fi0 + q0des)*sin(q0des + q1des) - 8*I2*l0*l1*m0^3*cos(q0des + q1des)*sin(fi0 + q0des) - 4*l0^2*l1*l2*m0*m1^3*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 4*l0^2*l1*l2*m0*m1^3*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 8*l0^2*l1*l2*m0^3*m1*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 8*l0^2*l1*l2*m0^3*m1*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 8*l0^2*l1*l2*m0^3*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 8*l0^2*l1*l2*m0^3*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 4*I2*l0*l1*m0*m1^2*cos(fi0 + q0des)*sin(q0des + q1des) - 4*I2*l0*l1*m0*m1^2*cos(q0des + q1des)*sin(fi0 + q0des) + 12*I2*l0*l1*m0^2*m1*cos(fi0 + q0des)*sin(q0des + q1des) - 12*I2*l0*l1*m0^2*m1*cos(q0des + q1des)*sin(fi0 + q0des) + 8*I2*l0*l1*m0^2*m2*cos(fi0 + q0des)*sin(q0des + q1des) - 8*I2*l0*l1*m0^2*m2*cos(q0des + q1des)*sin(fi0 + q0des) + 2*l0*l1^2*l2*m0^2*m1^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 2*l0*l1^2*l2*m0^2*m1^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 4*l0*l1^2*l2*m0^2*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 4*l0*l1^2*l2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 12*l0^2*l1*l2*m0^2*m1^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 12*l0^2*l1*l2*m0^2*m1^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 4*l0^2*l1*l2*m0^2*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 4*l0^2*l1*l2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 8*I1*l0*l2*m0*m1^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 8*I1*l0*l2*m0*m1^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 16*I1*l0*l2*m0^2*m1*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 16*I1*l0*l2*m0^2*m1*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 4*I1*l0*l2*m0*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 4*I1*l0*l2*m0*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 12*I1*l0*l2*m0^2*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 12*I1*l0*l2*m0^2*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 16*I0*l1*l2*m0*m1^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 16*I0*l1*l2*m0*m1^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 20*I0*l1*l2*m0^2*m1*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 20*I0*l1*l2*m0^2*m1*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 4*I0*l1*l2*m0*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 4*I0*l1*l2*m0*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 12*I0*l1*l2*m0^2*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 12*I0*l1*l2*m0^2*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 2*I0*l1*l2*m1*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 2*I0*l1*l2*m1*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 6*I0*l1*l2*m1^2*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 6*I0*l1*l2*m1^2*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 2*l0*l1*l2^2*m0^3*m2*cos(fi0 + q0des)*sin(q0des + q1des) - 2*l0*l1*l2^2*m0^3*m2*cos(q0des + q1des)*sin(fi0 + q0des) + 2*l0*l1^2*l2*m0^3*m1*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 2*l0*l1^2*l2*m0^3*m1*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 8*l0*l1^2*l2*m0^3*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 8*l0*l1^2*l2*m0^3*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 2*l0*l1^2*l2*m0*m1^3*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 2*l0*l1^2*l2*m0*m1^3*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) - 4*l0*l1^2*l2*m0^3*m1*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 4*l0*l1^2*l2*m0^3*m1*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) - 8*l0*l1^2*l2*m0^3*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 8*l0*l1^2*l2*m0^3*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) + 4*I2*l0*l1*m0*m1*m2*cos(fi0 + q0des)*sin(q0des + q1des) - 4*I2*l0*l1*m0*m1*m2*cos(q0des + q1des)*sin(fi0 + q0des) + 2*l0*l1*l2^2*m0^2*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(q2des) - 2*l0*l1*l2^2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(q2des) + 4*l0^2*l1*l2*m0^2*m1^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(fi0 - q1des) - 4*l0^2*l1*l2*m0^2*m1^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(fi0 - q1des) + 4*l0^2*l1*l2*m0^2*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(fi0 - q1des) - 4*l0^2*l1*l2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(fi0 - q1des) + 4*l0*l1^2*l2*m0^3*m2*cos(fi0 + q0des)*sin(q0des + q1des)*cos(q2des) - 4*l0*l1^2*l2*m0^3*m2*cos(q0des + q1des)*sin(fi0 + q0des)*cos(q2des) - 6*l0*l1^2*l2*m0^2*m1^2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 6*l0*l1^2*l2*m0^2*m1^2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) - 4*l0*l1^2*l2*m0^2*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 4*l0*l1^2*l2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) + l0*l1^2*l2*m0*m1*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - l0*l1^2*l2*m0*m1*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 2*l0*l1^2*l2*m0*m1^2*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 2*l0*l1^2*l2*m0*m1^2*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 11*l0*l1^2*l2*m0^2*m1*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 11*l0*l1^2*l2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 4*l0*l1*l2^2*m0^3*m2*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*sin(q0des + q1des) + 4*l0*l1*l2^2*m0^3*m2*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*cos(q0des + q1des) - 2*l0^2*l1*l2*m0*m1*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 2*l0^2*l1*l2*m0*m1*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 6*l0^2*l1*l2*m0*m1^2*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 6*l0^2*l1*l2*m0*m1^2*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 16*l0^2*l1*l2*m0^2*m1*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 16*l0^2*l1*l2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 12*I1*l0*l2*m0*m1*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 12*I1*l0*l2*m0*m1*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 18*I0*l1*l2*m0*m1*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 18*I0*l1*l2*m0*m1*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 4*l0*l1*l2^2*m0^3*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(q2des) - 4*l0*l1*l2^2*m0^3*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(q2des) - 2*l0*l1*l2^2*m0^2*m2^2*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*sin(q0des + q1des) + 2*l0*l1*l2^2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*cos(q0des + q1des) + 4*l0^2*l1*l2*m0^3*m1*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(fi0 - q1des) - 4*l0^2*l1*l2*m0^3*m1*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(fi0 - q1des) + 8*l0^2*l1*l2*m0^3*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(fi0 - q1des) - 8*l0^2*l1*l2*m0^3*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(fi0 - q1des) + 4*l0^2*l1*l2*m0^3*m2*cos(q1des - fi0 + q2des)*cos(fi0 + q0des)*sin(q0des + q1des) - 4*l0^2*l1*l2*m0^3*m2*cos(q1des - fi0 + q2des)*cos(q0des + q1des)*sin(fi0 + q0des) + l0*l1*l2^2*m0*m1^2*m2*cos(fi0 + q0des)*sin(q0des + q1des) - l0*l1*l2^2*m0*m1^2*m2*cos(q0des + q1des)*sin(fi0 + q0des) + 3*l0*l1*l2^2*m0^2*m1*m2*cos(fi0 + q0des)*sin(q0des + q1des) - 3*l0*l1*l2^2*m0^2*m1*m2*cos(q0des + q1des)*sin(fi0 + q0des) + l0*l1*l2^2*m0*m1*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(q2des) - l0*l1*l2^2*m0*m1*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(q2des) + 2*l0*l1*l2^2*m0*m1^2*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(q2des) - 2*l0*l1*l2^2*m0*m1^2*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(q2des) + 6*l0*l1*l2^2*m0^2*m1*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(q2des) - 6*l0*l1*l2^2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(q2des) + 10*l0^2*l1*l2*m0^2*m1*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(fi0 - q1des) - 10*l0^2*l1*l2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(fi0 - q1des) + 2*l0^2*l1*l2*m0^2*m1*m2*cos(q1des - fi0 + q2des)*cos(fi0 + q0des)*sin(q0des + q1des) - 2*l0^2*l1*l2*m0^2*m1*m2*cos(q1des - fi0 + q2des)*cos(q0des + q1des)*sin(fi0 + q0des) - 2*l0*l1^2*l2*m0*m1*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 2*l0*l1^2*l2*m0*m1*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) - 5*l0*l1^2*l2*m0*m1^2*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 5*l0*l1^2*l2*m0*m1^2*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) - 14*l0*l1^2*l2*m0^2*m1*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 14*l0*l1^2*l2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) + l0*l1^2*l2*m0*m1^2*m2*cos(fi0 + q0des)*sin(q0des + q1des)*cos(q2des) - l0*l1^2*l2*m0*m1^2*m2*cos(q0des + q1des)*sin(fi0 + q0des)*cos(q2des) + 4*l0*l1^2*l2*m0^2*m1*m2*cos(fi0 + q0des)*sin(q0des + q1des)*cos(q2des) - 4*l0*l1^2*l2*m0^2*m1*m2*cos(q0des + q1des)*sin(fi0 + q0des)*cos(q2des) - l0*l1*l2^2*m0*m1*m2^2*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*sin(q0des + q1des) + l0*l1*l2^2*m0*m1*m2^2*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*cos(q0des + q1des) - 2*l0*l1*l2^2*m0*m1^2*m2*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*sin(q0des + q1des) + 2*l0*l1*l2^2*m0*m1^2*m2*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*cos(q0des + q1des) - 6*l0*l1*l2^2*m0^2*m1*m2*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*sin(q0des + q1des) + 6*l0*l1*l2^2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*cos(q0des + q1des));
                          ((m0 + m1 + m2)*(8*I0*l2*m0^2*vx*cos(q0des + q1des + q2des) - 8*I2*l1*m0^2*vy*sin(q0des + q1des) - 4*I2*l1*m1^2*vy*sin(q0des + q1des) - 8*I2*l0*m0^2*vy*sin(fi0 + q0des) + 8*I0*l2*m1^2*vx*cos(q0des + q1des + q2des) + 8*I1*l2*m0^2*vx*cos(q0des + q1des + q2des) + 4*I0*l2*m2^2*vx*cos(q0des + q1des + q2des) + 8*I1*l2*m1^2*vx*cos(q0des + q1des + q2des) + 4*I1*l2*m2^2*vx*cos(q0des + q1des + q2des) + 8*I0*l2*m0^2*vy*sin(q0des + q1des + q2des) + 8*I0*l2*m1^2*vy*sin(q0des + q1des + q2des) + 8*I1*l2*m0^2*vy*sin(q0des + q1des + q2des) + 4*I0*l2*m2^2*vy*sin(q0des + q1des + q2des) + 8*I1*l2*m1^2*vy*sin(q0des + q1des + q2des) + 4*I1*l2*m2^2*vy*sin(q0des + q1des + q2des) - 8*I2*l0*m0^2*vx*cos(fi0 + q0des) - 8*I2*l1*m0^2*vx*cos(q0des + q1des) - 4*I2*l1*m1^2*vx*cos(q0des + q1des) - 12*I2*l1*m0*m1*vx*cos(q0des + q1des) - 8*I2*l1*m0*m2*vx*cos(q0des + q1des) - 4*I2*l1*m1*m2*vx*cos(q0des + q1des) - 8*I2*l0*m0*m1*vy*sin(fi0 + q0des) - 8*I2*l0*m0*m2*vy*sin(fi0 + q0des) - 12*I2*l1*m0*m1*vy*sin(q0des + q1des) - 8*I2*l1*m0*m2*vy*sin(q0des + q1des) - 4*I2*l1*m1*m2*vy*sin(q0des + q1des) + 8*k*l0*l2*m0^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 8*k*l0*l2*m0^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 8*k*l1*l2*m0^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) - 8*k*l1*l2*m0^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 4*k*l1*l2*m1^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) - 4*k*l1*l2*m1^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 8*l0^2*l2*m0*m1^2*vx*cos(q0des + q1des + q2des) + 8*l0^2*l2*m0^2*m1*vx*cos(q0des + q1des + q2des) + 4*l0^2*l2*m0*m2^2*vx*cos(q0des + q1des + q2des) + 8*l0^2*l2*m0^2*m2*vx*cos(q0des + q1des + q2des) + 2*l1^2*l2*m0*m1^2*vx*cos(q0des + q1des + q2des) + 2*l1^2*l2*m0^2*m1*vx*cos(q0des + q1des + q2des) + 4*l1^2*l2*m0*m2^2*vx*cos(q0des + q1des + q2des) + 8*l1^2*l2*m0^2*m2*vx*cos(q0des + q1des + q2des) + l1^2*l2*m1*m2^2*vx*cos(q0des + q1des + q2des) + 2*l1^2*l2*m1^2*m2*vx*cos(q0des + q1des + q2des) + 8*l0^2*l2*m0*m1^2*vy*sin(q0des + q1des + q2des) + 8*l0^2*l2*m0^2*m1*vy*sin(q0des + q1des + q2des) + 4*l0^2*l2*m0*m2^2*vy*sin(q0des + q1des + q2des) + 8*l0^2*l2*m0^2*m2*vy*sin(q0des + q1des + q2des) + 2*l1^2*l2*m0*m1^2*vy*sin(q0des + q1des + q2des) + 2*l1^2*l2*m0^2*m1*vy*sin(q0des + q1des + q2des) + 4*l1^2*l2*m0*m2^2*vy*sin(q0des + q1des + q2des) + 8*l1^2*l2*m0^2*m2*vy*sin(q0des + q1des + q2des) + l1^2*l2*m1*m2^2*vy*sin(q0des + q1des + q2des) + 2*l1^2*l2*m1^2*m2*vy*sin(q0des + q1des + q2des) + 16*I0*l2*m0*m1*vx*cos(q0des + q1des + q2des) + 12*I0*l2*m0*m2*vx*cos(q0des + q1des + q2des) + 16*I1*l2*m0*m1*vx*cos(q0des + q1des + q2des) + 12*I0*l2*m1*m2*vx*cos(q0des + q1des + q2des) + 12*I1*l2*m0*m2*vx*cos(q0des + q1des + q2des) + 12*I1*l2*m1*m2*vx*cos(q0des + q1des + q2des) + 16*I0*l2*m0*m1*vy*sin(q0des + q1des + q2des) + 12*I0*l2*m0*m2*vy*sin(q0des + q1des + q2des) + 16*I1*l2*m0*m1*vy*sin(q0des + q1des + q2des) + 12*I0*l2*m1*m2*vy*sin(q0des + q1des + q2des) + 12*I1*l2*m0*m2*vy*sin(q0des + q1des + q2des) + 12*I1*l2*m1*m2*vy*sin(q0des + q1des + q2des) - 2*l0*l2^2*m0^2*m2*vx*cos(fi0 + q0des) - 2*l1*l2^2*m0^2*m2*vx*cos(q0des + q1des) - l1*l2^2*m1^2*m2*vx*cos(q0des + q1des) - 2*l0*l2^2*m0^2*m2*vy*sin(fi0 + q0des) - 2*l1*l2^2*m0^2*m2*vy*sin(q0des + q1des) - l1*l2^2*m1^2*m2*vy*sin(q0des + q1des) - 8*I2*l0*m0*m1*vx*cos(fi0 + q0des) - 8*I2*l0*m0*m2*vx*cos(fi0 + q0des) - 3*l1*l2^2*m0*m1*m2*vy*sin(q0des + q1des) + 2*l1*l2^2*m0*m2^2*vx*cos(q0des + q1des + q2des)*cos(q2des) + 4*l1*l2^2*m0^2*m2*vx*cos(q0des + q1des + q2des)*cos(q2des) + l1*l2^2*m1*m2^2*vx*cos(q0des + q1des + q2des)*cos(q2des) + 2*l1*l2^2*m1^2*m2*vx*cos(q0des + q1des + q2des)*cos(q2des) + 8*k*l0*l2*m0*m1*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 8*k*l0*l2*m0*m1*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 4*k*l0*l2*m0*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 4*k*l0*l2*m0*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 2*l1*l2^2*m0*m2^2*vy*sin(q0des + q1des + q2des)*cos(q2des) + 4*l1*l2^2*m0^2*m2*vy*sin(q0des + q1des + q2des)*cos(q2des) + l1*l2^2*m1*m2^2*vy*sin(q0des + q1des + q2des)*cos(q2des) + 2*l1*l2^2*m1^2*m2*vy*sin(q0des + q1des + q2des)*cos(q2des) + 12*k*l1*l2*m0*m1*cos(q0des + q1des + q2des)*sin(q0des + q1des) - 12*k*l1*l2*m0*m1*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 4*k*l1*l2*m0*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) - 4*k*l1*l2*m0*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 2*k*l1*l2*m1*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) - 2*k*l1*l2*m1*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 4*l0^2*l2*m0^2*m2*vx*cos(q1des - fi0 + q2des)*cos(fi0 + q0des) - 4*l0^2*l2*m0^2*m2*vy*cos(q1des - fi0 + q2des)*sin(fi0 + q0des) + 12*l0^2*l2*m0*m1*m2*vx*cos(q0des + q1des + q2des) + 11*l1^2*l2*m0*m1*m2*vx*cos(q0des + q1des + q2des) + 12*l0^2*l2*m0*m1*m2*vy*sin(q0des + q1des + q2des) + 11*l1^2*l2*m0*m1*m2*vy*sin(q0des + q1des + q2des) - 4*l1^2*l2*m0^2*m2*vx*cos(q0des + q1des)*cos(q2des) - l1^2*l2*m1^2*m2*vx*cos(q0des + q1des)*cos(q2des) + 2*l0*l2^2*m0*m2^2*vx*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des) + 4*l0*l2^2*m0^2*m2*vx*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des) - 4*l1^2*l2*m0^2*m2*vy*sin(q0des + q1des)*cos(q2des) - l1^2*l2*m1^2*m2*vy*sin(q0des + q1des)*cos(q2des) + 2*l0*l2^2*m0*m2^2*vy*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des) + 4*l0*l2^2*m0^2*m2*vy*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des) - 2*l0*l2^2*m0*m1*m2*vx*cos(fi0 + q0des) - 3*l1*l2^2*m0*m1*m2*vx*cos(q0des + q1des) - 2*l0*l2^2*m0*m1*m2*vy*sin(fi0 + q0des) + 6*l1*l2^2*m0*m1*m2*vx*cos(q0des + q1des + q2des)*cos(q2des) + 6*l1*l2^2*m0*m1*m2*vy*sin(q0des + q1des + q2des)*cos(q2des) + 8*l0*l1*l2*m0*m1^2*vx*cos(q0des + q1des + q2des)*cos(fi0 - q1des) + 8*l0*l1*l2*m0^2*m1*vx*cos(q0des + q1des + q2des)*cos(fi0 - q1des) + 8*l0*l1*l2*m0*m2^2*vx*cos(q0des + q1des + q2des)*cos(fi0 - q1des) + 16*l0*l1*l2*m0^2*m2*vx*cos(q0des + q1des + q2des)*cos(fi0 - q1des) - 4*l0*l1*l2*m0^2*m2*vx*cos(q1des - fi0 + q2des)*cos(q0des + q1des) + 8*l0*l1*l2*m0*m1^2*vy*sin(q0des + q1des + q2des)*cos(fi0 - q1des) + 8*l0*l1*l2*m0^2*m1*vy*sin(q0des + q1des + q2des)*cos(fi0 - q1des) + 8*l0*l1*l2*m0*m2^2*vy*sin(q0des + q1des + q2des)*cos(fi0 - q1des) + 16*l0*l1*l2*m0^2*m2*vy*sin(q0des + q1des + q2des)*cos(fi0 - q1des) - 4*l0*l1*l2*m0^2*m2*vy*cos(q1des - fi0 + q2des)*sin(q0des + q1des) - 4*l0*l1*l2*m0^2*m2*vx*cos(fi0 + q0des)*cos(q2des) - 4*l1^2*l2*m0*m1*m2*vx*cos(q0des + q1des)*cos(q2des) - 4*l0*l1*l2*m0^2*m2*vy*sin(fi0 + q0des)*cos(q2des) + 4*l0*l2^2*m0*m1*m2*vx*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des) - 4*l1^2*l2*m0*m1*m2*vy*sin(q0des + q1des)*cos(q2des) + 4*l0*l2^2*m0*m1*m2*vy*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des) + 20*l0*l1*l2*m0*m1*m2*vx*cos(q0des + q1des + q2des)*cos(fi0 - q1des) - 2*l0*l1*l2*m0*m1*m2*vx*cos(q1des - fi0 + q2des)*cos(q0des + q1des) + 20*l0*l1*l2*m0*m1*m2*vy*sin(q0des + q1des + q2des)*cos(fi0 - q1des) - 2*l0*l1*l2*m0*m1*m2*vy*cos(q1des - fi0 + q2des)*sin(q0des + q1des) - 2*l0*l1*l2*m0*m1*m2*vx*cos(fi0 + q0des)*cos(q2des) - 2*l0*l1*l2*m0*m1*m2*vy*sin(fi0 + q0des)*cos(q2des)))/(8*I1*l0*l2*m0^3*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 8*I1*l0*l2*m0^3*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 8*I0*l1*l2*m0^3*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 8*I0*l1*l2*m0^3*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 4*I0*l1*l2*m1^3*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 4*I0*l1*l2*m1^3*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 8*I2*l0*l1*m0^3*cos(fi0 + q0des)*sin(q0des + q1des) - 8*I2*l0*l1*m0^3*cos(q0des + q1des)*sin(fi0 + q0des) - 4*l0^2*l1*l2*m0*m1^3*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 4*l0^2*l1*l2*m0*m1^3*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 8*l0^2*l1*l2*m0^3*m1*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 8*l0^2*l1*l2*m0^3*m1*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 8*l0^2*l1*l2*m0^3*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 8*l0^2*l1*l2*m0^3*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 4*I2*l0*l1*m0*m1^2*cos(fi0 + q0des)*sin(q0des + q1des) - 4*I2*l0*l1*m0*m1^2*cos(q0des + q1des)*sin(fi0 + q0des) + 12*I2*l0*l1*m0^2*m1*cos(fi0 + q0des)*sin(q0des + q1des) - 12*I2*l0*l1*m0^2*m1*cos(q0des + q1des)*sin(fi0 + q0des) + 8*I2*l0*l1*m0^2*m2*cos(fi0 + q0des)*sin(q0des + q1des) - 8*I2*l0*l1*m0^2*m2*cos(q0des + q1des)*sin(fi0 + q0des) + 2*l0*l1^2*l2*m0^2*m1^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 2*l0*l1^2*l2*m0^2*m1^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 4*l0*l1^2*l2*m0^2*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 4*l0*l1^2*l2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 12*l0^2*l1*l2*m0^2*m1^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 12*l0^2*l1*l2*m0^2*m1^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 4*l0^2*l1*l2*m0^2*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 4*l0^2*l1*l2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 8*I1*l0*l2*m0*m1^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 8*I1*l0*l2*m0*m1^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 16*I1*l0*l2*m0^2*m1*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 16*I1*l0*l2*m0^2*m1*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 4*I1*l0*l2*m0*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 4*I1*l0*l2*m0*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 12*I1*l0*l2*m0^2*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 12*I1*l0*l2*m0^2*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 16*I0*l1*l2*m0*m1^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 16*I0*l1*l2*m0*m1^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 20*I0*l1*l2*m0^2*m1*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 20*I0*l1*l2*m0^2*m1*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 4*I0*l1*l2*m0*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 4*I0*l1*l2*m0*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 12*I0*l1*l2*m0^2*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 12*I0*l1*l2*m0^2*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 2*I0*l1*l2*m1*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 2*I0*l1*l2*m1*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 6*I0*l1*l2*m1^2*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 6*I0*l1*l2*m1^2*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 2*l0*l1*l2^2*m0^3*m2*cos(fi0 + q0des)*sin(q0des + q1des) - 2*l0*l1*l2^2*m0^3*m2*cos(q0des + q1des)*sin(fi0 + q0des) + 2*l0*l1^2*l2*m0^3*m1*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 2*l0*l1^2*l2*m0^3*m1*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 8*l0*l1^2*l2*m0^3*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 8*l0*l1^2*l2*m0^3*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 2*l0*l1^2*l2*m0*m1^3*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 2*l0*l1^2*l2*m0*m1^3*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) - 4*l0*l1^2*l2*m0^3*m1*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 4*l0*l1^2*l2*m0^3*m1*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) - 8*l0*l1^2*l2*m0^3*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 8*l0*l1^2*l2*m0^3*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) + 4*I2*l0*l1*m0*m1*m2*cos(fi0 + q0des)*sin(q0des + q1des) - 4*I2*l0*l1*m0*m1*m2*cos(q0des + q1des)*sin(fi0 + q0des) + 2*l0*l1*l2^2*m0^2*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(q2des) - 2*l0*l1*l2^2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(q2des) + 4*l0^2*l1*l2*m0^2*m1^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(fi0 - q1des) - 4*l0^2*l1*l2*m0^2*m1^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(fi0 - q1des) + 4*l0^2*l1*l2*m0^2*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(fi0 - q1des) - 4*l0^2*l1*l2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(fi0 - q1des) + 4*l0*l1^2*l2*m0^3*m2*cos(fi0 + q0des)*sin(q0des + q1des)*cos(q2des) - 4*l0*l1^2*l2*m0^3*m2*cos(q0des + q1des)*sin(fi0 + q0des)*cos(q2des) - 6*l0*l1^2*l2*m0^2*m1^2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 6*l0*l1^2*l2*m0^2*m1^2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) - 4*l0*l1^2*l2*m0^2*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 4*l0*l1^2*l2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) + l0*l1^2*l2*m0*m1*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - l0*l1^2*l2*m0*m1*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 2*l0*l1^2*l2*m0*m1^2*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 2*l0*l1^2*l2*m0*m1^2*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) + 11*l0*l1^2*l2*m0^2*m1*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 11*l0*l1^2*l2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 4*l0*l1*l2^2*m0^3*m2*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*sin(q0des + q1des) + 4*l0*l1*l2^2*m0^3*m2*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*cos(q0des + q1des) - 2*l0^2*l1*l2*m0*m1*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 2*l0^2*l1*l2*m0*m1*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 6*l0^2*l1*l2*m0*m1^2*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 6*l0^2*l1*l2*m0*m1^2*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) - 16*l0^2*l1*l2*m0^2*m1*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 16*l0^2*l1*l2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 12*I1*l0*l2*m0*m1*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des) - 12*I1*l0*l2*m0*m1*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des) - 18*I0*l1*l2*m0*m1*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des) + 18*I0*l1*l2*m0*m1*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des) + 4*l0*l1*l2^2*m0^3*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(q2des) - 4*l0*l1*l2^2*m0^3*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(q2des) - 2*l0*l1*l2^2*m0^2*m2^2*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*sin(q0des + q1des) + 2*l0*l1*l2^2*m0^2*m2^2*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*cos(q0des + q1des) + 4*l0^2*l1*l2*m0^3*m1*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(fi0 - q1des) - 4*l0^2*l1*l2*m0^3*m1*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(fi0 - q1des) + 8*l0^2*l1*l2*m0^3*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(fi0 - q1des) - 8*l0^2*l1*l2*m0^3*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(fi0 - q1des) + 4*l0^2*l1*l2*m0^3*m2*cos(q1des - fi0 + q2des)*cos(fi0 + q0des)*sin(q0des + q1des) - 4*l0^2*l1*l2*m0^3*m2*cos(q1des - fi0 + q2des)*cos(q0des + q1des)*sin(fi0 + q0des) + l0*l1*l2^2*m0*m1^2*m2*cos(fi0 + q0des)*sin(q0des + q1des) - l0*l1*l2^2*m0*m1^2*m2*cos(q0des + q1des)*sin(fi0 + q0des) + 3*l0*l1*l2^2*m0^2*m1*m2*cos(fi0 + q0des)*sin(q0des + q1des) - 3*l0*l1*l2^2*m0^2*m1*m2*cos(q0des + q1des)*sin(fi0 + q0des) + l0*l1*l2^2*m0*m1*m2^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(q2des) - l0*l1*l2^2*m0*m1*m2^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(q2des) + 2*l0*l1*l2^2*m0*m1^2*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(q2des) - 2*l0*l1*l2^2*m0*m1^2*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(q2des) + 6*l0*l1*l2^2*m0^2*m1*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(q2des) - 6*l0*l1*l2^2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(q2des) + 10*l0^2*l1*l2*m0^2*m1*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des)*cos(fi0 - q1des) - 10*l0^2*l1*l2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des)*cos(fi0 - q1des) + 2*l0^2*l1*l2*m0^2*m1*m2*cos(q1des - fi0 + q2des)*cos(fi0 + q0des)*sin(q0des + q1des) - 2*l0^2*l1*l2*m0^2*m1*m2*cos(q1des - fi0 + q2des)*cos(q0des + q1des)*sin(fi0 + q0des) - 2*l0*l1^2*l2*m0*m1*m2^2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 2*l0*l1^2*l2*m0*m1*m2^2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) - 5*l0*l1^2*l2*m0*m1^2*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 5*l0*l1^2*l2*m0*m1^2*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) - 14*l0*l1^2*l2*m0^2*m1*m2*cos(q0des + q1des + q2des)*sin(q0des + q1des)*cos(fi0 - q1des) + 14*l0*l1^2*l2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(q0des + q1des)*cos(fi0 - q1des) + l0*l1^2*l2*m0*m1^2*m2*cos(fi0 + q0des)*sin(q0des + q1des)*cos(q2des) - l0*l1^2*l2*m0*m1^2*m2*cos(q0des + q1des)*sin(fi0 + q0des)*cos(q2des) + 4*l0*l1^2*l2*m0^2*m1*m2*cos(fi0 + q0des)*sin(q0des + q1des)*cos(q2des) - 4*l0*l1^2*l2*m0^2*m1*m2*cos(q0des + q1des)*sin(fi0 + q0des)*cos(q2des) - l0*l1*l2^2*m0*m1*m2^2*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*sin(q0des + q1des) + l0*l1*l2^2*m0*m1*m2^2*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*cos(q0des + q1des) - 2*l0*l1*l2^2*m0*m1^2*m2*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*sin(q0des + q1des) + 2*l0*l1*l2^2*m0*m1^2*m2*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*cos(q0des + q1des) - 6*l0*l1*l2^2*m0^2*m1*m2*cos(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*sin(q0des + q1des) + 6*l0*l1*l2^2*m0^2*m1*m2*sin(q0des + q1des + q2des)*cos(q1des - fi0 + q2des)*cos(q0des + q1des));
                          (((vy*(2*l2*m0*sin(q0des + q1des + q2des) + 2*l2*m1*sin(q0des + q1des + q2des) + l2*m2*sin(q0des + q1des + q2des) + 2*l0*m0*sin(fi0 + q0des) + 2*l1*m0*sin(q0des + q1des) + l1*m1*sin(q0des + q1des)))/(2*(m0 + m1 + m2)) + (vx*(2*l2*m0*cos(q0des + q1des + q2des) + 2*l2*m1*cos(q0des + q1des + q2des) + l2*m2*cos(q0des + q1des + q2des) + 2*l0*m0*cos(fi0 + q0des) + 2*l1*m0*cos(q0des + q1des) + l1*m1*cos(q0des + q1des)))/(2*(m0 + m1 + m2)))*(((2*l2*m0*cos(q0des + q1des + q2des) + 2*l2*m1*cos(q0des + q1des + q2des) + l2*m2*cos(q0des + q1des + q2des) + 2*l1*m0*cos(q0des + q1des) + l1*m1*cos(q0des + q1des))*(4*I0*m0 + 4*I0*m1 + 4*I1*m0 + 4*I0*m2 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + 4*l0^2*m0*m1 + 4*l0^2*m0*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 4*l0*l1*m0*cos(fi0 - q1des)*(m1 + 2*m2) + 4*l0*l2*m0*m2*cos(q1des - fi0 + q2des) + 4*l1*l2*m0*m2*cos(q2des) + 2*l1*l2*m1*m2*cos(q2des)))/(8*(m0 + m1 + m2)^2) - ((2*l2*m0*cos(q0des + q1des + q2des) + 2*l2*m1*cos(q0des + q1des + q2des) + l2*m2*cos(q0des + q1des + q2des) + 2*l0*m0*cos(fi0 + q0des) + 2*l1*m0*cos(q0des + q1des) + l1*m1*cos(q0des + q1des))*(4*I1*m0 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l0*l1*m0*cos(fi0 - q1des)*(m1 + 2*m2) + 2*l0*l2*m0*m2*cos(q1des - fi0 + q2des) + 4*l1*l2*m0*m2*cos(q2des) + 2*l1*l2*m1*m2*cos(q2des)))/(8*(m0 + m1 + m2)^2)) + ((vy*(4*I0*m0 + 4*I0*m1 + 4*I1*m0 + 4*I0*m2 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + 4*l0^2*m0*m1 + 4*l0^2*m0*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 4*l0*l1*m0*cos(fi0 - q1des)*(m1 + 2*m2) + 4*l0*l2*m0*m2*cos(q1des - fi0 + q2des) + 4*l1*l2*m0*m2*cos(q2des) + 2*l1*l2*m1*m2*cos(q2des)))/(4*(m0 + m1 + m2)) - (k*(2*l2*m0*cos(q0des + q1des + q2des) + 2*l2*m1*cos(q0des + q1des + q2des) + l2*m2*cos(q0des + q1des + q2des) + 2*l0*m0*cos(fi0 + q0des) + 2*l1*m0*cos(q0des + q1des) + l1*m1*cos(q0des + q1des)))/(2*(m0 + m1 + m2)))*((l0*l1*m0^2*cos(fi0 + q0des)*sin(q0des + q1des))/(m0 + m1 + m2)^2 - (l0*l1*m0^2*cos(q0des + q1des)*sin(fi0 + q0des))/(m0 + m1 + m2)^2 - (l0*l2*m0^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des))/(m0 + m1 + m2)^2 + (l0*l2*m0^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des))/(m0 + m1 + m2)^2 + (l0*l1*m0*m1*cos(fi0 + q0des)*sin(q0des + q1des))/(2*(m0 + m1 + m2)^2) - (l0*l1*m0*m1*cos(q0des + q1des)*sin(fi0 + q0des))/(2*(m0 + m1 + m2)^2) - (l0*l2*m0*m1*cos(q0des + q1des + q2des)*sin(fi0 + q0des))/(m0 + m1 + m2)^2 + (l0*l2*m0*m1*sin(q0des + q1des + q2des)*cos(fi0 + q0des))/(m0 + m1 + m2)^2 - (l0*l2*m0*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des))/(2*(m0 + m1 + m2)^2) + (l0*l2*m0*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des))/(2*(m0 + m1 + m2)^2)))/((((2*l2*m0*cos(q0des + q1des + q2des) + 2*l2*m1*cos(q0des + q1des + q2des) + l2*m2*cos(q0des + q1des + q2des))*(4*I0*m0 + 4*I0*m1 + 4*I1*m0 + 4*I0*m2 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + 4*l0^2*m0*m1 + 4*l0^2*m0*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 4*l0*l1*m0*cos(fi0 - q1des)*(m1 + 2*m2) + 4*l0*l2*m0*m2*cos(q1des - fi0 + q2des) + 4*l1*l2*m0*m2*cos(q2des) + 2*l1*l2*m1*m2*cos(q2des)))/(8*(m0 + m1 + m2)^2) - ((2*l2*m0*cos(q0des + q1des + q2des) + 2*l2*m1*cos(q0des + q1des + q2des) + l2*m2*cos(q0des + q1des + q2des) + 2*l0*m0*cos(fi0 + q0des) + 2*l1*m0*cos(q0des + q1des) + l1*m1*cos(q0des + q1des))*(4*I2*m0 + 4*I2*m1 + 4*I2*m2 + l2^2*m0*m2 + l2^2*m1*m2 + l1*l2*m2*cos(q2des)*(2*m0 + m1) + 2*l0*l2*m0*m2*cos(q1des - fi0 + q2des)))/(8*(m0 + m1 + m2)^2))*((l0*l1*m0^2*cos(fi0 + q0des)*sin(q0des + q1des))/(m0 + m1 + m2)^2 - (l0*l1*m0^2*cos(q0des + q1des)*sin(fi0 + q0des))/(m0 + m1 + m2)^2 - (l0*l2*m0^2*cos(q0des + q1des + q2des)*sin(fi0 + q0des))/(m0 + m1 + m2)^2 + (l0*l2*m0^2*sin(q0des + q1des + q2des)*cos(fi0 + q0des))/(m0 + m1 + m2)^2 + (l0*l1*m0*m1*cos(fi0 + q0des)*sin(q0des + q1des))/(2*(m0 + m1 + m2)^2) - (l0*l1*m0*m1*cos(q0des + q1des)*sin(fi0 + q0des))/(2*(m0 + m1 + m2)^2) - (l0*l2*m0*m1*cos(q0des + q1des + q2des)*sin(fi0 + q0des))/(m0 + m1 + m2)^2 + (l0*l2*m0*m1*sin(q0des + q1des + q2des)*cos(fi0 + q0des))/(m0 + m1 + m2)^2 - (l0*l2*m0*m2*cos(q0des + q1des + q2des)*sin(fi0 + q0des))/(2*(m0 + m1 + m2)^2) + (l0*l2*m0*m2*sin(q0des + q1des + q2des)*cos(fi0 + q0des))/(2*(m0 + m1 + m2)^2)) - (((2*l2*m0*cos(q0des + q1des + q2des) + 2*l2*m1*cos(q0des + q1des + q2des) + l2*m2*cos(q0des + q1des + q2des) + 2*l1*m0*cos(q0des + q1des) + l1*m1*cos(q0des + q1des))*(4*I0*m0 + 4*I0*m1 + 4*I1*m0 + 4*I0*m2 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + 4*l0^2*m0*m1 + 4*l0^2*m0*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 4*l0*l1*m0*cos(fi0 - q1des)*(m1 + 2*m2) + 4*l0*l2*m0*m2*cos(q1des - fi0 + q2des) + 4*l1*l2*m0*m2*cos(q2des) + 2*l1*l2*m1*m2*cos(q2des)))/(8*(m0 + m1 + m2)^2) - ((2*l2*m0*cos(q0des + q1des + q2des) + 2*l2*m1*cos(q0des + q1des + q2des) + l2*m2*cos(q0des + q1des + q2des) + 2*l0*m0*cos(fi0 + q0des) + 2*l1*m0*cos(q0des + q1des) + l1*m1*cos(q0des + q1des))*(4*I1*m0 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l0*l1*m0*cos(fi0 - q1des)*(m1 + 2*m2) + 2*l0*l2*m0*m2*cos(q1des - fi0 + q2des) + 4*l1*l2*m0*m2*cos(q2des) + 2*l1*l2*m1*m2*cos(q2des)))/(8*(m0 + m1 + m2)^2))*(((2*l2*m0*sin(q0des + q1des + q2des) + 2*l2*m1*sin(q0des + q1des + q2des) + l2*m2*sin(q0des + q1des + q2des))*(2*l2*m0*cos(q0des + q1des + q2des) + 2*l2*m1*cos(q0des + q1des + q2des) + l2*m2*cos(q0des + q1des + q2des) + 2*l0*m0*cos(fi0 + q0des) + 2*l1*m0*cos(q0des + q1des) + l1*m1*cos(q0des + q1des)))/(4*(m0 + m1 + m2)^2) - ((2*l2*m0*cos(q0des + q1des + q2des) + 2*l2*m1*cos(q0des + q1des + q2des) + l2*m2*cos(q0des + q1des + q2des))*(2*l2*m0*sin(q0des + q1des + q2des) + 2*l2*m1*sin(q0des + q1des + q2des) + l2*m2*sin(q0des + q1des + q2des) + 2*l0*m0*sin(fi0 + q0des) + 2*l1*m0*sin(q0des + q1des) + l1*m1*sin(q0des + q1des)))/(4*(m0 + m1 + m2)^2)))];
                
                ddq_des = (dq_des - stateVariables(4:6, step-1))/timeStep;   
            end
            desiredState = vertcat(q_des, dq_des);
            system.refreshMatricesForRK4(desiredState);
            
            usol = system.H \ (system.M * ddq_des + system.C);
            system.refreshMatrices();
            u = usol;
            
        end
    end
    
end

