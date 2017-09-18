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
            
           vx = -0.1*2*sin(2*timeVector(step));
           vy = 0.1*2*cos(2*timeVector(step));
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
                dq_des = [(16*I1*l2*m0^2*vx*cos(fi0 + q0des + q1des + q2des) - 8*k*l1*l2*m1^2*sin(q2des) - 16*I2*l1*m0^2*vy*sin(fi0 + q0des + q1des) - 8*I2*l1*m1^2*vy*sin(fi0 + q0des + q1des) - l1*vx*cos(fi0 + q0des + q1des)*(2*m0 + m1)*(8*I2*(m0 + m1 + m2) - l2^2*m2^2) - 16*k*l1*l2*m0^2*sin(q2des) + 16*I1*l2*m1^2*vx*cos(fi0 + q0des + q1des + q2des) + 8*I1*l2*m2^2*vx*cos(fi0 + q0des + q1des + q2des) + 16*I1*l2*m0^2*vy*sin(fi0 + q0des + q1des + q2des) + 16*I1*l2*m1^2*vy*sin(fi0 + q0des + q1des + q2des) + 8*I1*l2*m2^2*vy*sin(fi0 + q0des + q1des + q2des) - 4*l1^2*l2*m0^2*m2*vx*cos(fi0 + q0des + q1des - q2des) + 2*l1*l2^2*m0*m2^2*vx*cos(fi0 + q0des + q1des + 2*q2des) + 4*l1*l2^2*m0^2*m2*vx*cos(fi0 + q0des + q1des + 2*q2des) - l1^2*l2*m1^2*m2*vx*cos(fi0 + q0des + q1des - q2des) + l1*l2^2*m1*m2^2*vx*cos(fi0 + q0des + q1des + 2*q2des) + 2*l1*l2^2*m1^2*m2*vx*cos(fi0 + q0des + q1des + 2*q2des) - 4*l1^2*l2*m0^2*m2*vy*sin(fi0 + q0des + q1des - q2des) + 2*l1*l2^2*m0*m2^2*vy*sin(fi0 + q0des + q1des + 2*q2des) + 4*l1*l2^2*m0^2*m2*vy*sin(fi0 + q0des + q1des + 2*q2des) - l1^2*l2*m1^2*m2*vy*sin(fi0 + q0des + q1des - q2des) + l1*l2^2*m1*m2^2*vy*sin(fi0 + q0des + q1des + 2*q2des) + 2*l1*l2^2*m1^2*m2*vy*sin(fi0 + q0des + q1des + 2*q2des) - 24*k*l1*l2*m0*m1*sin(q2des) - 8*k*l1*l2*m0*m2*sin(q2des) - 4*k*l1*l2*m1*m2*sin(q2des) + 2*l1*l2^2*m0*m2^2*vy*sin(fi0 + q0des + q1des) + l1*l2^2*m1*m2^2*vy*sin(fi0 + q0des + q1des) - 24*I2*l1*m0*m1*vy*sin(fi0 + q0des + q1des) - 16*I2*l1*m0*m2*vy*sin(fi0 + q0des + q1des) - 8*I2*l1*m1*m2*vy*sin(fi0 + q0des + q1des) + 4*l1^2*l2*m0*m1^2*vx*cos(fi0 + q0des + q1des + q2des) + 4*l1^2*l2*m0^2*m1*vx*cos(fi0 + q0des + q1des + q2des) + 8*l1^2*l2*m0*m2^2*vx*cos(fi0 + q0des + q1des + q2des) + 12*l1^2*l2*m0^2*m2*vx*cos(fi0 + q0des + q1des + q2des) + 2*l1^2*l2*m1*m2^2*vx*cos(fi0 + q0des + q1des + q2des) + 3*l1^2*l2*m1^2*m2*vx*cos(fi0 + q0des + q1des + q2des) + 4*l1^2*l2*m0*m1^2*vy*sin(fi0 + q0des + q1des + q2des) + 4*l1^2*l2*m0^2*m1*vy*sin(fi0 + q0des + q1des + q2des) + 8*l1^2*l2*m0*m2^2*vy*sin(fi0 + q0des + q1des + q2des) + 12*l1^2*l2*m0^2*m2*vy*sin(fi0 + q0des + q1des + q2des) + 2*l1^2*l2*m1*m2^2*vy*sin(fi0 + q0des + q1des + q2des) + 3*l1^2*l2*m1^2*m2*vy*sin(fi0 + q0des + q1des + q2des) + 32*I1*l2*m0*m1*vx*cos(fi0 + q0des + q1des + q2des) + 24*I1*l2*m0*m2*vx*cos(fi0 + q0des + q1des + q2des) + 24*I1*l2*m1*m2*vx*cos(fi0 + q0des + q1des + q2des) + 32*I1*l2*m0*m1*vy*sin(fi0 + q0des + q1des + q2des) + 24*I1*l2*m0*m2*vy*sin(fi0 + q0des + q1des + q2des) + 24*I1*l2*m1*m2*vy*sin(fi0 + q0des + q1des + q2des) + 4*l0*l1*l2*m0*m1^2*vx*cos(fi0 + q0des + 2*q1des + q2des) + 4*l0*l1*l2*m0^2*m1*vx*cos(fi0 + q0des + 2*q1des + q2des) + 4*l0*l1*l2*m0*m2^2*vx*cos(fi0 + q0des + 2*q1des + q2des) + 4*l0*l1*l2*m0^2*m2*vx*cos(fi0 + q0des + 2*q1des + q2des) - 4*l1^2*l2*m0*m1*m2*vx*cos(fi0 + q0des + q1des - q2des) + 6*l1*l2^2*m0*m1*m2*vx*cos(fi0 + q0des + q1des + 2*q2des) + 4*l0*l1*l2*m0*m1^2*vy*sin(fi0 + q0des + 2*q1des + q2des) + 4*l0*l1*l2*m0^2*m1*vy*sin(fi0 + q0des + 2*q1des + q2des) + 4*l0*l1*l2*m0*m2^2*vy*sin(fi0 + q0des + 2*q1des + q2des) + 4*l0*l1*l2*m0^2*m2*vy*sin(fi0 + q0des + 2*q1des + q2des) - 4*l1^2*l2*m0*m1*m2*vy*sin(fi0 + q0des + q1des - q2des) + 6*l1*l2^2*m0*m1*m2*vy*sin(fi0 + q0des + q1des + 2*q2des) + 4*l0*l1*l2*m0*m1^2*vx*cos(fi0 + q0des + q2des) + 4*l0*l1*l2*m0^2*m1*vx*cos(fi0 + q0des + q2des) + 4*l0*l1*l2*m0*m2^2*vx*cos(fi0 + q0des + q2des) + 8*l0*l1*l2*m0^2*m2*vx*cos(fi0 + q0des + q2des) + 4*l0*l1*l2*m0*m1^2*vy*sin(fi0 + q0des + q2des) + 4*l0*l1*l2*m0^2*m1*vy*sin(fi0 + q0des + q2des) + 4*l0*l1*l2*m0*m2^2*vy*sin(fi0 + q0des + q2des) + 8*l0*l1*l2*m0^2*m2*vy*sin(fi0 + q0des + q2des) - 4*l0*l1*l2*m0^2*m2*vy*sin(fi0 + q0des - q2des) + 18*l1^2*l2*m0*m1*m2*vx*cos(fi0 + q0des + q1des + q2des) + 18*l1^2*l2*m0*m1*m2*vy*sin(fi0 + q0des + q1des + q2des) + 8*l0*l1*l2*m0*m1*m2*vx*cos(fi0 + q0des + 2*q1des + q2des) + 8*l0*l1*l2*m0*m1*m2*vy*sin(fi0 + q0des + 2*q1des + q2des) + 10*l0*l1*l2*m0*m1*m2*vx*cos(fi0 + q0des + q2des) + 10*l0*l1*l2*m0*m1*m2*vy*sin(fi0 + q0des + q2des) - 2*l0*l1*l2*m0*m1*m2*vy*sin(fi0 + q0des - q2des) - 2*l0*l1*l2*m0*m2*vx*cos(fi0 + q0des - q2des)*(2*m0 + m1))/(2*(l2*(l0*m0*(sin(q1des + q2des)*((2*m0*m2 - m1^2)*l1^2 + 4*I1*(2*m0 + 2*m1 + m2)) + 2*l0*l1*m0*sin(2*q1des + q2des)*(m1 + m2)) - 2*l1*sin(q2des)*(m0*(m0*(3*m1 + m2) + m1*(2*m1 + m2))*l0^2 + I0*(2*m0 + m1)*(2*m0 + 2*m1 + m2)) + l0*l1^2*m0*sin(q1des - q2des)*(m1 + m2)*(2*m0 + m1)) - l0*l1*m0*sin(q1des)*(2*m0 + m1)*(- m2*l2^2 + 4*I2)));
                         -(16*I0*l2*m0^2*vx*cos(fi0 + q0des + q1des + q2des) - 16*k*l0*l2*m0^2*sin(q1des + q2des) - 16*k*l1*l2*m0^2*sin(q2des) - 8*k*l1*l2*m1^2*sin(q2des) - 16*I2*l1*m0^2*vy*sin(fi0 + q0des + q1des) - 8*I2*l1*m1^2*vy*sin(fi0 + q0des + q1des) - l1*vx*cos(fi0 + q0des + q1des)*(2*m0 + m1)*(8*I2*(m0 + m1 + m2) - l2^2*m2^2) - 16*I2*l0*m0^2*vy*sin(fi0 + q0des) + 16*I0*l2*m1^2*vx*cos(fi0 + q0des + q1des + q2des) + 16*I1*l2*m0^2*vx*cos(fi0 + q0des + q1des + q2des) + 8*I0*l2*m2^2*vx*cos(fi0 + q0des + q1des + q2des) + 16*I1*l2*m1^2*vx*cos(fi0 + q0des + q1des + q2des) + 8*I1*l2*m2^2*vx*cos(fi0 + q0des + q1des + q2des) + 16*I0*l2*m0^2*vy*sin(fi0 + q0des + q1des + q2des) + 16*I0*l2*m1^2*vy*sin(fi0 + q0des + q1des + q2des) + 16*I1*l2*m0^2*vy*sin(fi0 + q0des + q1des + q2des) + 8*I0*l2*m2^2*vy*sin(fi0 + q0des + q1des + q2des) + 16*I1*l2*m1^2*vy*sin(fi0 + q0des + q1des + q2des) + 8*I1*l2*m2^2*vy*sin(fi0 + q0des + q1des + q2des) - 2*l0*m0*vx*cos(fi0 + q0des)*(8*I2*(m0 + m1 + m2) - l2^2*m2^2) - 16*I2*l0*m0*m1*vy*sin(fi0 + q0des) - 16*I2*l0*m0*m2*vy*sin(fi0 + q0des) - 4*l1^2*l2*m0^2*m2*vx*cos(fi0 + q0des + q1des - q2des) + 2*l1*l2^2*m0*m2^2*vx*cos(fi0 + q0des + q1des + 2*q2des) + 4*l1*l2^2*m0^2*m2*vx*cos(fi0 + q0des + q1des + 2*q2des) - l1^2*l2*m1^2*m2*vx*cos(fi0 + q0des + q1des - q2des) + l1*l2^2*m1*m2^2*vx*cos(fi0 + q0des + q1des + 2*q2des) + 2*l1*l2^2*m1^2*m2*vx*cos(fi0 + q0des + q1des + 2*q2des) - 16*k*l0*l2*m0*m1*sin(q1des + q2des) - 8*k*l0*l2*m0*m2*sin(q1des + q2des) - 4*l1^2*l2*m0^2*m2*vy*sin(fi0 + q0des + q1des - q2des) + 2*l1*l2^2*m0*m2^2*vy*sin(fi0 + q0des + q1des + 2*q2des) + 4*l1*l2^2*m0^2*m2*vy*sin(fi0 + q0des + q1des + 2*q2des) - l1^2*l2*m1^2*m2*vy*sin(fi0 + q0des + q1des - q2des) + l1*l2^2*m1*m2^2*vy*sin(fi0 + q0des + q1des + 2*q2des) + 2*l1*l2^2*m1^2*m2*vy*sin(fi0 + q0des + q1des + 2*q2des) - 24*k*l1*l2*m0*m1*sin(q2des) - 8*k*l1*l2*m0*m2*sin(q2des) - 4*k*l1*l2*m1*m2*sin(q2des) + 2*l1*l2^2*m0*m2^2*vy*sin(fi0 + q0des + q1des) + l1*l2^2*m1*m2^2*vy*sin(fi0 + q0des + q1des) - 24*I2*l1*m0*m1*vy*sin(fi0 + q0des + q1des) - 16*I2*l1*m0*m2*vy*sin(fi0 + q0des + q1des) - 8*I2*l1*m1*m2*vy*sin(fi0 + q0des + q1des) - 4*l0^2*l2*m0^2*m2*vx*cos(fi0 + q0des - q1des - q2des) + 2*l0*l2^2*m0*m2^2*vx*cos(fi0 + q0des + 2*q1des + 2*q2des) + 4*l0*l2^2*m0^2*m2*vx*cos(fi0 + q0des + 2*q1des + 2*q2des) - 4*l0^2*l2*m0^2*m2*vy*sin(fi0 + q0des - q1des - q2des) + 2*l0*l2^2*m0*m2^2*vy*sin(fi0 + q0des + 2*q1des + 2*q2des) + 4*l0*l2^2*m0^2*m2*vy*sin(fi0 + q0des + 2*q1des + 2*q2des) + 16*l0^2*l2*m0*m1^2*vx*cos(fi0 + q0des + q1des + q2des) + 16*l0^2*l2*m0^2*m1*vx*cos(fi0 + q0des + q1des + q2des) + 8*l0^2*l2*m0*m2^2*vx*cos(fi0 + q0des + q1des + q2des) + 12*l0^2*l2*m0^2*m2*vx*cos(fi0 + q0des + q1des + q2des) + 4*l1^2*l2*m0*m1^2*vx*cos(fi0 + q0des + q1des + q2des) + 4*l1^2*l2*m0^2*m1*vx*cos(fi0 + q0des + q1des + q2des) + 8*l1^2*l2*m0*m2^2*vx*cos(fi0 + q0des + q1des + q2des) + 12*l1^2*l2*m0^2*m2*vx*cos(fi0 + q0des + q1des + q2des) + 2*l1^2*l2*m1*m2^2*vx*cos(fi0 + q0des + q1des + q2des) + 3*l1^2*l2*m1^2*m2*vx*cos(fi0 + q0des + q1des + q2des) + 16*l0^2*l2*m0*m1^2*vy*sin(fi0 + q0des + q1des + q2des) + 16*l0^2*l2*m0^2*m1*vy*sin(fi0 + q0des + q1des + q2des) + 8*l0^2*l2*m0*m2^2*vy*sin(fi0 + q0des + q1des + q2des) + 12*l0^2*l2*m0^2*m2*vy*sin(fi0 + q0des + q1des + q2des) + 4*l1^2*l2*m0*m1^2*vy*sin(fi0 + q0des + q1des + q2des) + 4*l1^2*l2*m0^2*m1*vy*sin(fi0 + q0des + q1des + q2des) + 8*l1^2*l2*m0*m2^2*vy*sin(fi0 + q0des + q1des + q2des) + 12*l1^2*l2*m0^2*m2*vy*sin(fi0 + q0des + q1des + q2des) + 2*l1^2*l2*m1*m2^2*vy*sin(fi0 + q0des + q1des + q2des) + 3*l1^2*l2*m1^2*m2*vy*sin(fi0 + q0des + q1des + q2des) + 32*I0*l2*m0*m1*vx*cos(fi0 + q0des + q1des + q2des) + 24*I0*l2*m0*m2*vx*cos(fi0 + q0des + q1des + q2des) + 32*I1*l2*m0*m1*vx*cos(fi0 + q0des + q1des + q2des) + 24*I0*l2*m1*m2*vx*cos(fi0 + q0des + q1des + q2des) + 24*I1*l2*m0*m2*vx*cos(fi0 + q0des + q1des + q2des) + 24*I1*l2*m1*m2*vx*cos(fi0 + q0des + q1des + q2des) + 2*l0*l2^2*m0*m2^2*vy*sin(fi0 + q0des) + 32*I0*l2*m0*m1*vy*sin(fi0 + q0des + q1des + q2des) + 24*I0*l2*m0*m2*vy*sin(fi0 + q0des + q1des + q2des) + 32*I1*l2*m0*m1*vy*sin(fi0 + q0des + q1des + q2des) + 24*I0*l2*m1*m2*vy*sin(fi0 + q0des + q1des + q2des) + 24*I1*l2*m0*m2*vy*sin(fi0 + q0des + q1des + q2des) + 24*I1*l2*m1*m2*vy*sin(fi0 + q0des + q1des + q2des) + 8*l0*l1*l2*m0*m1^2*vx*cos(fi0 + q0des + 2*q1des + q2des) + 8*l0*l1*l2*m0^2*m1*vx*cos(fi0 + q0des + 2*q1des + q2des) + 8*l0*l1*l2*m0*m2^2*vx*cos(fi0 + q0des + 2*q1des + q2des) + 12*l0*l1*l2*m0^2*m2*vx*cos(fi0 + q0des + 2*q1des + q2des) - 4*l1^2*l2*m0*m1*m2*vx*cos(fi0 + q0des + q1des - q2des) + 6*l1*l2^2*m0*m1*m2*vx*cos(fi0 + q0des + q1des + 2*q2des) + 8*l0*l1*l2*m0*m1^2*vy*sin(fi0 + q0des + 2*q1des + q2des) + 8*l0*l1*l2*m0^2*m1*vy*sin(fi0 + q0des + 2*q1des + q2des) + 8*l0*l1*l2*m0*m2^2*vy*sin(fi0 + q0des + 2*q1des + q2des) + 12*l0*l1*l2*m0^2*m2*vy*sin(fi0 + q0des + 2*q1des + q2des) - 4*l1^2*l2*m0*m1*m2*vy*sin(fi0 + q0des + q1des - q2des) + 6*l1*l2^2*m0*m1*m2*vy*sin(fi0 + q0des + q1des + 2*q2des) + 8*l0*l1*l2*m0*m1^2*vx*cos(fi0 + q0des + q2des) + 8*l0*l1*l2*m0^2*m1*vx*cos(fi0 + q0des + q2des) + 8*l0*l1*l2*m0*m2^2*vx*cos(fi0 + q0des + q2des) + 12*l0*l1*l2*m0^2*m2*vx*cos(fi0 + q0des + q2des) + 8*l0*l1*l2*m0*m1^2*vy*sin(fi0 + q0des + q2des) + 8*l0*l1*l2*m0^2*m1*vy*sin(fi0 + q0des + q2des) + 8*l0*l1*l2*m0*m2^2*vy*sin(fi0 + q0des + q2des) + 12*l0*l1*l2*m0^2*m2*vy*sin(fi0 + q0des + q2des) + 4*l0*l2^2*m0*m1*m2*vx*cos(fi0 + q0des + 2*q1des + 2*q2des) + 4*l0*l2^2*m0*m1*m2*vy*sin(fi0 + q0des + 2*q1des + 2*q2des) - 8*l0*l1*l2*m0^2*m2*vx*cos(fi0 + q0des - q2des) - 8*l0*l1*l2*m0^2*m2*vy*sin(fi0 + q0des - q2des) + 24*l0^2*l2*m0*m1*m2*vx*cos(fi0 + q0des + q1des + q2des) + 18*l1^2*l2*m0*m1*m2*vx*cos(fi0 + q0des + q1des + q2des) + 24*l0^2*l2*m0*m1*m2*vy*sin(fi0 + q0des + q1des + q2des) + 18*l1^2*l2*m0*m1*m2*vy*sin(fi0 + q0des + q1des + q2des) + 18*l0*l1*l2*m0*m1*m2*vx*cos(fi0 + q0des + 2*q1des + q2des) + 18*l0*l1*l2*m0*m1*m2*vy*sin(fi0 + q0des + 2*q1des + q2des) + 18*l0*l1*l2*m0*m1*m2*vx*cos(fi0 + q0des + q2des) + 18*l0*l1*l2*m0*m1*m2*vy*sin(fi0 + q0des + q2des) - 4*l0*l1*l2*m0*m1*m2*vx*cos(fi0 + q0des - q2des) - 4*l0*l1*l2*m0*m1*m2*vy*sin(fi0 + q0des - q2des))/(2*(l2*(l0*m0*(sin(q1des + q2des)*((2*m0*m2 - m1^2)*l1^2 + 4*I1*(2*m0 + 2*m1 + m2)) + 2*l0*l1*m0*sin(2*q1des + q2des)*(m1 + m2)) - 2*l1*sin(q2des)*(m0*(m0*(3*m1 + m2) + m1*(2*m1 + m2))*l0^2 + I0*(2*m0 + m1)*(2*m0 + 2*m1 + m2)) + l0*l1^2*m0*sin(q1des - q2des)*(m1 + m2)*(2*m0 + m1)) - l0*l1*m0*sin(q1des)*(2*m0 + m1)*(- m2*l2^2 + 4*I2)));
                         (8*I0*l1*m0^2*vx*cos(fi0 + q0des + q1des) - 8*I2*l0*m0^2*vy*sin(fi0 + q0des) - 8*k*l0*l2*m0^2*sin(q1des + q2des) - l0*m0*vx*cos(fi0 + q0des)*(8*I1*(m0 + m1 + m2) + 8*I2*(m0 + m1 + m2) - l1^2*m1^2 - l2^2*m2^2 + 4*l1^2*m0*m2) - 8*k*l0*l1*m0^2*sin(q1des) - 8*I1*l0*m0^2*vy*sin(fi0 + q0des) + 4*I0*l1*m1^2*vx*cos(fi0 + q0des + q1des) + 8*I0*l1*m0^2*vy*sin(fi0 + q0des + q1des) + 4*I0*l1*m1^2*vy*sin(fi0 + q0des + q1des) + 8*I0*l2*m0^2*vx*cos(fi0 + q0des + q1des + q2des) + 8*I0*l2*m1^2*vx*cos(fi0 + q0des + q1des + q2des) + 4*I0*l2*m2^2*vx*cos(fi0 + q0des + q1des + q2des) + 8*I0*l2*m0^2*vy*sin(fi0 + q0des + q1des + q2des) + 8*I0*l2*m1^2*vy*sin(fi0 + q0des + q1des + q2des) + 4*I0*l2*m2^2*vy*sin(fi0 + q0des + q1des + q2des) - 8*I1*l0*m0*m1*vy*sin(fi0 + q0des) - 8*I1*l0*m0*m2*vy*sin(fi0 + q0des) - 8*I2*l0*m0*m1*vy*sin(fi0 + q0des) - 8*I2*l0*m0*m2*vy*sin(fi0 + q0des) - 8*k*l0*l2*m0*m1*sin(q1des + q2des) - 4*k*l0*l2*m0*m2*sin(q1des + q2des) + 4*l0^2*l1*m0*m1^2*vx*cos(fi0 + q0des + q1des) + 6*l0^2*l1*m0^2*m1*vx*cos(fi0 + q0des + q1des) + 4*l0^2*l1*m0^2*m2*vx*cos(fi0 + q0des + q1des) - 4*k*l0*l1*m0*m1*sin(q1des) + 4*l0^2*l1*m0*m1^2*vy*sin(fi0 + q0des + q1des) + 6*l0^2*l1*m0^2*m1*vy*sin(fi0 + q0des + q1des) + 4*l0^2*l1*m0^2*m2*vy*sin(fi0 + q0des + q1des) + 12*I0*l1*m0*m1*vx*cos(fi0 + q0des + q1des) + 8*I0*l1*m0*m2*vx*cos(fi0 + q0des + q1des) + 4*I0*l1*m1*m2*vx*cos(fi0 + q0des + q1des) + 12*I0*l1*m0*m1*vy*sin(fi0 + q0des + q1des) + 8*I0*l1*m0*m2*vy*sin(fi0 + q0des + q1des) + 4*I0*l1*m1*m2*vy*sin(fi0 + q0des + q1des) - 2*l0^2*l2*m0^2*m2*vx*cos(fi0 + q0des - q1des - q2des) + l0*l2^2*m0*m2^2*vx*cos(fi0 + q0des + 2*q1des + 2*q2des) + 2*l0*l2^2*m0^2*m2*vx*cos(fi0 + q0des + 2*q1des + 2*q2des) - 2*l0^2*l2*m0^2*m2*vy*sin(fi0 + q0des - q1des - q2des) + l0*l2^2*m0*m2^2*vy*sin(fi0 + q0des + 2*q1des + 2*q2des) + 2*l0*l2^2*m0^2*m2*vy*sin(fi0 + q0des + 2*q1des + 2*q2des) + l0*l1^2*m0*m1^2*vx*cos(fi0 + q0des + 2*q1des) + 2*l0*l1^2*m0^2*m1*vx*cos(fi0 + q0des + 2*q1des) + 4*l0*l1^2*m0^2*m2*vx*cos(fi0 + q0des + 2*q1des) - 2*l0^2*l1*m0^2*m1*vy*sin(fi0 + q0des - q1des) + l0*l1^2*m0*m1^2*vy*sin(fi0 + q0des + 2*q1des) + 2*l0*l1^2*m0^2*m1*vy*sin(fi0 + q0des + 2*q1des) - 4*l0^2*l1*m0^2*m2*vy*sin(fi0 + q0des - q1des) + 4*l0*l1^2*m0^2*m2*vy*sin(fi0 + q0des + 2*q1des) + 8*l0^2*l2*m0*m1^2*vx*cos(fi0 + q0des + q1des + q2des) + 8*l0^2*l2*m0^2*m1*vx*cos(fi0 + q0des + q1des + q2des) + 4*l0^2*l2*m0*m2^2*vx*cos(fi0 + q0des + q1des + q2des) + 6*l0^2*l2*m0^2*m2*vx*cos(fi0 + q0des + q1des + q2des) + 8*l0^2*l2*m0*m1^2*vy*sin(fi0 + q0des + q1des + q2des) + 8*l0^2*l2*m0^2*m1*vy*sin(fi0 + q0des + q1des + q2des) + 4*l0^2*l2*m0*m2^2*vy*sin(fi0 + q0des + q1des + q2des) + 6*l0^2*l2*m0^2*m2*vy*sin(fi0 + q0des + q1des + q2des) - 2*l0^2*l1*m0^2*vx*cos(fi0 + q0des - q1des)*(m1 + 2*m2) + 16*I0*l2*m0*m1*vx*cos(fi0 + q0des + q1des + q2des) + 12*I0*l2*m0*m2*vx*cos(fi0 + q0des + q1des + q2des) + 12*I0*l2*m1*m2*vx*cos(fi0 + q0des + q1des + q2des) + l0*l1^2*m0*m1^2*vy*sin(fi0 + q0des) - 4*l0*l1^2*m0^2*m2*vy*sin(fi0 + q0des) + l0*l2^2*m0*m2^2*vy*sin(fi0 + q0des) + 16*I0*l2*m0*m1*vy*sin(fi0 + q0des + q1des + q2des) + 12*I0*l2*m0*m2*vy*sin(fi0 + q0des + q1des + q2des) + 12*I0*l2*m1*m2*vy*sin(fi0 + q0des + q1des + q2des) + 2*l0*l1*l2*m0*m1^2*vx*cos(fi0 + q0des + 2*q1des + q2des) + 2*l0*l1*l2*m0^2*m1*vx*cos(fi0 + q0des + 2*q1des + q2des) + 2*l0*l1*l2*m0*m2^2*vx*cos(fi0 + q0des + 2*q1des + q2des) + 6*l0*l1*l2*m0^2*m2*vx*cos(fi0 + q0des + 2*q1des + q2des) + 2*l0*l1*l2*m0*m1^2*vy*sin(fi0 + q0des + 2*q1des + q2des) + 2*l0*l1*l2*m0^2*m1*vy*sin(fi0 + q0des + 2*q1des + q2des) + 2*l0*l1*l2*m0*m2^2*vy*sin(fi0 + q0des + 2*q1des + q2des) + 6*l0*l1*l2*m0^2*m2*vy*sin(fi0 + q0des + 2*q1des + q2des) + 2*l0*l1*l2*m0*m1^2*vx*cos(fi0 + q0des + q2des) + 2*l0*l1*l2*m0^2*m1*vx*cos(fi0 + q0des + q2des) + 2*l0*l1*l2*m0*m2^2*vx*cos(fi0 + q0des + q2des) + 4*l0^2*l1*m0*m1*m2*vx*cos(fi0 + q0des + q1des) + 2*l0*l1*l2*m0*m1^2*vy*sin(fi0 + q0des + q2des) + 2*l0*l1*l2*m0^2*m1*vy*sin(fi0 + q0des + q2des) + 2*l0*l1*l2*m0*m2^2*vy*sin(fi0 + q0des + q2des) + 4*l0^2*l1*m0*m1*m2*vy*sin(fi0 + q0des + q1des) + 2*l0*l2^2*m0*m1*m2*vx*cos(fi0 + q0des + 2*q1des + 2*q2des) + 2*l0*l2^2*m0*m1*m2*vy*sin(fi0 + q0des + 2*q1des + 2*q2des) - 2*l0*l1*l2*m0^2*m2*vx*cos(fi0 + q0des - q2des) + 2*l0*l1^2*m0*m1*m2*vx*cos(fi0 + q0des + 2*q1des) - 2*l0*l1*l2*m0^2*m2*vy*sin(fi0 + q0des - q2des) + 2*l0*l1^2*m0*m1*m2*vy*sin(fi0 + q0des + 2*q1des) + 12*l0^2*l2*m0*m1*m2*vx*cos(fi0 + q0des + q1des + q2des) + 12*l0^2*l2*m0*m1*m2*vy*sin(fi0 + q0des + q1des + q2des) + 6*l0*l1*l2*m0*m1*m2*vx*cos(fi0 + q0des + 2*q1des + q2des) + 6*l0*l1*l2*m0*m1*m2*vy*sin(fi0 + q0des + 2*q1des + q2des) + 3*l0*l1*l2*m0*m1*m2*vx*cos(fi0 + q0des + q2des) + 3*l0*l1*l2*m0*m1*m2*vy*sin(fi0 + q0des + q2des) - l0*l1*l2*m0*m1*m2*vx*cos(fi0 + q0des - q2des) - l0*l1*l2*m0*m1*m2*vy*sin(fi0 + q0des - q2des))/(l2*(l0*m0*(sin(q1des + q2des)*((2*m0*m2 - m1^2)*l1^2 + 4*I1*(2*m0 + 2*m1 + m2)) + 2*l0*l1*m0*sin(2*q1des + q2des)*(m1 + m2)) - 2*l1*sin(q2des)*(m0*(m0*(3*m1 + m2) + m1*(2*m1 + m2))*l0^2 + I0*(2*m0 + m1)*(2*m0 + 2*m1 + m2)) + l0*l1^2*m0*sin(q1des - q2des)*(m1 + m2)*(2*m0 + m1)) - l0*l1*m0*sin(q1des)*(2*m0 + m1)*(- m2*l2^2 + 4*I2))];
                
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

