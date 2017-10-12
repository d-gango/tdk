classdef robotModel < handle
    properties
        angMom
        q
        qd
        u
        M
        C
        H
        %symbolic variables
        q_sym
        qd_sym
        qdd_sym
        M_sym
        C_sym
        H_sym
        u_sym
        endEffector_sym
        endEffectorVel_sym
    end
    
    properties (Constant)
        m0 = 40
        m1 = 4
        m2 = 3
        l0 = 1
        l1 = 1
        l2 = 1
        I0 = 6.667
        I1 = 0.333
        I2 = 0.25
        fi0 = pi/6
    end
    
    methods
        % constructor
        function obj = robotModel(state)    
            obj.q = zeros(3,1);
            obj.qd = zeros(3,1);
            obj.q(1) = state(1);
            obj.q(2) = state(2);
            obj.q(3) = state(3);
            obj.qd(1) = state(4);
            obj.qd(2) = state(5);
            obj.qd(3) = state(6);
            obj.u = [1; 0];
            refreshMatrices(obj);
            
            makeSymbolicVariables(obj);
            obj.angMom = angularMomentum(obj);
        end
        
        % calculate system matrices
        function refreshMatrices(obj)
            m0 = obj.m0;
            m1 = obj.m1;
            m2 = obj.m2;
            l0 = obj.l0;
            l1 = obj.l1;
            l2 = obj.l2;
            I0 = obj.I0;
            I1 = obj.I1;
            I2 = obj.I2;
            fi0 = obj.fi0;
            q0 = obj.q(1);
            q1 = obj.q(2);
            q2 = obj.q(3);
            dq0 = obj.qd(1);
            dq1 = obj.qd(2);
            dq2 = obj.qd(3);
            
            m11 = (4*I0*m0 + 4*I0*m1 + 4*I1*m0 + 4*I0*m2 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + 4*l0^2*m0*m1 + 4*l0^2*m0*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 4*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 4*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2))/(4*(m0 + m1 + m2));
            m12 = (4*I1*m0 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2))/(4*(m0 + m1 + m2));
            m13 = (4*I2*m0 + 4*I2*m1 + 4*I2*m2 + l2^2*m0*m2 + l2^2*m1*m2 + l1*l2*m2*cos(q2)*(2*m0 + m1) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2))/(4*(m0 + m1 + m2));
            m21 = (4*I1*m0 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2))/(4*(m0 + m1 + m2));
            m22 = (4*I2*m0 + 4*I2*m1 + 4*I2*m2 + 4*I1*(m0 + m1 + m2) + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l1*l2*m2*cos(q2)*(2*m0 + m1))/(4*(m0 + m1 + m2));
            m23 = (m2*(m0 + m1)*l2^2 + l1*m2*cos(q2)*(2*m0 + m1)*l2 + 4*I2*(m0 + m1 + m2))/(4*(m0 + m1 + m2));
            m31 = (4*I2*m0 + 4*I2*m1 + 4*I2*m2 + l2^2*m0*m2 + l2^2*m1*m2 + l1*l2*m2*cos(q2)*(2*m0 + m1) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2))/(4*(m0 + m1 + m2));
            m32 = (m2*(m0 + m1)*l2^2 + l1*m2*cos(q2)*(2*m0 + m1)*l2 + 4*I2*(m0 + m1 + m2))/(4*(m0 + m1 + m2));
            m33 = (m2*(m0 + m1)*l2^2 + 4*I2*(m0 + m1 + m2))/(4*(m0 + m1 + m2));
            
            c1 = (dq0*(4*dq1*l0*m0*(l1*sin(fi0 - q1)*(m1 + 2*m2) - l2*m2*sin(q1 - fi0 + q2)) - 2*dq2*l2*m2*(l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2))) - dq2^2*l2*m2*(l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2)) + 2*dq1^2*l0*m0*(l1*sin(fi0 - q1)*(m1 + 2*m2) - l2*m2*sin(q1 - fi0 + q2)) - 2*dq1*dq2*l2*m2*(l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2)))/(4*(m0 + m1 + m2));
            c2 = -(2*l0*m0*(l1*sin(fi0 - q1)*(m1 + 2*m2) - l2*m2*sin(q1 - fi0 + q2))*dq0^2 + 2*dq2*l1*l2*m2*sin(q2)*(2*m0 + m1)*dq0 + dq2*l1*l2*m2*sin(q2)*(2*dq1 + dq2)*(2*m0 + m1))/(4*(m0 + m1 + m2));
            c3 = (l2*m2*((l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2))*dq0^2 + 2*l1*sin(q2)*(2*m0 + m1)*dq0*dq1 + l1*sin(q2)*(2*m0 + m1)*dq1^2))/(4*(m0 + m1 + m2));

                        
            obj.M = [m11, m12, m13;...
                     m21, m22, m23;...
                     m31, m32, m33];
                
            obj.C = [c1;...
                     c2;...
                     c3];
         
            obj.H = [0, 0;
                     1, 0;
                     0, 1];
        end
        
         function refreshMatricesForRK4(obj, state)
            m0 = obj.m0;
            m1 = obj.m1;
            m2 = obj.m2;
            l0 = obj.l0;
            l1 = obj.l1;
            l2 = obj.l2;
            I0 = obj.I0;
            I1 = obj.I1;
            I2 = obj.I2;
            fi0 = obj.fi0;
            q0 = state(1);
            q1 = state(2);
            q2 = state(3);
            dq0 = state(4);
            dq1 = state(5);
            dq2 = state(6);
            
            m11 = (4*I0*m0 + 4*I0*m1 + 4*I1*m0 + 4*I0*m2 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + 4*l0^2*m0*m1 + 4*l0^2*m0*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 4*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 4*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2))/(4*(m0 + m1 + m2));
            m12 = (4*I1*m0 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2))/(4*(m0 + m1 + m2));
            m13 = (4*I2*m0 + 4*I2*m1 + 4*I2*m2 + l2^2*m0*m2 + l2^2*m1*m2 + l1*l2*m2*cos(q2)*(2*m0 + m1) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2))/(4*(m0 + m1 + m2));
            m21 = (4*I1*m0 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2))/(4*(m0 + m1 + m2));
            m22 = (4*I2*m0 + 4*I2*m1 + 4*I2*m2 + 4*I1*(m0 + m1 + m2) + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l1*l2*m2*cos(q2)*(2*m0 + m1))/(4*(m0 + m1 + m2));
            m23 = (m2*(m0 + m1)*l2^2 + l1*m2*cos(q2)*(2*m0 + m1)*l2 + 4*I2*(m0 + m1 + m2))/(4*(m0 + m1 + m2));
            m31 = (4*I2*m0 + 4*I2*m1 + 4*I2*m2 + l2^2*m0*m2 + l2^2*m1*m2 + l1*l2*m2*cos(q2)*(2*m0 + m1) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2))/(4*(m0 + m1 + m2));
            m32 = (m2*(m0 + m1)*l2^2 + l1*m2*cos(q2)*(2*m0 + m1)*l2 + 4*I2*(m0 + m1 + m2))/(4*(m0 + m1 + m2));
            m33 = (m2*(m0 + m1)*l2^2 + 4*I2*(m0 + m1 + m2))/(4*(m0 + m1 + m2));
            
            c1 = (dq0*(4*dq1*l0*m0*(l1*sin(fi0 - q1)*(m1 + 2*m2) - l2*m2*sin(q1 - fi0 + q2)) - 2*dq2*l2*m2*(l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2))) - dq2^2*l2*m2*(l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2)) + 2*dq1^2*l0*m0*(l1*sin(fi0 - q1)*(m1 + 2*m2) - l2*m2*sin(q1 - fi0 + q2)) - 2*dq1*dq2*l2*m2*(l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2)))/(4*(m0 + m1 + m2));
            c2 = -(2*l0*m0*(l1*sin(fi0 - q1)*(m1 + 2*m2) - l2*m2*sin(q1 - fi0 + q2))*dq0^2 + 2*dq2*l1*l2*m2*sin(q2)*(2*m0 + m1)*dq0 + dq2*l1*l2*m2*sin(q2)*(2*dq1 + dq2)*(2*m0 + m1))/(4*(m0 + m1 + m2));
            c3 = (l2*m2*((l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2))*dq0^2 + 2*l1*sin(q2)*(2*m0 + m1)*dq0*dq1 + l1*sin(q2)*(2*m0 + m1)*dq1^2))/(4*(m0 + m1 + m2));

                        
            obj.M = [m11, m12, m13;...
                     m21, m22, m23;...
                     m31, m32, m33];
                
            obj.C = [c1;...
                     c2;...
                     c3];
         
            obj.H = [0, 0;
                     1, 0;
                     0, 1];
        end
        
        % calculate system derivatives (equation of motion)
        function deriv = derivatives(obj, state,~)
            obj.refreshMatricesForRK4(state);
            qdd = obj.M \ (obj.H*obj.u - obj.C);
            deriv = vertcat(obj.qd, qdd);
            obj.refreshMatrices();
        end
        
        % integrates 1 timestep
        function integrate(obj)
            state = obj.getStateVariables();
            newState = RK4(@obj.derivatives, state);
            
            obj.q = newState(1:3);
            obj.qd = newState(4:6);
        end
        
        % state variables getter
        function s = getStateVariables(obj)
            s = vertcat(obj.q, obj.qd);
        end
        
        % animates the simulation
        function animate(obj)
            global stateVariables timeVector
            
            figure
            axesHandle = gca;
            xlim(axesHandle, [-5, 5] );
            ylim(axesHandle, [-5, 5] );
            axis equal;
            
            m0 = obj.m0;
            m1 = obj.m1;
            m2 = obj.m2;
            l0 = obj.l0;
            l1 = obj.l1;
            l2 = obj.l2;
            I0 = obj.I0;
            I1 = obj.I1;
            I2 = obj.I2;
            fi0 = obj.fi0;
            
            q0 = stateVariables(1,1);
            q1 = stateVariables(2,1);
            q2 = stateVariables(3,1);
            
            baseCenter = [ -(l1*cos(q0 + q1)*(m1 + 2*m2) + l2*m2*cos(q0 + q1 + q2) + 2*l0*cos(fi0 + q0)*(m1 + m2))/(2*(m0 + m1 + m2)), -(l1*sin(q0 + q1)*(m1 + 2*m2) + l2*m2*sin(q0 + q1 + q2) + 2*l0*sin(fi0 + q0)*(m1 + m2))/(2*(m0 + m1 + m2))];
            baseCorners = [ baseCenter(1) - l0*cos(fi0+q0), baseCenter(2) - l0*sin(fi0+q0);...
                            baseCenter(1) + l0*cos(fi0-q0), baseCenter(2) - l0*sin(fi0-q0);...
                            baseCenter(1) + l0*cos(fi0+q0), baseCenter(2) + l0*sin(fi0+q0);...
                            baseCenter(1) - l0*cos(fi0-q0), baseCenter(2) + l0*sin(fi0-q0)];
            baseFaces = [1, 2, 3, 4];
            base = patch('Vertices',baseCorners, 'Faces', baseFaces, 'FaceColor', 'r');
            
            link1x = [ -(l1*cos(q0 + q1)*(m1 + 2*m2) + l2*m2*cos(q0 + q1 + q2) - 2*l0*m0*cos(fi0 + q0))/(2*(m0 + m1 + m2)), (l1*cos(q0 + q1)*(2*m0 + m1) - l2*m2*cos(q0 + q1 + q2) + 2*l0*m0*cos(fi0 + q0))/(2*(m0 + m1 + m2))];
            link1y = [ -(l1*sin(q0 + q1)*(m1 + 2*m2) + l2*m2*sin(q0 + q1 + q2) - 2*l0*m0*sin(fi0 + q0))/(2*(m0 + m1 + m2)), (l1*sin(q0 + q1)*(2*m0 + m1) - l2*m2*sin(q0 + q1 + q2) + 2*l0*m0*sin(fi0 + q0))/(2*(m0 + m1 + m2))];
            link1 = line('XData', link1x, 'YData', link1y);
            
            link2x = [ link1x(2), (l1*cos(q0 + q1)*(2*m0 + m1) + l2*cos(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*cos(fi0 + q0))/(2*(m0 + m1 + m2)) ];
            link2y = [ link1y(2), (l1*sin(q0 + q1)*(2*m0 + m1) + l2*sin(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*sin(fi0 + q0))/(2*(m0 + m1 + m2)) ];
            link2 = line('XData', link2x, 'YData', link2y);
            
            r = 0.1;
            centerOfMass = [-r, -r, 2*r, 2*r];
            centerCircle = rectangle('Position', centerOfMass, 'Curvature', [1, 1], 'FaceColor', 'b');
            
            time = text(4,-4,num2str(timeVector(1)));
            
            for i=1:50:length(timeVector)
                q0 = stateVariables(1,i);
                q1 = stateVariables(2,i);
                q2 = stateVariables(3,i);
                
                baseCenter = [ -(l1*cos(q0 + q1)*(m1 + 2*m2) + l2*m2*cos(q0 + q1 + q2) + 2*l0*cos(fi0 + q0)*(m1 + m2))/(2*(m0 + m1 + m2)), -(l1*sin(q0 + q1)*(m1 + 2*m2) + l2*m2*sin(q0 + q1 + q2) + 2*l0*sin(fi0 + q0)*(m1 + m2))/(2*(m0 + m1 + m2))];
                baseCorners = [ baseCenter(1) - l0*cos(fi0+q0), baseCenter(2) - l0*sin(fi0+q0);...
                                baseCenter(1) + l0*cos(fi0-q0), baseCenter(2) - l0*sin(fi0-q0);...
                                baseCenter(1) + l0*cos(fi0+q0), baseCenter(2) + l0*sin(fi0+q0);...
                                baseCenter(1) - l0*cos(fi0-q0), baseCenter(2) + l0*sin(fi0-q0)];
                link1x = [ -(l1*cos(q0 + q1)*(m1 + 2*m2) + l2*m2*cos(q0 + q1 + q2) - 2*l0*m0*cos(fi0 + q0))/(2*(m0 + m1 + m2)), (l1*cos(q0 + q1)*(2*m0 + m1) - l2*m2*cos(q0 + q1 + q2) + 2*l0*m0*cos(fi0 + q0))/(2*(m0 + m1 + m2))];
                link1y = [ -(l1*sin(q0 + q1)*(m1 + 2*m2) + l2*m2*sin(q0 + q1 + q2) - 2*l0*m0*sin(fi0 + q0))/(2*(m0 + m1 + m2)), (l1*sin(q0 + q1)*(2*m0 + m1) - l2*m2*sin(q0 + q1 + q2) + 2*l0*m0*sin(fi0 + q0))/(2*(m0 + m1 + m2))];
                link2x = [ link1x(2), (l1*cos(q0 + q1)*(2*m0 + m1) + l2*cos(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*cos(fi0 + q0))/(2*(m0 + m1 + m2)) ];
                link2y = [ link1y(2), (l1*sin(q0 + q1)*(2*m0 + m1) + l2*sin(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*sin(fi0 + q0))/(2*(m0 + m1 + m2)) ];
            
                set(base,'Vertices', baseCorners);
                set(link1, 'XData', link1x, 'YData', link1y);
                set(link2, 'XData', link2x, 'YData', link2y);
                set(time,'String', num2str(timeVector(i)));
                drawnow
                %pause(0.0001);
            end
            
            
            
        end
        
        % calculates the angulat momentum of the system
        function K = angularMomentum(obj)
            m0 = obj.m0;
            m1 = obj.m1;
            m2 = obj.m2;
            l0 = obj.l0;
            l1 = obj.l1;
            l2 = obj.l2;
            I0 = obj.I0;
            I1 = obj.I1;
            I2 = obj.I2;
            fi0 = obj.fi0;
            q0 = obj.q(1);
            q1 = obj.q(2);
            q2 = obj.q(3);
            dq0 = obj.qd(1);
            dq1 = obj.qd(2);
            dq2 = obj.qd(3);

            K = (dq1*(4*I1*m0 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2)) + dq0*(4*I0*m0 + 4*I0*m1 + 4*I1*m0 + 4*I0*m2 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + 4*l0^2*m0*m1 + 4*l0^2*m0*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 4*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 4*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2)) + dq2*(4*I2*m0 + 4*I2*m1 + 4*I2*m2 + l2^2*m0*m2 + l2^2*m1*m2 + l1*l2*m2*cos(q2)*(2*m0 + m1) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2)))/(4*(m0 + m1 + m2));
        end
        
        % calculates the kinetic energy of the system
        function T = kineticEnergy(obj)
             m0 = obj.m0;
            m1 = obj.m1;
            m2 = obj.m2;
            l0 = obj.l0;
            l1 = obj.l1;
            l2 = obj.l2;
            I0 = obj.I0;
            I1 = obj.I1;
            I2 = obj.I2;
            fi0 = obj.fi0;
            q0 = obj.q(1);
            q1 = obj.q(2);
            q2 = obj.q(3);
            dq0 = obj.qd(1);
            dq1 = obj.qd(2);
            dq2 = obj.qd(3);

            T = (I1*(dq0 + dq1)^2)/2 + (I0*dq0^2)/2 + (I2*(dq0 + dq1 + dq2)^2)/2 + (m1*((l1*cos(q0 + q1)*(dq0 + dq1)*(m0 - m2) + 2*dq0*l0*m0*cos(fi0 + q0) - l2*m2*cos(q0 + q1 + q2)*(dq0 + dq1 + dq2))^2 + (l1*sin(q0 + q1)*(dq0 + dq1)*(m0 - m2) + 2*dq0*l0*m0*sin(fi0 + q0) - l2*m2*sin(q0 + q1 + q2)*(dq0 + dq1 + dq2))^2))/(8*(m0 + m1 + m2)^2) + (m2*((l1*cos(q0 + q1)*(dq0 + dq1)*(2*m0 + m1) + 2*dq0*l0*m0*cos(fi0 + q0) + l2*cos(q0 + q1 + q2)*(m0 + m1)*(dq0 + dq1 + dq2))^2 + (l1*sin(q0 + q1)*(dq0 + dq1)*(2*m0 + m1) + 2*dq0*l0*m0*sin(fi0 + q0) + l2*sin(q0 + q1 + q2)*(m0 + m1)*(dq0 + dq1 + dq2))^2))/(8*(m0 + m1 + m2)^2) + (m0*((2*dq0*l0*cos(fi0 + q0)*(m1 + m2) + l1*cos(q0 + q1)*(dq0 + dq1)*(m1 + 2*m2) + l2*m2*cos(q0 + q1 + q2)*(dq0 + dq1 + dq2))^2 + (2*dq0*l0*sin(fi0 + q0)*(m1 + m2) + l1*sin(q0 + q1)*(dq0 + dq1)*(m1 + 2*m2) + l2*m2*sin(q0 + q1 + q2)*(dq0 + dq1 + dq2))^2))/(8*(m0 + m1 + m2)^2);
        end
        
        function pos = endEffectorPos(obj)
            global step stateVariables
            m0 = obj.m0;
            m1 = obj.m1;
            m2 = obj.m2;
            l0 = obj.l0;
            l1 = obj.l1;
            l2 = obj.l2;
            I0 = obj.I0;
            I1 = obj.I1;
            I2 = obj.I2;
            fi0 = obj.fi0;
            q0 = stateVariables(1,step);
            q1 = stateVariables(2,step);
            q2 = stateVariables(3,step);
            
            pos = [ (l1*cos(q0 + q1)*(2*m0 + m1) + l2*cos(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*cos(fi0 + q0))/(2*(m0 + m1 + m2));
                    (l1*sin(q0 + q1)*(2*m0 + m1) + l2*sin(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*sin(fi0 + q0))/(2*(m0 + m1 + m2))];

        end
        
        function pos = calcEndEffectorPos(obj, q)
            m0 = obj.m0;
            m1 = obj.m1;
            m2 = obj.m2;
            l0 = obj.l0;
            l1 = obj.l1;
            l2 = obj.l2;
            I0 = obj.I0;
            I1 = obj.I1;
            I2 = obj.I2;
            fi0 = obj.fi0;
            q0 = q(1);
            q1 = q(2);
            q2 = q(3);
            
            pos = [ (l1*cos(q0 + q1)*(2*m0 + m1) + l2*cos(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*cos(fi0 + q0))/(2*(m0 + m1 + m2));
                    (l1*sin(q0 + q1)*(2*m0 + m1) + l2*sin(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*sin(fi0 + q0))/(2*(m0 + m1 + m2))];
        end
        
        % initializes symbolic variables
        function makeSymbolicVariables(obj)
            % state variables
            q_sym = sym('q', [3 1]);
            qd_sym = sym('qd', [3 1]);
            qdd_sym = sym('qdd', [3 1]);
            u_sym = sym('u', [2 1]);
            
            m0 = obj.m0;
            m1 = obj.m1;
            m2 = obj.m2;
            l0 = obj.l0;
            l1 = obj.l1;
            l2 = obj.l2;
            I0 = obj.I0;
            I1 = obj.I1;
            I2 = obj.I2;
            fi0 = obj.fi0;
            q0 = q_sym(1);
            q1 = q_sym(2);
            q2 = q_sym(3);
            dq0 = qd_sym(1);
            dq1 = qd_sym(2);
            dq2 = qd_sym(3);
            
            m11 = (4*I0*m0 + 4*I0*m1 + 4*I1*m0 + 4*I0*m2 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + 4*l0^2*m0*m1 + 4*l0^2*m0*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 4*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 4*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2))/(4*(m0 + m1 + m2));
            m12 = (4*I1*m0 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2))/(4*(m0 + m1 + m2));
            m13 = (4*I2*m0 + 4*I2*m1 + 4*I2*m2 + l2^2*m0*m2 + l2^2*m1*m2 + l1*l2*m2*cos(q2)*(2*m0 + m1) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2))/(4*(m0 + m1 + m2));
            m21 = (4*I1*m0 + 4*I1*m1 + 4*I2*m0 + 4*I1*m2 + 4*I2*m1 + 4*I2*m2 + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l0*l1*m0*cos(fi0 - q1)*(m1 + 2*m2) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2) + 4*l1*l2*m0*m2*cos(q2) + 2*l1*l2*m1*m2*cos(q2))/(4*(m0 + m1 + m2));
            m22 = (4*I2*m0 + 4*I2*m1 + 4*I2*m2 + 4*I1*(m0 + m1 + m2) + l1^2*m0*m1 + 4*l1^2*m0*m2 + l1^2*m1*m2 + l2^2*m0*m2 + l2^2*m1*m2 + 2*l1*l2*m2*cos(q2)*(2*m0 + m1))/(4*(m0 + m1 + m2));
            m23 = (m2*(m0 + m1)*l2^2 + l1*m2*cos(q2)*(2*m0 + m1)*l2 + 4*I2*(m0 + m1 + m2))/(4*(m0 + m1 + m2));
            m31 = (4*I2*m0 + 4*I2*m1 + 4*I2*m2 + l2^2*m0*m2 + l2^2*m1*m2 + l1*l2*m2*cos(q2)*(2*m0 + m1) + 2*l0*l2*m0*m2*cos(q1 - fi0 + q2))/(4*(m0 + m1 + m2));
            m32 = (m2*(m0 + m1)*l2^2 + l1*m2*cos(q2)*(2*m0 + m1)*l2 + 4*I2*(m0 + m1 + m2))/(4*(m0 + m1 + m2));
            m33 = (m2*(m0 + m1)*l2^2 + 4*I2*(m0 + m1 + m2))/(4*(m0 + m1 + m2));
            
            c1 = (dq0*(4*dq1*l0*m0*(l1*sin(fi0 - q1)*(m1 + 2*m2) - l2*m2*sin(q1 - fi0 + q2)) - 2*dq2*l2*m2*(l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2))) - dq2^2*l2*m2*(l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2)) + 2*dq1^2*l0*m0*(l1*sin(fi0 - q1)*(m1 + 2*m2) - l2*m2*sin(q1 - fi0 + q2)) - 2*dq1*dq2*l2*m2*(l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2)))/(4*(m0 + m1 + m2));
            c2 = -(2*l0*m0*(l1*sin(fi0 - q1)*(m1 + 2*m2) - l2*m2*sin(q1 - fi0 + q2))*dq0^2 + 2*dq2*l1*l2*m2*sin(q2)*(2*m0 + m1)*dq0 + dq2*l1*l2*m2*sin(q2)*(2*dq1 + dq2)*(2*m0 + m1))/(4*(m0 + m1 + m2));
            c3 = (l2*m2*((l1*sin(q2)*(2*m0 + m1) + 2*l0*m0*sin(q1 - fi0 + q2))*dq0^2 + 2*l1*sin(q2)*(2*m0 + m1)*dq0*dq1 + l1*sin(q2)*(2*m0 + m1)*dq1^2))/(4*(m0 + m1 + m2));

                        
            obj.M_sym = [m11, m12, m13;...
                         m21, m22, m23;...
                         m31, m32, m33];
                
            obj.C_sym = [c1;...
                         c2;...
                         c3];
         
            obj.H_sym = obj.H;
            
            obj.endEffector_sym = [ (l1*cos(q0 + q1)*(2*m0 + m1) + l2*cos(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*cos(fi0 + q0))/(2*(m0 + m1 + m2));
                                    (l1*sin(q0 + q1)*(2*m0 + m1) + l2*sin(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*sin(fi0 + q0))/(2*(m0 + m1 + m2))];
                                
            obj.endEffectorVel_sym = [ -(l1*sin(q0 + q1)*(dq0 + dq1)*(2*m0 + m1) + 2*dq0*l0*m0*sin(fi0 + q0) + l2*sin(q0 + q1 + q2)*(2*m0 + 2*m1 + m2)*(dq0 + dq1 + dq2))/(2*(m0 + m1 + m2));
                                        (l1*cos(q0 + q1)*(dq0 + dq1)*(2*m0 + m1) + 2*dq0*l0*m0*cos(fi0 + q0) + l2*cos(q0 + q1 + q2)*(2*m0 + 2*m1 + m2)*(dq0 + dq1 + dq2))/(2*(m0 + m1 + m2))];

            obj.q_sym = q_sym;
            obj.qd_sym = qd_sym;
            obj.qdd_sym = qdd_sym;
            obj.u_sym = u_sym;
            
        end
        
        function q = jointCoordinates(obj,x,y)
            m0 = obj.m0;
            m1 = obj.m1;
            m2 = obj.m2;
            l0 = obj.l0;
            l1 = obj.l1;
            l2 = obj.l2;
            I0 = obj.I0;
            I1 = obj.I1;
            I2 = obj.I2;
            fi0 = obj.fi0;
            q0 = obj.q(1);
            q1 = sym('q1');
            q2 = sym('q2');
            
            pos = [ (l1*cos(q0 + q1)*(2*m0 + m1) + l2*cos(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*cos(fi0 + q0))/(2*(m0 + m1 + m2));
                    (l1*sin(q0 + q1)*(2*m0 + m1) + l2*sin(q0 + q1 + q2)*(2*m0 + 2*m1 + m2) + 2*l0*m0*sin(fi0 + q0))/(2*(m0 + m1 + m2))];
            [s1,s2]=vpasolve([ pos(1) == x, pos(2) == y], [q1,q2]);
            
            if size(s1,2) > 1
                if abs(s1(1)-obj.q(2)) < abs(s1(2)-obj.q(2))
                    choice = 1;
                else
                    choice = 2;
                end
            else
                choice = 1;
            end
            
            q1sol = wrapTo2Pi(s1(choice));
            q2sol = wrapTo2Pi(s2(choice));
            
            q = [q0; q1sol; q2sol];           
        end
    end
    
end