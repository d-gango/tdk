classdef LaplaceController < handle
    
    
    properties (SetAccess = private, GetAccess = public)

         t
         s
         t0
         q0
         qd0
         qdd0
         Q0
         Q1
         Q2
         trajectory
         EQfinal
         prevError
         P = [50;50]
         D = [20;10]
    end
    
    methods
        function obj = LaplaceController(system)
            global desiredPosition maxStep timeStep laplaceSolution
            laplaceSolution = zeros(5, maxStep);
            %desired trajectory
            obj.t = sym('t');
            obj.t0 = sym('t0');
            
            obj.trajectory = [2.0222 + 0.05*(obj.t+obj.t0) * cos((obj.t+obj.t0));...
                              0.7154 + 0.05*(obj.t+obj.t0) * sin((obj.t+obj.t0))];
            for i = 1:maxStep+1
                desiredPosition(:,i) = double(subs(obj.trajectory, [obj.t, obj.t0], [(i-1)*timeStep, 0]));
            end
            
            %mozg�segyenlet �s szervo-k�nyszer
            eq1 = system.M_sym * system.qdd_sym + system.C_sym - system.H *system.u_sym;
            eq2 = system.endEffector_sym - obj.trajectory;
            
            eq = vertcat(eq1, eq2);
            
            %kezdeti �rt�kek v�ltoz�i
            obj.q0 = sym('q0', [3 1]);
            obj.qd0 = sym('qd0', [3 1]);
            obj.qdd0 = sym('qdd0', [3 1]);
            
            %lineariz�l�s
            F0 = subs(eq, [system.q_sym, system.qd_sym, system.qdd_sym], [obj.q0, obj.qd0, obj.qdd0]);
            Fq = subs(jacobian(eq, system.q_sym), [system.q_sym, system.qd_sym, system.qdd_sym], [obj.q0, obj.qd0, obj.qdd0]) * (system.q_sym - obj.q0);
            Fdq = subs(jacobian(eq, system.qd_sym), [system.q_sym, system.qd_sym, system.qdd_sym], [obj.q0, obj.qd0, obj.qdd0]) * (system.qd_sym - obj.qd0);
            Fddq = subs(jacobian(eq, system.qdd_sym), [system.q_sym, system.qd_sym, system.qdd_sym], [obj.q0, obj.qd0, obj.qdd0]) * (system.qdd_sym - obj.qdd0);
            
            eqlin = F0 + Fq + Fdq + Fddq;
            
            %szimbolikus id�f�gg� v�ltoz�k �s deriv�ltjaik,  u1 �s u2 KONSTANS
            syms q0(t) q1(t) q2(t) u1 u2
            q = [q0; q1; q2];
            dq = diff(q, obj.t);
            ddq = diff(q, obj.t, 2);
  %         obj.xd = sym('xd(t)');
  
            %id�f�gg� v�ltoz�k behelyettes�t�se
            eqlin_t = subs(eqlin, [system.q_sym, system.qd_sym, system.qdd_sym], [q, dq, ddq]);
            
            %Laplace-transzform�ci�
            obj.s = sym('s');
            EQ = laplace(eqlin_t, obj.t, obj.s);
            
            %kezdeti �rt�kek rendbehoz�sa
            EQinit = subs(EQ, 'q0(0)', obj.q0(1));
            EQinit = subs(EQinit, 'q1(0)', obj.q0(2));
            EQinit = subs(EQinit, 'q2(0)', obj.q0(3));
            EQinit = subs(EQinit, 'D(q0)(0)', obj.qd0(1));
            EQinit = subs(EQinit, 'D(q1)(0)', obj.qd0(2));
            EQinit = subs(EQinit, 'D(q2)(0)', obj.qd0(3));
            
            %Laplace-traszform�lt f�ggv�nyek behelyettes�t�se
            obj.Q0 = sym('Q0');
            obj.Q1 = sym('Q1');
            obj.Q2 = sym('Q2');
            
            EQsub = subs(EQinit, 'laplace(q0(t), t, s)', obj.Q0);
            EQsub = subs(EQsub, 'laplace(q1(t), t, s)', obj.Q1);
            EQsub = subs(EQsub, 'laplace(q2(t), t, s)', obj.Q2);
            
            %el��rt p�lya behelyettes�t�se
            %EQsub = subs(EQsub, 'laplace(xd(t), t, s)', laplace(obj.trajectory, t, s))
            

            obj.EQfinal = collect(EQsub, [obj.Q0, obj.Q1, obj.Q2, system.u_sym(1), system.u_sym(2)]);
             
        end
        
        function u_ret = getU(obj, system)
            global timeStep stateVariables timeVector step laplaceSolution
            %kezdeti �rt�kek ment�se
            state = stateVariables(:, step);
            deriv = system.derivatives(state);
         
            %kezdeti �rt�kek behelyettes�t�se
            eq = subs(obj.EQfinal, [obj.q0; obj.qd0; obj.qdd0; obj.t0], [state(1:3); state(4:6); deriv(4:6); timeVector(step)]);
            eq = collect(eq,[obj.Q0, obj.Q1, obj.Q2, system.u_sym(1), system.u_sym(2)]);
            
            %egyenletrendszer megold�sa
            [Q0sol, Q1sol, Q2sol, U1sol, U2sol] = vpasolve(eq,[obj.Q0, obj.Q1, obj.Q2, system.u_sym(1), system.u_sym(2)]);
            
            solution = [ilaplace(Q0sol, obj.s, obj.t);...
                       ilaplace(Q1sol, obj.s, obj.t);...
                       ilaplace(Q2sol, obj.s, obj.t);...
                       ilaplace(U1sol, obj.s, obj.t);...
                       ilaplace(U2sol, obj.s, obj.t)];
            sol = real(double(limit(solution, obj.t, 0, 'right')));
            laplaceSolution(:,step) = sol;
            
            %hibasz�m�t�s a PD-hez
           % disp(state(1));
            error = [sol(2) - state(2);...
                     sol(3) - state(3)];
          %  disp(error);
            if step == 1
                errorDot = [0; 0];
                obj.prevError = error;
            else 
                errorDot = (error - obj.prevError) / timeStep;
                obj.prevError = error;
            end
                   
            u_ret = sol(4:5) + obj.P  .* error + obj.D .* errorDot;
        end            
    end
    
end

