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
         P = [70;50]
         D = [30;15]%[15;10]
    end
    
    methods
        function obj = LaplaceController(system)
            global desiredPosition maxStep timeStep laplaceSolution
            laplaceSolution = zeros(5, maxStep);
            %desired trajectory
            obj.t = sym('t');
            obj.t0 = sym('t0');
            
            obj.trajectory = [2 + 0.05*(obj.t + obj.t0)*cos(obj.t + obj.t0);...
                              1 + 0.05*(obj.t + obj.t0)*sin(obj.t + obj.t0)];
            for i = 1:maxStep+1
                desiredPosition(:,i) = double(subs(obj.trajectory, [obj.t, obj.t0], [(i-1)*timeStep, 0]));
            end
            
            %mozgásegyenlet és szervo-kényszer
            eq1 = system.M_sym * system.qdd_sym + system.C_sym - system.H *system.u_sym;
            eq2 = system.endEffector_sym - obj.trajectory;
            
            eq = vertcat(eq1, eq2);
            
            %kezdeti értékek változói
            obj.q0 = sym('q0', [3 1]);
            obj.qd0 = sym('qd0', [3 1]);
            obj.qdd0 = sym('qdd0', [3 1]);
            
            %linearizálás
            F0 = subs(eq, [system.q_sym, system.qd_sym, system.qdd_sym], [obj.q0, obj.qd0, obj.qdd0]);
            Fq = subs(jacobian(eq, system.q_sym), [system.q_sym, system.qd_sym, system.qdd_sym], [obj.q0, obj.qd0, obj.qdd0]) * (system.q_sym - obj.q0);
            Fdq = subs(jacobian(eq, system.qd_sym), [system.q_sym, system.qd_sym, system.qdd_sym], [obj.q0, obj.qd0, obj.qdd0]) * (system.qd_sym - obj.qd0);
            Fddq = subs(jacobian(eq, system.qdd_sym), [system.q_sym, system.qd_sym, system.qdd_sym], [obj.q0, obj.qd0, obj.qdd0]) * (system.qdd_sym - obj.qdd0);
            
            eqlin = F0 + Fq + Fdq + Fddq;
            
            %szimbolikus idõfüggõ változók és deriváltjaik,  u1 és u2 KONSTANS
            syms q0(t) q1(t) q2(t) u1 u2
            q = [q0; q1; q2];
            dq = diff(q, obj.t);
            ddq = diff(q, obj.t, 2);
  %         obj.xd = sym('xd(t)');
  
            %idõfüggõ változók behelyettesítése
            eqlin_t = subs(eqlin, [system.q_sym, system.qd_sym, system.qdd_sym], [q, dq, ddq]);
            
            %Laplace-transzformáció
            obj.s = sym('s');
            EQ = laplace(eqlin_t, obj.t, obj.s);
            
            %kezdeti értékek rendbehozása
            EQinit = subs(EQ, 'q0(0)', obj.q0(1));
            EQinit = subs(EQinit, 'q1(0)', obj.q0(2));
            EQinit = subs(EQinit, 'q2(0)', obj.q0(3));
            EQinit = subs(EQinit, 'D(q0)(0)', obj.qd0(1));
            EQinit = subs(EQinit, 'D(q1)(0)', obj.qd0(2));
            EQinit = subs(EQinit, 'D(q2)(0)', obj.qd0(3));
            
            %Laplace-traszformált függvények behelyettesítése
            obj.Q0 = sym('Q0');
            obj.Q1 = sym('Q1');
            obj.Q2 = sym('Q2');
            
            EQsub = subs(EQinit, 'laplace(q0(t), t, s)', obj.Q0);
            EQsub = subs(EQsub, 'laplace(q1(t), t, s)', obj.Q1);
            EQsub = subs(EQsub, 'laplace(q2(t), t, s)', obj.Q2);
            
            %elõírt pálya behelyettesítése
            %EQsub = subs(EQsub, 'laplace(xd(t), t, s)', laplace(obj.trajectory, t, s))
            

            obj.EQfinal = collect(EQsub, [obj.Q0, obj.Q1, obj.Q2, system.u_sym(1), system.u_sym(2)]);
             
        end
        
        function u_ret = getU(obj, system)
            global timeStep stateVariables timeVector step laplaceSolution
            %kezdeti értékek mentése
            state = stateVariables(:, step);
            deriv = system.derivatives(state);
         
            %kezdeti értékek behelyettesítése
            eq = subs(obj.EQfinal, [obj.q0; obj.qd0; obj.qdd0; obj.t0], [state(1:3); state(4:6); deriv(4:6); timeVector(step)]);
            eq = collect(eq,[obj.Q0, obj.Q1, obj.Q2, system.u_sym(1), system.u_sym(2)]);
            
            %egyenletrendszer megoldása
            [Q0sol, Q1sol, Q2sol, U1sol, U2sol] = vpasolve(eq,[obj.Q0, obj.Q1, obj.Q2, system.u_sym(1), system.u_sym(2)]);
            
            solution = [ilaplace(Q0sol, obj.s, obj.t);...
                       ilaplace(Q1sol, obj.s, obj.t);...
                       ilaplace(Q2sol, obj.s, obj.t);...
                       ilaplace(U1sol, obj.s, obj.t);...
                       ilaplace(U2sol, obj.s, obj.t)];
            sol = real(double(limit(solution, obj.t, 0, 'right')));
            sol(1:3) = wrapTo2Pi(sol(1:3));
            laplaceSolution(:,step) = sol;
            
            %hibaszámítás a PD-hez
           % disp(state(1));
            error = [sol(2) - state(2);...
                     sol(3) - state(3)];
            if abs(error(1)) > pi
                error(1) = wrapToPi(sol(2)) - wrapToPi(state(2));
            end
            if abs(error(2)) > pi
                error(2) = wrapToPi(sol(3)) - wrapToPi(state(3));
            end
            
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

