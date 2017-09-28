model = robotModel([0 0 0 0 0 0]);

xy_laplace = zeros(2,size(laplaceSolution,2));
for i = 1:size(laplaceSolution,2)
    xy_laplace(:,i) = model.calcEndEffectorPos(laplaceSolution(1:3,i));
end

error0 = laplaceSolution(1,:) - stateVariables(1,1:end-1);
error1 = laplaceSolution(2,:) - stateVariables(2,1:end-1);
error2 = laplaceSolution(3,:) - stateVariables(3,1:end-1);
u_pd = uVector - laplaceSolution(4:5,:);

figure
plot(xy_laplace(1,:), xy_laplace(2,:), desiredPosition(1,:), desiredPosition(2,:))
title('end effector position')
legend('solution', 'desired')
grid on
axis equal

figure
plot(timeVector(1:end-1), error0)
title('q0 error')
grid on

figure
plot(timeVector(1:end-1), error1)
title('q1 error')
grid on

figure
plot(timeVector(1:end-1), error2)
title('q2 error')
grid on

figure
plot(timeVector(1:end-1), laplaceSolution(4,:))
title('u1 laplace')
grid on


figure
plot(timeVector(1:end-1), laplaceSolution(5,:))
title('u2 laplace')
grid on

figure
plot(timeVector(1:end-1), u_pd(1,:))
title('u1 PD')
grid on

figure
plot(timeVector(1:end-1), u_pd(2,:))
title('u2 PD')
grid on