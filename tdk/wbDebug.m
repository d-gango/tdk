q_error = wrapToPi(controller.q_des- stateVariables(1:3,1:end-1));

figure
plot( timeVector(1:end-1), q_error(1,:), timeVector(1:end-1), q_error(2,:),timeVector(1:end-1), q_error(3,:));
title('q errors')
legend('q0','q1','q2')