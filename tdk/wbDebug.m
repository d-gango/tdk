figure
plot( timeVector(1:end-1), controller.q_error(1,:), timeVector(1:end-1), controller.q_error(2,:));
title('q errors')
legend('q1','q2')

figure
plot( timeVector(1:end-1), controller.q_des(1,:), timeVector(1:end-1), controller.q_des(2,:), timeVector(1:end-1), controller.q_des(3,:));
title('q_d')
legend('q0','q1','q2')

figure
plot( timeVector(1:end), stateVariables(1,:), timeVector(1:end), stateVariables(2,:), timeVector(1:end), stateVariables(3,:));
title('q')
legend('q0','q1','q2')