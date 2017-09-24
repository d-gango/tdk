clear all
close all

global timeStep stateVariables timeVector maxStep step desiredPosition
%initialization
timeStep = 0.01;
maxStep = 450;
stateVariables = zeros(6, maxStep+1);
timeVector = zeros(1, maxStep+1);
angularMomentum = zeros(1, maxStep+1);
kineticEnergy = zeros(1, maxStep+1);
uVector = zeros(2, maxStep+1);
endEffectorPosition = zeros(2, maxStep+1);
desiredPosition = zeros(2, maxStep+1);
errorNorm = zeros(1, maxStep+1);
step = 1;

model = robotModel([0 pi/3 -pi/2 0 0 0]);

stateVariables(:,1) = model.getStateVariables();
angularMomentum(1) = model.angularMomentum();
kineticEnergy(1) = model.kineticEnergy();
endEffectorPosition(:,1) = model.endEffectorPos();

controller = LaplaceController(model);
errorNorm(1) = norm(desiredPosition(:,1) - endEffectorPosition(:,1));
for i = 1:maxStep
    % calculate u
    u = controller.getU(model);
    model.u = u;
    % save the input
    uVector(:,step) = u;
    
    % step the simulation
    model.integrate();
    step = step+1;
    
    %save the variables
    timeVector(step) = (step-1)*timeStep;
    stateVariables(:,step) = model.getStateVariables();
    angularMomentum(step) = model.angularMomentum();
    kineticEnergy(step) = model.kineticEnergy();
    endEffectorPosition(:,step) = model.endEffectorPos();
    errorNorm(i) = norm(desiredPosition(:,i) - endEffectorPosition(:,i));
    
    disp(timeVector(step));
end

roundedEndEff = round(endEffectorPosition,4);
model.animate();

figure
plot(timeVector, angularMomentum)
title('angular momentum')

figure
plot(timeVector, kineticEnergy)
title('kinetic energy')

figure
plot(roundedEndEff(1,:), roundedEndEff(2,:), desiredPosition(1,:), desiredPosition(2,:))
title('end effector position')
legend('actual', 'desired')
grid on
axis equal

figure
plot(timeVector, errorNorm);
title('error')
