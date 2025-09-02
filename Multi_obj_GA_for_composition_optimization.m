close all
clear all
clc

FitnessFunction = @input_p; % Function handle to the fitness function
numberOfVariables = 2; % Number of decision variables
lb = [0.06,0.06]; % Lower bound
ub = [0.28,0.28]; % Upper bound
A = []; % No linear inequality constraints
b = []; % No linear inequality constraints
Aeq = []; % No linear equality constraints
beq = []; % No linear equality constraints
nonlcon= @co; %non linear constaints
options = optimoptions('gamultiobj','PopulationSize',50,'ParetoFraction', 0.35,'FunctionTolerance',1E-4);
[x,Fval,exitFlag,Output] = gamultiobj(FitnessFunction,numberOfVariables,A,b,Aeq,beq,lb,ub,nonlcon,options);

Omega=1./Fval(:,1); Delta=Fval(:,2), Density=Fval(:,3);
plot3(abs(Fval(:,1)),abs(Fval(:,2)),Fval(:,3),'*')
hold on
grid on
xlabel('Obj1')
ylabel('Obj2')
zlabel('Obj3')

delH=4.*(0.6.*x(:,1).*(29.51)+0.6.*x(:,2).*(-3.445)+0.6.*(0.4-x(:,1)-x(:,2)).*13.542+ x(:,1).*x(:,2).*(-9.194)+ x(:,1).*(0.4-x(:,1)-x(:,2)).*1.975+ x(:,2).*(0.4-x(:,1)-x(:,2)).*(-22.324));
delS=(-8.3144.*(0.6.*log(0.6)+x(:,1).*log(x(:,1))+x(:,2).*log(x(:,2))+(0.4-x(:,1)-x(:,2)).*log(0.4-x(:,1)-x(:,2))));
VEC=0.6.*2+x(:,1).*5+x(:,2).*12+(0.4-x(:,1)-x(:,2)).*4;
Lambda=delS(:,1)./(Fval(:,2)).^2;

table(x,Omega,Delta, Density, delH,delS,VEC,Lambda)