function dydt = neuron(t,y)
global phi A omega alpha delta mu
dydt=zeros(2,1);
dydt(1)=y(1)*(1-phi)-(1/3)*y(1).^3-y(2)+A*cos(omega*t);
dydt(2)=mu*(y(1)+alpha-delta*y(2));
