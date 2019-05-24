%% loading given data
Q11 = csvread('Q11.csv',1,0);
Q21 = csvread('Q21.csv',1,0);
%% storing aerodynamic data
k_dat  = [Q11(:,1) Q21(:,1)];
Qr_dat = [Q11(:,2) Q21(:,2)];
Qi_dat = [Q11(:,3) Q21(:,3)];
Q_dat = [(Qr_dat(:,1)+i*Qi_dat(:,1)) (Qr_dat(:,2)+i*Qi_dat(:,2))];

n1 = length(k_dat); % num of row
n2 = length(k_dat(1,:)); % num of col
%% constants
beta1 = 0.2*ones(n1,1);
beta2 = 0.2*ones(n1,1);
%% start program
Q_tild = zeros(n1,1); % approximated solution

for index = 1:n2
  b1 = ones(n1,1);
  b2 = i*k_dat(:,index);
  b3 = -k_dat(:,index).^2;
  b4 = (i*beta1.*k_dat(:,index) + k_dat(:,index).^2)./(k_dat(:,index).^2 + beta1.^2);
  b5 = (i*beta2.*k_dat(:,index) + k_dat(:,index).^2)./(k_dat(:,index).^2 + beta2.^2);
  % find conjugate complex of B components
  b1_ = conj(b1); b2_ = conj(b2); b3_ = conj(b3); b4_ = conj(b4); b5_ = conj(b5);
  
  Bconj(:,1) = [b1_.'; b2_.'; b3_.'; b4_.'; b5_.']*b1;
  Bconj(:,2) = [b1_.'; b2_.'; b3_.'; b4_.'; b5_.']*b2;
  Bconj(:,3) = [b1_.'; b2_.'; b3_.'; b4_.'; b5_.']*b3;
  Bconj(:,4) = [b1_.'; b2_.'; b3_.'; b4_.'; b5_.']*b4;
  Bconj(:,5) = [b1_.'; b2_.'; b3_.'; b4_.'; b5_.']*b5; 
  rhs = [b1_.'; b2_.'; b3_.'; b4_.'; b5_.']*Q_dat(:,index);

  C = Bconj\rhs;
  Q_tild = [b1 b2 b3 b4 b5]*C;
  err_r = abs(real(Q_dat(:,index)) - real(Q_tild));
  err_i = abs(imag(Q_dat(:,index)) - imag(Q_tild));
  err  = sum(err_r + err_i);

  % plotting
  ... % more plotting options
  figure; % generalized aerodynamic forces
  plot(Qr_dat(:,index),Qi_dat(:,index)); grid on; hold on;
  plot(real(Q_tild), imag(Q_tild), '+');
  ... % more plotting options
  set(gca, 'FontSize', 12, 'FontName', 'TlwgTypewriter');
  ... % more plotting options

  ... % more plotting options
  set(get(gca, 'XLabel'), 'String', 'Q_{Re}');
  set(get(gca, 'YLabel'), 'String', 'Q_{Im}');
  legend('Q_{data}','Q_{approx}', 'Orientation', 'Vertical', 'Location', 'northeast');
end

