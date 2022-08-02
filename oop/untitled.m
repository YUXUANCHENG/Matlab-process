figure(1);hold on; box on;
set(gcf,'color','w');
%scatter([0.215,0.148],[0.820,1.198],140,'<','filled')
% scatter([0.0955,0.092,0.0306,0.0042,0.0014,0.000154],[1.579,1.620,1.766,2.254,2.418,2.815],140,'<','filled')
% scatter([0.049,0.085,0.147,0.213],[1.02,0.895,0.85,0.80],140,'<','filled')
%scatter([0.049,0.085,0.147,0.213],[0,1.01,0,0],140,'<','filled')

% scatter([0.1593],[0.994],140,'>','filled')

% % scatter([0.00516,0.012,0.0065,0.0013],[1.49,1.679,2.14,2.76],140,'>','filled')
scatter([1.1546,1.112,1.00007,1.0000013]-1,[1.49,1.679,2.14,2.76],140,'>','filled','DisplayName','calA0 = 1.15')
% 
scatter([1.0871,1.0473,1.0000434]-1,[1.643,1.778,2.236],140,'>','filled','DisplayName','calA0 = 1.08')
scatter([1.0399,1.009,1.0000324]-1,[1.747,1.8612,2.32],140,'>','filled','DisplayName','calA0 = 1.03')


% scatter([0.0761],[1.78],140,'>','filled')
% scatter([0.0855,0.0321,0.0042,0.0013],[1.6,1.6,2.48,2.48],140,'>','filled','DisplayName','calA0 = 1.0')
scatter([1.01628,1.00193]-1,[1.78,1.886],140,'>','filled','DisplayName','calA0 = 1.0')
scatter([1.01628,1.00193]-1,[1.81,1.92],140,'>','filled','DisplayName','calA0 = 1.0, same non-dim Ks')
% % 
% scatter(arrayfun(@trans2calA,[0.0002,0.0006,0.002,0.006,0.02,0.06,0.2])-1,[3.003,2.837,2.553,2.237,1.826,1.095,0.879],140,'o','filled','DisplayName','Soft particle simulation')
scatter([1.01675,1.01648]-1,[1.429,1.12],140,'s','filled','DisplayName','calA0 = 1.0, reduce friction')
scatter([1.00000119,1.00000112]-1,[2.67,2.6],140,'s','filled','DisplayName','calA0 = 1.15, reduce friction')
scatter([1.00000119,1.00000112]-1,[3.15,3.5],140,'p','filled','DisplayName','calA0 = 1.15, increase N')
scatter([1.00000119]-1,[2.37],140,'p','filled','DisplayName','calA0 = 1.15, frictionless')

% legend('calA0 = 1.15','calA0 = 1.08','calA0 = 1.03','calA0 = 1.0','calA0 = 1.0, same non-dim Ks','calA0 = 1.0, reduce friction','calA0 = 1.15, reduce friction')
% legend('calA0 = 1.15','calA0 = 1.08','calA0 = 1.03','calA0 = 1.0')

% scatter([0.006090,0.004718,0.002978,0.0018956],[1.76,1.92,2.10,2.37],140,'>','filled')
% scatter([0.0081,0.0067,0.0056],[1.97,2.06,2.18],140,'v','filled')
% scatter(arrayfun(@trans2calA,[0.0081,0.0067,0.0056])-1,[1.97,2.06,2.18],140,'v','filled')
% scatter([0.0125,0.0104,0.0061],[1.45,1.82,2.06],140,'^','filled')
% scatter(arrayfun(@trans2calA,[0.0125,0.0104,0.0061])-1,[1.45,1.82,2.06],140,'^','filled')
% % scatter([0.0006876,0.00057567],[2.99,3.04],140,'s','filled')
% % scatter([0.0000018552,0.0000006897],[3.06,3.36],140,'d','filled')
% scatter(arrayfun(@trans2calA,[0.0000018552,0.0000006897])-1,[3.06,3.36],140,'d','filled')
% % scatter([0.0002,0.0006,0.002,0.006,0.02,0.06,0.2],[3.003,2.837,2.553,2.237,1.826,1.095,0.879],140,'o','filled')
scatter(arrayfun(@trans2calA,[0.0002,0.0006,0.002,0.006,0.02,0.06,0.2])-1,[3.003,2.837,2.553,2.237,1.826,1.095,0.879],140,'o','filled','DisplayName','Soft particle simulation')
legend
xline(1.15-1,'--r','DisplayName','calA = 1.15')
ax = gca;
ax.XScale = "log";
%ax.YScale = "log";
%xlabel('delta/d');
xlabel('calA - 1');
ylabel('a');

% 0,8E-6 5,1e-10 9,1e-12
%%
figure(2);hold on; box on;
set(gcf,'color','w');
plot(linspace(1.1,1.22,5),[0.003,0.0027,0.0021,0.0016,0.0012].^2);
plot(linspace(1,1.22,8),[0.0004,0.0008,0.001,0.0014,0.0022,0.0016,0.0016,0.001].^2);
plot(linspace(1,1.18,8),[0.0006,0.0006,0.0004,0.0004,0.0004,0.0004,0.0004,0.0004].^2);
plot(linspace(1,1.18,8),[0.0006,0.0006,0.0006,0.0008,0.0008,0.0006,0.0006,0.0006].^2);
% ylim([0 0.003])
xlabel('calA');
ylabel('Teff');

figure(3);hold on; box on;
set(gcf,'color','w');
scatter([0.838,0.921,0.945,0.958,0.964,0.966,0.965,0.967],[0.0004,0.0008,0.001,0.0014,0.0022,0.0016,0.0016,0.001].^2);
scatter([0.828,0.872,0.8748,0.8749,0.8693,0.8678,0.8670,0.8681],[0.0006,0.0006,0.0004,0.0004,0.0004,0.0004,0.0004,0.0004].^2);
scatter([0.823,0.830,0.837,0.8378,0.841,0.834,0.837,0.834],[0.0006,0.0006,0.0006,0.0008,0.0008,0.0006,0.0006,0.0006].^2);
% ylim([0 0.003])
xlabel('phiJ');
ylabel('Teff');

figure(4);hold on; box on;
set(gcf,'color','w');
plot(linspace(1,1.18,8),[0.02515,0.02656,0.02918,0.03188,0.03483,0.03772,0.04054,0.0433].^2);
xlabel('calA');
ylabel('U');
figure(5);hold on; box on;
set(gcf,'color','w');
plot(linspace(1,1.18,8),[8e-6,6e-7,3e-8,3e-9,1e-10,7e-11,3e-11,1e-12]);
ax = gca;
ax.YScale = "log";
xlabel('calA');
ylabel('U');

% function calA = trans2calA(delta)
%     angle = acos(1 - 2 * delta);
%     calA = (sin(angle) + (pi - angle))^2/((pi - angle)*pi);
% end

%%
figure(5);hold on; box on;
set(gcf,'color','w');
m= 1 * pi /4;
scatter([16, 32, 64], [m, m, m].*[1.37, 1.22, 0.91]/1,140,'d','filled','DisplayName','Fn = G')
scatter([16, 32, 64], [m, m, m].*[5.05, 4.2, 1.98]/5,140,'>','filled','DisplayName','Fn = 5G')
legend;
xlabel('NV');
ylabel('u');

%%
function calA = trans2calA(delta)
    b = 1 - delta;
    a = 1/b;
    per = ellipsePerimeter(a,b);
    calA = per^2/(4*pi^2);
end

function per = ellipsePerimeter(a,b,varargin)
% ELLIPSEPERIMETER Compute the perimeter of an ellipse.
%   per = ELLIPSEPERIMETER(a,b) Computes the perimeter PER of an ellipse 
%   with semi-axes A and B with an infinite series solution set up until 
%   order 20. 
%   per = ELLIPSEPERIMETER(a,b,order) Computes the perimter PER with
%   semi-axes A and B with an infinite series solution set up until the
%   ORDER selected.
%
%   ######################################################################
%   EXAMPLES:
%   >> per = ellipsePerimeter(1,1);
%   Gives as a result: per = 6.283185307179586 = 2*pi
%

p = inputParser;
validationFcn = @(x) isnumeric(x) && sum(x > 0)/size(x,2) == 1;
addRequired(p,'a',validationFcn);
addRequired(p,'b',validationFcn);
defaultOrderValue = 20;
addOptional(p,'order',defaultOrderValue,validationFcn);
parse(p,a,b,varargin{:});
order = p.Results.order;
% Program start
h = (a-b).^2./(a+b).^2;
s = 1;
for i=1:order
    s = s + (h.^i * (gamma(0.5+1)/gamma(0.5-i+1)/gamma(i+1))^2);
end
per = pi*(a+b).*s;

end
