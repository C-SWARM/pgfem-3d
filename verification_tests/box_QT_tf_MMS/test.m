% Author: Phillip Hughes
% Date: July 27, 2015

clc;
clear;

discretization9 = load('discretization9.000000e-01.in');
discretization7 = load('discretization7.000000e-01.in');
discretization5 = load('discretization5.000000e-01.in');
discretization3 = load('discretization3.000000e-01.in');
error1 = load('error3.txt');
error2 = load('error5.txt');
error3 = load('error7.txt');
error4 = load('error9.txt');

discretization = [discretization3, discretization5,...
    discretization7, discretization9];
errors = [error1(1,1), error2(1,1), error3(1,1), error4(1,1)];


% loglog(discretization, errors,'ro')
% hold on

p = polyfit(log10(discretization),log10(errors), 1);
% ymatrix = 10.^(p(1)*log10(discretization) + p(2));
% plot(discretization,ymatrix,'--r')
% xlabel('h (discretation size)')
% ylabel('Error')
% legend('Error',sprintf('slope = %d',p(1)), 'Location','SouthWest')
% hold off

slope = p(1);

if slope < 2.25 && slope > 1.75
	disp('second order convergence');
else
    error('Does not reach second order convergence');
end
