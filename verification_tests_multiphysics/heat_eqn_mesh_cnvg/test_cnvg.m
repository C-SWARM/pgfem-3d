%% clean up
clear; clc;

out_dir = 'out';

filebase = 'box';
NP = 16;

tests = [50,75,100,125,150];   % test ids
testno = numel(tests); % number of tests
stepno = 9;
t = 0.001;

X  = cell(testno, NP);   % x,y,z coordinate from mesh
T  = cell(testno, NP);   % temperature
Ta = cell(testno, NP);   % temperature (analytic)
e  = zeros(testno, 1);   % error at each node
n  = zeros(testno, 1);   % total number of node for unknown


% pre-computed mesh size h = [average minmum]
h = [3.150330e-02 2.005973e-02  %  50
     4.694265e-02 2.872470e-02  %  75
     6.312517e-02 4.335686e-02  % 100
     7.821258e-02 4.724483e-02  % 125
     9.205791e-02 5.870671e-02];% 150

for ia=1: testno
  test_case = tests(ia);
  %% read mesh
  in_dir  = sprintf('%s_%dCPU_%.3d', filebase, NP, test_case);
  err = 0;
  for np=0: NP-1
    f_mesh = sprintf('%s/%s_%d.in', in_dir, filebase, np);
    fp = fopen(f_mesh);
    mesh_size = fscanf(fp, '%d', 3); % mesh_size(1) = nodeno
                                     % mesh_size(2) = nsd
                                     % mesh_size(3) = elemno
    % read values not used here
    fscanf(fp, '%e', 6 + mesh_size(3));
    
    % read coordinates
    x = fscanf(fp, '%e', 9*mesh_size(1));
    m = reshape(x, [9, mesh_size(1)])';

    id = m(:, 2) ~= np;

    nbc = fscanf(fp, '%d', 1);
    if(nbc>0)
      BC  = fscanf(fp, '%d', 2*nbc);
      bc  = reshape(BC, [2, nbc])';
    
      id(bc(:, 1)+1) = 1;
    end

    m(id, :) = [];
    X{ia, np+1} = m(:, 4:6);
    fclose(fp);

    % read temperature
    f_data = sprintf('%s/mesh_cnvg/%.3d/%s_%dCPU/restart/heat_transfer/STEP_%.6d/%s_%d_%d.res', out_dir,test_case,filebase,NP,stepno,filebase,np,stepno);
    T{ia, np+1} = load(f_data);
    T{ia, np+1}(id, :) = [];

    Ta{ia, np+1} = exp(-2*pi*pi*t)*sin(pi*m(:, 4)).*sin(pi*m(:, 5));
    dT = Ta{ia, np+1} - T{ia, np+1}(:, 2);
    nodeno = size(m, 1);
    n(ia) = n(ia) + nodeno;
    err = err + dT'*dT;
  end
  e(ia) = sqrt(err/n(ia));
end

% check convergence order
x = log10(h(:, 1));
y = log10(e);

A = [x, ones(size(x))];
b = y;
c = inv(A'*A)*(A'*b);
tol = 0.2;
if(2-tol<c(1) && c(1)<2+tol)
  fprintf('Test is successful: second order convergence is achieved.\n');
else
   error('Does not reach second order convergence.');
end


