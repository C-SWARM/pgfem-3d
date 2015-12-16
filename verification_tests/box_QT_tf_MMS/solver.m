% Author: Phillip Hughes
% Date: July 23, 2015

function solver(timestep, time)

file1 = fopen('original_box_QT_0.in.st','r'); % reads file
file2 = fopen('box_QT_0.in.st','wt'); % writes to new file

tline = fgets(file1);
fprintf(file2,'%s\n',tline); % copies the second to the 8th line

for i = 1:3
    tline = fgets(file1);
end

fprintf(file2,'%d\n', timestep);

for j = 1:timestep+1
    A(j) = ((j-1)/timestep)*time;
end

fprintf(file2, '%d ', A);
fprintf(file2, '\n');
fprintf(file2, '\n');
fprintf(file2, '%i\n', 1);
fprintf(file2, '%i\n', timestep-1);
fprintf(file2, '\n');
fprintf(file2, '%i\n', 0);

