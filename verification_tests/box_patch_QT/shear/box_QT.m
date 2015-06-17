% Author: Phillip Hughes
% Date: June 15, 2015
% Title: Shear
    
clc;
clear;

files = dir; % determines the files in the directory

for id = 2:length(files) - 5 % determines the number of files in the directory...
    % minus the files: box_QT_0.in.st, box_QT.out.header, box_QT.t3d,...
    % makeset.pl, and this file
    k = 2^(id - 1);
    if isdir(sprintf('box_QT_%iCPU',k)) == 1
        cd (sprintf('box_QT_%iCPU',k)) % changes to the correct directory if found
        for i = 1:k  
            number = i - 1;
            file1 = fopen(sprintf('box_QT_%i.in',number),'r'); % reads file
            file2 = fopen(sprintf('box_QT%i.in',number),'wt'); % writes to new file
            %% initial values
            tline = fgets(file1);
            first_line = str2num(tline);
            fprintf(file2, '%d %d %d\n',first_line);
            nodes = first_line(1); % determines the number of nodes in input file
            element_number = first_line(3); % determines the element number 
            displacement = .1; % the amount of displacement of the sheared object

            for n = 1:7
                tline = fgets(file1);
                fprintf(file2,'%s',tline); % copies the second to the 8th line 
            end

            B = []; % refreshes the matrix every iteration    
            for j = 1:nodes
                tline = fgets(file1);
                A = str2num(tline); 
                if isempty(A) == 1        
                else
                    for k = 1:9
                        B(j,k) = A(k); % remembers many numbers including the x,y,and z coordinate
                    end
                end
            end
    
            fprintf(file2, '%4d %4d %4d %17.8e %16.8e %16.8e %2d %4d %1d\n', B');
            fprintf(file2,'\n');

            count = 0; 
            for l = 1:nodes
                if B(l,4) == 0 || B(l,5) == 0 || B(l,6) == 0 || B(l,4) == 1 ||...
                    B(l,5) == 1 || B(l,6) == 1
                    shear(1,l) = B(l,5)*displacement; % determines the amount of shear in the y-direction
                    count = count + 1; % iteration counter
                else
                end
            end

            fprintf(file2, '%d\n', count); % number of supported entities

            for m = 1:count
                bc = [m-1, -m, 1, 1]; % [geom id] [x comp. supporting type] [y comp. supporting type] [z comp. supporting type]
                fprintf(file2, '%d %d %d %d\n', bc);
            end
    
            fprintf(file2, '\n');
            fprintf(file2, '%d\n', count); % [# of prescribed]
            fprintf(file2, '%e  ', shear); % [value 1] [value 2] [value #]
            fprintf(file2, '\n');
            fprintf(file2, '\n');

            tline = fgets(file1);
            tline = fgets(file1);
            original_count = str2double(tline); 
            count2 = 0;
            while ischar(tline)
                tline = fgets(file1);
                if count2 > original_count + 3 % skips the values that have already been replaced
                    fprintf(file2,'%s',tline); % prints out the rest of the values
                end
                count2 = count2 + 1; % iteration counter
            end
           delete(sprintf('box_QT_%i.in',number)); % deletes the original *.in file
           movefile(sprintf('box_QT%i.in',number),sprintf('box_QT_%i.in',number)); % file being written to renamed to box_QT_...
        end
    else
    end
    cd .. % goes back up to original directory
end