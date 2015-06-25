% Author: Phillip Hughes
% Date: June 15, 2015
% Title: Shear
    
clc;
clear;

filename = 'box_LT';
files = dir; % determines the number of files in the directory

for id = 1:length(files) % determines the number of files in the directory...
    k = 2^(id - 1);

    if isdir(sprintf(strcat(filename,'_%iCPU'),k)) == 1
        cd (sprintf(strcat(filename,'_%iCPU'),k)) % changes to the correct directory if found
        for i = 1:k  
            number = i - 1;
            file1 = fopen(sprintf(strcat(filename,'_%i.in'),number),'r'); % reads file
            file2 = fopen(sprintf(strcat(filename,'%i.in'),number),'wt'); % writes to new file
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
    
            fprintf(file2,'%4d %4d %4d   %4.8e   %4.8e   %4.8e %2d %4d %d\n', B');
            fprintf(file2,'\n');

            count = 0;

	    shear=[];
	    bcmat = [];	 
            for l = 1:nodes
                if B(l,4) == 0 || B(l,5) == 0 || B(l,6) == 0 || B(l,4) == 1 ||...
                    B(l,5) == 1 || B(l,6) == 1
			count = count + 1; % iteration counter 
                   	bc = [B(l,3) -count 1 1]; % [geom id] [x comp. supporting type] [y comp. supporting type] [z comp. supporting type]
		      %  tline = fgets(file1);
			for o = 1:4
			    bcmat(count,o) = bc(o); % matrix that remembers the values in bc
			end
			shear(1,count) = B(l,5)*displacement; % determines the amount of shear in the y-direction
                 else
		%	    bcmat(l,1) = l - 1;
		%	tline = fgets(file1);
		%	tlinemat = str2num(tline)
		% 	for o = 1:4
		%	    bcmat(l,o) = tlinemat(o)
		%	end
                end
            end

            fprintf(file2, '%d\n', count); % number of supported entities
            fprintf(file2, '%4d  %d %d %d\n', bcmat');
               
            fprintf(file2, '\n');
            fprintf(file2, '%d\n', count); % [# of prescribed]
            fprintf(file2, '%e  ', shear); % [value 1] [value 2] [value #]
            fprintf(file2, '\n');
            fprintf(file2, '\n');

            tline = fgets(file1);
            tline = fgets(file1);
            count2 = 0;
	   original_count = str2double(tline);
            while ischar(tline)
                tline = fgets(file1);
                if count2 > original_count + 3 % skips the values that have already been replaced
                    fprintf(file2,'%s',tline); % prints out the rest of the values
                end
                count2 = count2 + 1; % iteration counter
            end
           delete(sprintf(strcat(filename,'_%i.in'),number)); % deletes the original *.in file
           movefile(sprintf(strcat(filename,'%i.in'),number),sprintf(strcat(filename,'_%i.in'),number)); % file being written to renamed to filename...
        end
    cd .. % goes back up to original directory
    else
    end
end
