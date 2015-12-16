% Author: Phillip Hughes
% Date: June 29, 2015
% Title: Radius

function radius(makeset)
    
filename = 'box_QT';
files = dir; % determines the number of files in the directory

k = 2;

    if isdir(sprintf(strcat(filename,'_%iCPU'),k)) == 1
        file2 = fopen(sprintf(strcat('discretization','%d.in'),makeset),'wt'); % writes to new file
        cd (sprintf(strcat(filename,'_%iCPU'),k)) % changes to the correct directory if found
        for counter = 1:k  
            number = counter - 1;
            file1 = fopen(sprintf(strcat(filename,'_%i.in'),number),'r'); % reads file
            %% initial values
            tline = fgets(file1);
            first_line = str2num(tline);
            nodes = first_line(1); % determines the number of nodes in input file
            element_number = first_line(3); % determines the element number 

            for n = 1:7
                tline = fgets(file1);
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
            
            x = B(:,4);
            y = B(:,5);
            z = B(:,6);
            
            f = [x, y, z];
    
%           scatter3(f(:,1),f(:,2),f(:,3))



            tline = fgets(file1);
            tline = fgets(file1);

            original_count = str2double(tline);

            for y = 1:original_count + 4
                tline = fgets(file1);
            end
  
            AB = [];
            for id = 1:element_number
               tline = fgets(file1);
                a = str2num(tline);
                for k = 1:length(a)
                    AB(id,k) = a(k);
                end              
            end
        
            for id = 1:element_number
                a = [];
                b = [];
                c = [];
                d = [];
                e = [];
                for i = 1:4
                    xpt(i) = f(AB(id,i+1) + 1, 1);
                    ypt(i) = f(AB(id,i+1) + 1, 2);
                    zpt(i) = f(AB(id,i+1) + 1, 3);
                end
            
                for i = 1:4
                    a(i,1) = xpt(i);
                    a(i,2) = ypt(i);
                    a(i,3) = zpt(i);
                    a(i,4) = 1;
                end
                m11 = det(a);
            
                for i = 1:4    
                    b(i,1) = (xpt(i))^2 + (ypt(i))^2 + (zpt(i))^2;
                    b(i,2) = ypt(i);
                    b(i,3) = zpt(i);
                    b(i,4) = 1;  
                end
                m12 = det(b);
            
                for i = 1:4  
                    c(i,1) = (xpt(i))^2 + (ypt(i))^2 + (zpt(i))^2;
                    c(i,2) = xpt(i);
                    c(i,3) = zpt(i);
                    c(i,4) = 1; 
                end   
                m13 = det(c);
        
                for i = 1:4
                    d(i,1) = (xpt(i))^2 + (ypt(i))^2 + (zpt(i))^2;
                    d(i,2) = xpt(i);
                    d(i,3) = ypt(i);
                    d(i,4) = 1; 
                end   
                m14 = det(d);
        
                for i = 1:4   
                    e(i,1) = (xpt(i))^2 + (ypt(i))^2 + (zpt(i))^2;
                    e(i,2) = xpt(i);
                    e(i,3) = ypt(i);
                    e(i,4) = zpt(i); 
                end   
                m15 = det(e);
        
                if m11 == 0
                    r = 0;
                else
                    c1 = .5*m12/m11;
                    c2 = -0.5*m13/m11;
                    c3 = 0.5*m14/m11;
                    r(id) = sqrt(c1^2 + c2^2 + c3^2 - m15/m11);
                end
            end
        
            ravg(counter) = sum(r)/length(r);
        end
        rsum = sum(ravg)/length(ravg);
        fprintf(file2,'%d',rsum);
        cd .. % goes back up to original directory
    else
    end
