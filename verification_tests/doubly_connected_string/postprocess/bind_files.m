
filelist = readdir (pwd);
counter = 0;
%timestep / node / displacement
for ii = 1:numel(filelist)
        if (findstr (filelist{ii},"box"))
                counter = counter + 1;
                data(counter,:,:) = load(filelist{ii});
                splitName = strsplit (filelist{ii}, "_");
                timeArray(counter) = str2num(splitName{3});
        end
end

sizeDat = size(data);
saveMatrix = zeros(sizeDat(2),sizeDat(1));
for n = 1: sizeDat(2)
	for t = 1:sizeDat(1)
		saveMatrix(n,timeArray(t) + 1) = data(t,n,3);
			
		
	end
end


savename = ['combined.txt'];
save ("-ascii", savename, "saveMatrix");

