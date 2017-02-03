#!/bin/bash

folderbase=STEP
filebase=box_0
n2=3499;
for i in `seq 0 $n2`
do

	if [ $i -lt 10 ]
	then 
		thisFolder=$folderbase\_0000$i	
	else
		if [ $i -lt 100 ]
        	then 
                	thisFolder=$folderbase\_000$i
		else
			if [ $i -lt 1000 ]
        		then
				thisFolder=$folderbase\_00$i
			else
				if [ $i -lt 10000 ]
				then
					thisFolder=$folderbase\_0$i	
				fi
			fi
		fi
		
        fi

	cd $thisFolder
	mv $filebase\_$i.res ../$filebase\_$i\_.res
	cd ..

done
