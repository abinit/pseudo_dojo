#! /bin/bash



for file in `ls -1 *.fhi`
do
 newfilename=`echo $file | sed 's/6\.fhi/GGA\.fhi/'`

 mv $file $newfilename

done
