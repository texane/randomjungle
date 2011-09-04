./bin/bash

chmod 775 ./bin/rjungle
chmod 775 ./bin/rjunglesparse
rm demo/output/*
clear
echo Analyzing IRIS data set \(150 samples / 5 variables\):
echo --------------------------------------------------
echo - Grow 100 trees
echo - Calculate Gini Index Importance
echo - Use one thread \(1 CPU\)
echo
echo
echo
echo Press ENTER to continue!
read cd
./bin/rjungle -f demo/input/iris.dat -D Species -t 100 -y 1 -i 1 -v -o demo/output/run1
echo
echo Finished!
echo
echo Press ENTER to continue!
read cd
clear
echo Results:
echo
echo Confusion matrix
echo ----------------
cat demo/output/run1.confusion
echo
echo
echo
echo Importance values
echo ----------------
cat demo/output/run1.importance
echo
echo Press ENTER to continue!
read cd
clear


echo Analyzing GAW15 data set \(3500 samples / 9198 SNPs\):
echo --------------------------------------------------
echo - Use rjunglesparse program for genom-wide data
echo - Read _compressed_ PLINK PED file
echo - Grow 10 trees
echo - Use mtry = 96
echo - Calculate Permutation Importance
echo - Use 2 threads if available \(2 CPUs\)
echo
echo
echo
echo Press ENTER to continue!
read cd
./bin/rjunglesparse -f demo/input/gaw15.ped.gz -U2 -p -t 10  -m96 -y 1 -i 4 -v -o demo/output/run2
echo
echo
echo
echo Results in directory demo/output, file name run2.*
echo
echo Press ENTER to continue!
read cd




# GAW15 data)
#    Support for generation of the simulated data (GAW15) was provided from
#    NIH grants 5RO1-HL049609-14, 1R01-AG021917-01A1, the University of
#    Minnesota, and the Minnesota Supercomputing Institute. The authors also
#    acknowledge the GAW grant, R01-GM031575.
