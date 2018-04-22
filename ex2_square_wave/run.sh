#! /bin/bash

make -j;
#for n in 4 8 16 32 64 128; # 256 512 1024;
#for n in 10 20 40 80 160 320;
for n in 100;
do
./main $n 0 1;
done

