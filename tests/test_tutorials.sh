#!/bin/bash
cd ../tutorials
g++ -std=c++11 -march=native -o3 lattice_base_generation_tutorial.cpp -lfplll -lmpfr -lgmp -o test_base
g++ -std=c++11 -march=native -o3 lll_reduction_tutorial.cpp -lfplll -lmpfr -lgmp -o test__lll_reduction
g++ -std=c++11 -march=native -o3 bkz_reduction_tutorial.cpp -lfplll -lmpfr -lgmp -o test_bkz_reduction
mv test_base ../tests/lattices
mv test_bkz_reduction ../tests/lattices
mv test__lll_reduction ../tests/lattices
cd ../tests/lattices
rm -f base_output
rm -f lll_output
rm -f bkz_output
touch base_output
touch lll_output
touch bkz_output
./test_base
./test__lll_reduction
./test_bkz_reduction
diff bkz_output BKZ_reduction_evaluation