3L
=====
Implementation of the 3L algorithm
-----

![изображение](https://user-images.githubusercontent.com/59603419/114199557-151a8580-995d-11eb-8afd-c5564e42ca02.png)

## Compile
MacOS 11 Terminal build example: clang++ -Xpreprocessor -fopenmp -lomp -std=c++17 main.cpp

*****

## Tests
* to generate tests use one of these options:
1. add pictures tests/functions/contur_pictures in .tiff format and run tests/functions/contur.py (TODO: os.walk in that folder);
2. add function to tests/functions/fun_create.py/run_2();
3. run tests/functions/generate_tests.py to automaticly generate tests.

* to draw tests file run tests/draw_all.py


