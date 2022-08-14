#PRESICION=DOUBLE_PRECISION
PRESICION=DEFAULT_IS_DOUBLE_PRECISION

gendata_compile:
	gfortran -D$(PRESICION) -o gen_data.exec  generate_data.F90

gendata:
	./gen_data.exec
	python convert2binary.py

training:
	python training.py

predict_compile:
	gfortran -c -D$(PRESICION) -o ml.o ml.F90
	gfortran -c -D$(PRESICION) -o predict.o predict.F90
	gfortran ml.o predict.o -o predict.exec

predict:
	./predict.exec

all:
	gfortran -D$(PRESICION) -o gen_data.exec  generate_data.F90
	./gen_data.exec
	python convert2binary.py
	python training.py
	gfortran -c -D$(PRESICION) -o ml.o ml.F90
	gfortran -c -D$(PRESICION) -o predict.o predict.F90
	gfortran ml.o predict.o -o predict.exec
	./predict.exec

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.exec
	rm -f *.dat
	rm -f *.npy
