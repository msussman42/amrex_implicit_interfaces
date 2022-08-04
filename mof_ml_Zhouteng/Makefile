PRESICION=DOUBLE_PRECISION

gendata_compile:
	gfortran -cpp -D$(PRESICION) -o gen_data.exec  generate_data.f90

gendata:
	./gen_data.exec
	python convert2binary.py

training:
	python training.py

predict_compile:
	gfortran -c -cpp -D$(PRESICION) -o ml.o ml.f90
	gfortran -c -cpp -D$(PRESICION) -o predict.o predict.f90
	gfortran ml.o predict.o -o predict.exec

predict:
	./predict.exec

all:
	gfortran -cpp -D$(PRESICION) -o gen_data.exec  generate_data.f90
	./gen_data.exec
	python convert2binary.py
	python training.py
	gfortran -c -cpp -D$(PRESICION) -o ml.o ml.f90
	gfortran -c -cpp -D$(PRESICION) -o predict.o predict.f90
	gfortran ml.o predict.o -o predict.exec
	./predict.exec

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.exec
	rm -f *.dat
	rm -f *.npy