CC = g++
INCLUDES = -I../..
CFLAGS = -g

define SRCI
../../types/filter/probabilistic_filter.cpp ../../types/data_loader/data_loader.cpp
endef

define SRCL 
../../core/filter/kalman_extended/kalman_extended.cpp ../../core/filter/kalman_unscented/kalman_unscented.cpp ../../core/input_filter/averaging_filter/averaging_filter.cpp ../../core/regularizer/boundary_projector/boundary_projector.cpp ../../core/updater/decay_updater/decay_updater.cpp
endef

define SRCU 
../../utilities/parsing/simple_parser/simple_parser.cpp ../../utilities/mathematics/fourier/fft23.cpp ../../utilities/mathematics/matrix/matrix.cpp ../../utilities/mathematics/special_functions/special_functions.cpp ../../utilities/mathematics/spline/spline.cpp
endef

define SRCM 
./model/polynomial_model.cpp
endef

all:
	$(CC) main.cpp $(CFLAGS) $(INCLUDES) $(SRCM) $(SRCL) $(SRCU) $(SRCI) -o pipeline