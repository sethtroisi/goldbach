# Copyright 2025 Seth Troisi
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

OPT     = -O3 -std=c++20 -g
OBJS	=
OUT	= goldbach
CC	= g++
CFLAGS	= $(OPT) -fopenmp -Wall -Werror -Wno-vla

LDFLAGS	= -lprimesieve


all: $(OUT)

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(DEFINES)

goldbach: goldbach.cpp $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

goldbach_gpu: goldbach_gpu.cu $(OBJS)
	nvcc -o $@ $^ $(OPT) $(LDFLAGS) \
		--generate-code arch=compute_61,code=sm_61 \
		-I CUDASieve/include/ -L CUDASieve -lcudasieve  \
		--compiler-options "-fopenmp -Wall -Werror -Wno-vla"

.PHONY: all clean

clean:
	rm -f $(OUT) *.o goldbach_gpu
