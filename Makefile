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
CFLAGS	= $(OPT) -Wall -Werror -Wno-vla -fopenmp

LDFLAGS	= -lprimesieve


all: $(OUT)

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(DEFINES)

goldbach: goldbach.cpp $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(DEFINES)

.PHONY: all clean

clean:
	rm -f $(OUT) *.o
