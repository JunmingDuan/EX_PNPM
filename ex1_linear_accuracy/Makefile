####################################################
# Generic makefile - 万能Makefile
# for compiling and linking C++ projects on Linux
# Author: George Foot  Modified:Jackie Lee
####################################################
### Customising
#
# Adjust the following if necessary; EXECUTABLE is the target
# executable's filename, and LIBS is a list of libraries to link in
# (e.g. alleg, stdcx, iostr, etc). You can override these on make's
# command line of course, if you prefer to do it that way.
#
#
EXECUTABLE := main    # 可执行文件名
LIBDIR:=              # 静态库目录
LIBS :=               # 静态库文件名
INCLUDES:=../include /home/jimmy/eigen        # 头文件目录
SRCDIR:=../src              # 除了当前目录外，其他的源代码文件目录
#
# # Now alter any implicit rules' variables if you like, e.g.:

CC:=g++
CFLAGS := -g -Wall -O3 -std=c++11 -fopenmp -lmx -leng
CPPFLAGS := $(CFLAGS)
CPPFLAGS += $(addprefix -I,$(INCLUDES))

RM-F := rm -f


SRCS := $(wildcard *.cpp) $(wildcard $(addsuffix /*.cpp, $(SRCDIR)))
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

.PHONY : all objs clean veryclean info

all: $(EXECUTABLE)

objs : $(OBJS)

clean :
	@$(RM-F) $(OBJS)
veryclean: clean
	@$(RM-F) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJS)
	$(CC) -o $(EXECUTABLE) $(OBJS) $(addprefix -L,$(LIBDIR)) $(addprefix -l,$(LIBS)) $(CFLAGS)

info:
	@echo $(SRCS)
	@echo $(OBJS)

