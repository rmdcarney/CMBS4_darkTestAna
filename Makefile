# Simple Makefile for multiple executeables + X other files/classes
# Main files .cxx, Classes in .cpp and .h
# Define INCLUDE, LIBS and EXE accordingly

MAKEFLAGS=--warn-undefined-variables

INCLUDE = -I./. -I./include $(shell root-config --cflags --libs)

# Define compiler and flags
CXX = g++
CFLAGS = -g -O0 -Wall -std=c++11 $(INCLUDE) -Wno-psabi
LDFLAGS =  -L/usr/local/lib 

# Dir
SRC_DIR := src
OBJ_DIR := obj
TRG_DIR := exec
BIN_DIR := bin

SRC = $(wildcard $(SRC_DIR)/*.cpp) 
TARGET = $(wildcard $(TRG_DIR)/*.cxx)
OBJ = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC))
OOBJ = $(patsubst $(TRG_DIR)/%.cxx, $(OBJ_DIR)/%.o, $(TARGET))
EXE = $(patsubst $(TRG_DIR)/%.cxx, $(BIN_DIR)/%, $(TARGET))

.PHONY: all
all: $(EXE) $(OBJ) $(OOBJ)

$(OBJ_DIR)/%.o: $(TRG_DIR)/%.cxx $(OBJ)
	@$(CXX) -c $(CFLAGS) $< -o $@
	@echo "[Compiling] $@"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp 
	@$(CXX) -c $(CFLAGS) $< -o $@
	@echo "[Compiling] $@"

$(BIN_DIR)/%: $(OBJ_DIR)/%.o
	@$(CXX) $(CFLAGS) $(LDFLAGS) $< $(OBJ) -o $@ 
	@echo "[Linking] $@"

#If running with sc compile with: 'make sc'
.PHONY: sc
sc: CFLAGS += -DSC
sc: all

.PHONY: clean
clean:
	rm -f $(OBJ) $(OOBJ)
	rm -f $(EXE)

