
.PHONY: clean help clean-all plot
# $@ is the recipe name
# $< is the first dependence
# $^ are all the dependencies

# Directories:
INC_DIR     = ./include
SRC_DIR     = ./src
BUILD_DIR   = build
OUTDATA_DIR = output_data
OBJ_DIR     = $(BUILD_DIR)/obj
LIB_DIR     = lib

# Excecutables
CPP         = g++
C           = gcc
RM          = rm -f

# Flags
CPP_FLAGS   = -I$(INC_DIR) -O3 -std=c++17 -Warray-bounds -Wall -Wpedantic -Wextra
C_FLAGS     = -I$(INC_DIR) -O3
ASM_FLAGS   = -I$(INC_DIR)
LIBS        =  -lgsl -lgslcblas# -lboost_iostreams -lboost_system -lboost_filesystem
L_FLAGS     = $(LIBS)

# CPP_FLAGS += -DDEBUG

# Header files
# _DEPS       =
# DEPS        = $(patsubst %,$(INC_DIR)/%,$(_DEPS))

# Souce files
SRCS_C      = $(wildcard $(SRC_DIR)/*.c)
SRCS_CPP    = $(wildcard $(SRC_DIR)/*.cpp)
# SRCS_ASM    = $(wildcard $(SRC_DIR)/*.asm)
MAINS       = $(wildcard *.c) \
							$(wildcard *.cpp)
# MAINS      += $(wildcard *.asm)

# Objects
MAIN_OBJS   = \
	$(patsubst %.c,$(OBJ_DIR)/%.o,$(notdir $(MAINS))) \
	$(patsubst %.cpp,$(OBJ_DIR)/%.o,$(notdir $(MAINS)))
# MAIN_OBJS  += $(patsubst %.asm,$(OBJ_DIR)/%.o,$(MAINS))

OBJS        = \
	$(patsubst %.c,$(OBJ_DIR)/%.o,$(notdir $(SRCS_C))) \
	$(patsubst %.cpp,$(OBJ_DIR)/%.o,$(notdir $(SRCS_CPP))) \
	$(patsubst %.asm,$(OBJ_DIR)/%.o,$(notdir $(SRCS_ASM)))

# OBJS       += $(MAIN_OBJS)

# Attempt to create the output directories
ifneq ($(BUILD_DIR),)
$(shell [ -d $(BUILD_DIR) ] || mkdir -p $(BUILD_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
endif
ifneq ($(OUTDATA_DIR),)
$(shell [ -d $(OUTDATA_DIR) ] || mkdir -p $(OUTDATA_DIR))
endif


# Rules
Lang: $(OBJS) $(OBJ_DIR)/Langevien.o
	$(CPP) -o $@ $^ $(CPP_FLAGS) $(L_FLAGS)

# Mains rules
$(OBJ_DIR)/%.o: %.cpp
	$(CPP) -c -o $@ $< $(CPP_FLAGS)

$(OBJ_DIR)/%.o: %.c
	$(C) -c -o $@ $< $(C_FLAGS)

$(OBJ_DIR)/%.o: %.asm
	$(C) -c -o $@ $< $(ASM_FLAGS)

# Obj rules
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) -c -o $@ $< $(CPP_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(C) -c -o $@ $< $(C_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.asm
	$(C) -c -o $@ $< $(ASM_FLAGS)

help:
	@echo ciao

clean:
	$(RM) $(OBJ_DIR)/*.o

clean-all:
	$(RM) -r $(BUILD_DIR)

debug:
	@echo $(OBJ_PDF)
