CC = gcc
CFLAGS = -Wall -g

SRC_DIR = src
INC_DIR = $(SRC_DIR)/include
BIN_DIR = bin
OBJ_DIR = $(BIN_DIR)/obj

OBJ_ENC = 	$(OBJ_DIR)/mars_enc.o		\
			$(OBJ_DIR)/image_io.o		\
			$(OBJ_DIR)/index_func.o		\
			$(OBJ_DIR)/coding_func.o	\
			$(OBJ_DIR)/miscell.o		\
        	$(OBJ_DIR)/split_func.o		\
			$(OBJ_DIR)/nn_search.o		\

OBJ_DEC = 	$(OBJ_DIR)/mars_dec.o		\
			$(OBJ_DIR)/miscell.o		\
			$(OBJ_DIR)/image_io.o		\


all: encoder decoder

encoder: $(OBJ_ENC)
	$(CC) $(CFLAGS) -I$(INC_DIR) $(OBJ_ENC) -o $(BIN_DIR)/enc -lm

decoder: $(OBJ_DEC)
	$(CC) $(CFLAGS) -I$(INC_DIR) $(OBJ_DEC) -o $(BIN_DIR)/dec -lm

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -lm
	mv *.o $(OBJ_DIR)

