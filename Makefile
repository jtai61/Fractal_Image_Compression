CC = gcc
CFLAGS = -Wall -g

SRC_DIR = src
INC_DIR = $(SRC_DIR)/include
BIN_DIR = bin

OBJ_ENC = 	$(BIN_DIR)/enc.o			\
			$(BIN_DIR)/image_io.o		\
			$(BIN_DIR)/index_func.o		\
			$(BIN_DIR)/coding_func.o	\
			$(BIN_DIR)/miscell.o		\
        	$(BIN_DIR)/split_func.o		\
			$(BIN_DIR)/nn_search.o		\

OBJ_DEC = 	$(BIN_DIR)/dec.o			\
			$(BIN_DIR)/miscell.o		\
			$(BIN_DIR)/image_io.o		\

OBJ_EVAL =	$(BIN_DIR)/eval.o			\
			$(BIN_DIR)/miscell.o		\
			$(BIN_DIR)/image_io.o		\


all: enc dec eval

enc: $(OBJ_ENC)
	$(CC) $(CFLAGS) -I$(INC_DIR) $(OBJ_ENC) -o $(BIN_DIR)/enc -lm

dec: $(OBJ_DEC)
	$(CC) $(CFLAGS) -I$(INC_DIR) $(OBJ_DEC) -o $(BIN_DIR)/dec -lm

eval: $(OBJ_EVAL)
	$(CC) $(CFLAGS) -I$(INC_DIR) $(OBJ_EVAL) -o $(BIN_DIR)/eval -lm


$(BIN_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -lm
	mv *.o $(BIN_DIR)

