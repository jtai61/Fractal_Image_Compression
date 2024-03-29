CC = gcc
CFLAGS = -Wall -g

SRC_DIR = src
BIN_DIR = build
OBJ_DIR = $(BIN_DIR)/obj

OBJ_ENC = 	$(OBJ_DIR)/enc.o			\
			$(OBJ_DIR)/imgio.o			\
			$(OBJ_DIR)/indexing.o		\
			$(OBJ_DIR)/coding.o			\
			$(OBJ_DIR)/miscell.o		\
        	$(OBJ_DIR)/split.o			\
			$(OBJ_DIR)/nns.o

OBJ_DEC = 	$(OBJ_DIR)/dec.o			\
			$(OBJ_DIR)/miscell.o		\
			$(OBJ_DIR)/imgio.o

OBJ_EVAL =	$(OBJ_DIR)/eval.o			\
			$(OBJ_DIR)/miscell.o		\
			$(OBJ_DIR)/imgio.o


all: enc dec eval

enc: $(OBJ_ENC)
	$(CC) $(CFLAGS) $(OBJ_ENC) -o $(BIN_DIR)/enc -lm

dec: $(OBJ_DEC)
	$(CC) $(CFLAGS) $(OBJ_DEC) -o $(BIN_DIR)/dec -lm

eval: $(OBJ_EVAL)
	$(CC) $(CFLAGS) $(OBJ_EVAL) -o $(BIN_DIR)/eval -lm


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -lm
	mv *.o $(OBJ_DIR)

