CC = gcc
CFLAGS = -Wall

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


all: enc dec

enc: $(OBJ_ENC)
	$(CC) $(CFLAGS) $(OBJ_ENC) -o $(BIN_DIR)/enc -lm

dec: $(OBJ_DEC)
	$(CC) $(CFLAGS) $(OBJ_DEC) -o $(BIN_DIR)/dec -lm


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -lm
	mv *.o $(OBJ_DIR)

