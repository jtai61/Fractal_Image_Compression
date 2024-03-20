CC = g++
CPP_FLAGS = -std=c++17 -Wall -g

SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = obj

$(shell if [ ! -e $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi)
$(shell if [ ! -e $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi)

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
	$(CC) $(CPP_FLAGS) $(OBJ_ENC) -o $(BIN_DIR)/enc

dec: $(OBJ_DEC)
	$(CC) $(CPP_FLAGS) $(OBJ_DEC) -o $(BIN_DIR)/dec

eval: $(OBJ_EVAL)
	$(CC) $(CPP_FLAGS) $(OBJ_EVAL) -o $(BIN_DIR)/eval


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CPP_FLAGS) -c $<
	mv *.o $(OBJ_DIR)

