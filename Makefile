# press-dens
TARGET	= press-dens


CC 		= g++
CFLAGS	= -std=c++11 -O3 -Wall
LFLAGS	= -std=c++11 -O3 -Wall


SRC_DIR	= src
OBJ_DIR	= obj
BIN_DIR = bin


SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)


$(BIN_DIR)/$(TARGET): $(OBJS)
	@[ -d $(BIN_DIR) ] || mkdir $(BIN_DIR)
	$(CC) $(LFLAGS) $(OBJS) -o $@

$(OBJS): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	@[ -d $(OBJ_DIR) ] || mkdir $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@
	

.PHONY:	clean
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
