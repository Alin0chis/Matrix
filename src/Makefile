CC = g++ -std=c++17
CC_FLAGS = -Wall -Werror -Wextra -Wpedantic
VALGRIND_FLAGS = --leak-check=full --show-reachable=yes --show-leak-kinds=all --track-origins=yes
SOURCE = s21_matrix_oop.cpp
OBJECT = $(SOURCE:.cpp=.o)
OS := $(shell uname -s)
TEST_FLAGS =-lgtest --coverage
TEST_FILENAME = test
GCOV_OUTPUT = ./gcov/gcov_test
GCOV_FLAGS = -fprofile-arcs -ftest-coverage --coverage
LINKFLAGS = -lstdc++ -lm

COLOR_GREEN := \033[0;32m
COLOR_WHITE := \033[0m

ifeq ($(OS), Darwin)
	GTEST_FLAGS = -I/opt/homebrew/opt/googletest/include -L/opt/homebrew/opt/googletest/lib -lgtest -lgtest_main -pthread
	OPEN_CMD = open
	LCOV_FLAG = --ignore-errors inconsistent
else
	GTEST_FLAGS = -I/usr/include/gtest -L/usr/lib -lgtest -lgtest_main -pthread
	OPEN_CMD = xdg-open
	LCOV_FLAG = 
endif

all: clean s21_matrix_oop.a

s21_matrix_oop.a:$(SOURCE)
	@$(CC) $(CC_FLAGS) -c $(SOURCE)
	@ar rcs s21_matrix_oop.a $(OBJECT)

test: clean test.cpp s21_matrix_oop.a
	@$(CC) $(CC_FLAGS) $(TEST_FILENAME).cpp s21_matrix_oop.a $(TEST_FLAGS) -o $(TEST_FILENAME)
	@./$(TEST_FILENAME)

clang:
	clang-format -style=Google -n *.cpp  *.h

gcov_report:
	mkdir -p gcov
	$(CC) $(CC_FLAGS) test.cpp $(SOURCE) -o $(GCOV_OUTPUT) $(GTEST_FLAGS) $(GCOV_FLAGS) $(LINKFLAGS)
	./gcov/gcov_test
ifeq ($(OS), Darwin)
	mv *.gcno *.gcda ./gcov
endif
	lcov --capture --directory ./gcov --output-file ./gcov/coverage.info $(LCOV_FLAG)
	lcov --remove ./gcov/coverage.info '/opt/homebrew/opt/googletest/*' '/usr/*' --output-file ./gcov/coverage.info $(LCOV_FLAG)
	genhtml -o ./gcov/report ./gcov/coverage.info $(LCOV_FLAG)
	$(OPEN_CMD) ./gcov/report/index.html

leaks_test: clean test
ifeq ($(OS), Darwin)
	@MallocStackLogging=1 leaks -quiet -atExit -- ./$(TEST_FILENAME)
else
	@valgrind $(VALGRIND_FLAGS) --log-file=RESULT_VALGRIND.txt ./$(TEST_FILENAME)
	@echo "result log file in 'RESULT_VALGRIND.txt'"
endif

clean :
	@rm -rf *.o *.gcda *.gcno *.a test RESULT_VALGRIND.txt gcov
	@echo "${COLOR_GREEN}Successfully cleaned${COLOR_WHITE}"