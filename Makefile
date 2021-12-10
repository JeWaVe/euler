clean:
	rm -f *.out

release:
	g++ --std=c++2a -O3 src/main.cpp -o app.out

tests:
	g++ -Isrc --std=c++2a -O3 tests/main.cpp -o tests.out

all: clean release tests

.PHONY: clean release tests