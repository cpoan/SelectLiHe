all: main
main: main.cc
	g++ -g -o main1 main.cc -lEG `root-config --cflags` `root-config --libs` -lTreePlayer
clean:
	rm ./*~ ./*.o ./main
