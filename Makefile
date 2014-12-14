all: qgen

OBJECTS = obj/qgen.o obj/qindivid.o obj/qobservstate.o obj/qrotoperator.o obj/mpicheck.o obj/mtrand.o
DFLAG =
CFLAG = -Iinclude/ -O3

qgen: $(OBJECTS)
	mkdir -p lib
	@ar rc lib/QGen.a $(OBJECTS)
	@echo -en "\033[30;1;41m QGen builded successfully! \033[0m\n"

obj/qgen.o: src/qgen.cpp
	mkdir -p obj
	@mpicxx -c $(CFLAG) $^ -o $@
	@echo -en "\033[30;1;46m $@ - done \033[0m\n"

obj/qindivid.o: src/qindivid.cpp
	@mkdir -p obj
	@mpicxx -c $(CFLAG) $^ -o $@
	@echo -en "\033[30;1;46m $@ - done \033[0m\n"

obj/qobservstate.o: src/qobservstate.cpp
	@mkdir -p obj
	@mpicxx -c $(CFLAG) $^ -o $@
	@echo -en "\033[30;1;46m $@ - done \033[0m\n"

obj/qrotoperator.o: src/qrotoperator.cpp
	@mkdir -p obj
	@mpicxx -c $(CFLAG) $^ -o $@
	@echo -en "\033[30;1;46m $@ - done \033[0m\n"

obj/mpicheck.o: src/mpicheck.cpp
	@mkdir -p obj
	@mpicxx -c $(CFLAG) $^ -o $@
	@echo -en "\033[30;1;46m $@ - done \033[0m\n"

obj/mtrand.o: src/mtrand.cpp
	@mkdir -p obj
	@mpicxx -c $(CFLAG) $^ -o $@
	@echo -en "\033[30;1;46m $@ - done \033[0m\n"

clean:
	rm -r -f lib
	rm -r -f obj