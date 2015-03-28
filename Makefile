CXX=g++
CXXFLAGS=-static -fopenmp -g -O3
LDFLAGS=-I/usr/include/eigen3/
TEX=pdflatex
BIB=bibtex

install:	title manual
	echo "Compiling the FLUKE binary..."; \
	cd src/; \
	$(CXX) $(CXXFLAGS) FLUKE.cpp -o FLUKE $(LDFLAGS); \
	mv FLUKE ../.

Dev:	install manual
	cd src/; \
	$(CC) -DDEVCOMP -o FLUKE; \
	mv FLUKE ../.

manual:	
	echo "Compiling the documentation..."; \
	cd src/; \
	$(TEX) manual; \
	$(BIB) manual; \
	$(TEX) manual; \
	$(TEX) manual; \
	$(TEX) manual; \
	$(TEX) manual; \
	$(BIB) manual; \
	$(TEX) manual; \
	$(TEX) manual; \
	$(TEX) manual; \
	mv manual.pdf ../doc/FLUKE_manual.pdf

title:	
	echo ""; \
	echo "###################################################"; \
	echo "#                                                 #"; \
	echo "# FLUKE: Fields Layered Under Kohn-Sham Electrons #"; \
	echo "#                                                 #"; \
	echo "###################################################"; \
	echo ""
