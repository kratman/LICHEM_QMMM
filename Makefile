###################################################
#                                                 #
# FLUKE: Fields Layered Under Kohn-Sham Electrons #
#                                                 #
###################################################

CXX=g++
CXXFLAGS=-static -fopenmp -g -O3
LDFLAGS=-I/usr/include/eigen3/
TEX=pdflatex
BIB=bibtex

install:	title
	cd src/; \
	$(CXX) $(CXXFLAGS) FLUKE.cpp -o FLUKE $(LDFLAGS); \
	mv FLUKE ../.

all:	install manual

Dev:	install manual
	cd src/; \
	$(CC) -DDEVCOMP -o FLUKE; \
	mv FLUKE ../.

manual:	
	cd src/; \
	$(TEX) manual > doclog.txt; \
	$(BIB) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(BIB) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	mv manual.pdf ../doc/FLUKE_manual.pdf; \
	rm -f manual.aux manual.bbl manual.blg; \
	rm -f manual.log manual.out manual.toc; \
	rm -f doclog.txt;

title:	
	echo ""; \
	echo "###################################################"; \
	echo "#                                                 #"; \
	echo "# FLUKE: Fields Layered Under Kohn-Sham Electrons #"; \
	echo "#                                                 #"; \
	echo "###################################################"; \
	echo ""
