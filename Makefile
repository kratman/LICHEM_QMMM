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
	@echo "Compiling the FLUKE binary..."
	cd src/; \
	$(CXX) $(CXXFLAGS) FLUKE.cpp -o FLUKE $(LDFLAGS)
	@mv src/FLUKE .; \
	echo ""

all:	install manual

Dev:	install
	@echo "Compiling the FLUKE development binary..."
	cd src/; \
	$(CXX) $(CXXFLAGS) -DDEVCOMP FLUKE.cpp -o FLUKE $(LDFLAGS)
	@mv src/FLUKE .; \
	echo ""

DevAll:	Dev manual

manual:	
	@echo "Compiling the documentation..."; \
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
	@echo ""; \
	echo "###################################################"; \
	echo "#                                                 #"; \
	echo "# FLUKE: Fields Layered Under Kohn-Sham Electrons #"; \
	echo "#                                                 #"; \
	echo "###################################################"; \
	echo ""
