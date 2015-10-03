###################################################
#                                                 #
#   LICHEM: Layered Interacting CHEmical Models   #
#                                                 #
#        Symbiotic Computational Chemistry        #
#                                                 #
###################################################

### Compiler settings ###

CXX=g++
CXXFLAGS=-static -fopenmp -O3
DEVFLAGS=-g -Wall
LDFLAGS=-I./src/ -I/usr/include/eigen3/

### Python settings ###

PYPATH=/usr/bin/python

### Manual settings ###

TEX=pdflatex
BIB=bibtex

###################################################

### Compile rules for users and devs ###

install:	title binary testexe manual compdone

Dev:	title devbin testexe manual stats compdone

clean:	title delbin compdone

###################################################

### Rules for building various parts of the code ###

binary:	
	@echo ""; \
	echo "### Compiling the LICHEM binary ###"
	$(CXX) $(CXXFLAGS) ./src/LICHEM.cpp -o lichem $(LDFLAGS)

devbin:	
	@echo ""; \
	echo "### Compiling the LICHEM development binary ###"
	$(CXX) $(CXXFLAGS) $(DEVFLAGS) ./src/LICHEM.cpp -o lichem $(LDFLAGS)

testexe:	
	@echo ""; \
	echo "### Creating test suite executable ###"
	@echo 'echo "#!$(PYPATH)" > ./tests/runtests'; \
	echo "!!$(PYPATH)" > ./tests/runtests; \
	echo "" >> ./tests/runtests
	cat ./src/runtests.py >> ./tests/runtests
	@sed -i 's/\#.*//g' ./tests/runtests; \
	sed -i 's/\s*$$//g' ./tests/runtests; \
	sed -i '/^$$/d' ./tests/runtests; \
	sed -i 's/\!\!/\#\!/g' ./tests/runtests; \
	chmod a+x ./tests/runtests

checksyntax:	
	@echo ""; \
	echo "### Checking for warnings and syntax errors ###"
	$(CXX) $(CXXFLAGS) $(DEVFLAGS) -fsyntax-only ./src/LICHEM.cpp -o lichem $(LDFLAGS)
	@echo ""; \
	echo "### Source code statistics ###"; \
	echo "Number of LICHEM source code files:"; \
	ls -al src/* | wc -l; \
	echo "Total length of LICHEM (lines):"; \
	cat src/* | wc -l; \

manual:	
	@echo ""; \
	echo "### Compiling the documentation ###"; \
	cd src/; \
	echo "$(TEX) manual"; \
	$(TEX) manual > doclog.txt; \
	$(BIB) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	echo "$(BIB) manual"; \
	$(BIB) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	mv manual.pdf ../doc/LICHEM_manual.pdf; \
	rm -f manual.aux manual.bbl manual.blg; \
	rm -f manual.log manual.out manual.toc; \
	rm -f doclog.txt

title:	
	@echo ""; \
	echo "###################################################"; \
	echo "#                                                 #"; \
	echo "#   LICHEM: Layered Interacting CHEmical Models   #"; \
	echo "#                                                 #"; \
	echo "#        Symbiotic Computational Chemistry        #"; \
	echo "#                                                 #"; \
	echo "###################################################"

stats:	
	@echo ""; \
        echo "### Source code statistics ###"; \
	echo "Number of LICHEM source code files:"; \
	ls -al src/* | wc -l; \
	echo "Total length of LICHEM (lines):"; \
	cat src/* | wc -l

compdone:	
	@echo ""; \
	echo "Done."; \
	echo ""

delbin:	
	@echo ""; \
	if grep -q "Jokes = 1" src/LICHEM_headers.h; then \
	echo '     ___'; \
	echo '    |_  |'; \
	echo '      \ \'; \
	echo '      |\ \'; \
	echo '      | \ \'; \
	echo '      \  \ \'; \
	echo '       \  \ \'; \
	echo '        \  \ \       <wrrr vroooom wrrr> '; \
	echo '         \__\ \________'; \
	echo '             |_________\'; \
	echo '             |__________|  ..,  ,.,. .,.,, ,..'; \
	echo ""; \
 	fi; \
        echo ""; \
	echo "Removing binary and manual..."; \
	rm -f lichem ./doc/LICHEM_manual.pdf ./tests/runtests
