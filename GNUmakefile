all: sphinx

sphinx:
	sphinx-build -b html source html

sphinx-latex:
	sphinx-build -b latex -d build/doctrees   source build/latex

doc:
	epydoc --name biopy --url http://code.google.com/p/biopy/ -o html/  --config conf --inheritance=listed --docformat plaintext -v --graph all

sphinx-full:
#	rm -rf html/
	python setup.py install  --home=./`python --version 2>&1 | sed -s "s/Python /version-/"` > /dev/null 2>&1 
	sphinx-build -b html source html
