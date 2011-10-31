
doc:
	epydoc --name biopy --url http://code.google.com/p/biopy/ -o html/  --config conf --inheritance=listed --docformat plaintext -v --graph all

sphinx:
	sphinx-build -b html source build/html

sphinx-full:
#	rm -rf build/html/
	python setup.py install  --home=./`python --version 2>&1 | sed -s "s/Python /version-/"` > /dev/null 2>&1 
	sphinx-build -b html source build/html
