#!/bin/bash

for i in snpdb variantdetails analysis upload reports; do
	echo "graphviz for '$i'";
	python3 manage.py graph_models $i > graphs/$i.dot;
	dot graphs/$i.dot -Tpng -o graphs/$i.png;
done


python3 manage.py graph_models snpdb variantdetails analysis upload reports > graphs/all.dot;
dot graphs/all.dot -Tpng -o graphs/all.png;


