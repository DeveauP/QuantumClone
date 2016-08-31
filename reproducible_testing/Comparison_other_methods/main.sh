#!/bin/bash

for param in depth CNA clones  contamination samples  mutations
	do
	echo $param
	bash scripts/parametered_testing.sh -p=$param
done
