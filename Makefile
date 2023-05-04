format:
	python3 -m autopep8 -r --in-place ensemble_analyser
	python3 -m black ensemble_analyser

check:
		python3 -m autopep8 -rdv ensemble_analyser
