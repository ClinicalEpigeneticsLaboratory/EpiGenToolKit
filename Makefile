all: black autopep8

black:
	black .
autopep8:
	autopep8 -i --recursive .
