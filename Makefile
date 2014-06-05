all: build

build:
	javac polyAlg/*.java distanceAlg1/*.java
	jar cfm gtp.jar manifest.txt polyAlg/*.class distanceAlg1/*.class

clean:
	rm polyAlg/*.class distanceAlg1/*.class