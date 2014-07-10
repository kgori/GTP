all: build

build:
	javac polyAlg/*.java distanceAlg1/*.java
	jar cfm gtp.jar manifest.txt polyAlg/*.class distanceAlg1/*.class

clean:
	@if [[ -f polyAlg/*.class ]]; then rm polyAlg/*.class; fi
	@if [[ -f distanceAlg1/*.class ]]; then rm distanceAlg1/*.class; fi