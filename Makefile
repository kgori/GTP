objs = distanceAlg1/*.java polyAlg/*.java
classes = $(objs:java=class)

all: build jar

build:
	javac $(objs)

jar:
	jar cfm gtp.jar manifest.txt $(classes)

clean:
	@rm -f $(classes)
	@rm -f gtp.jar