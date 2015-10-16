# bcbio.prioritize -- Create executable script

all:
	rm -f target/*.jar
	lein uberjar
	cat bin/bcbio-prioritize.template target/*-standalone.jar > bin/bcbio-prioritize
	chmod +x bin/bcbio-prioritize
