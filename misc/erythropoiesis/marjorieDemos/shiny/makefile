ID = 06apr2020
IMAGE_NAME = tms_brand_105_$(ID)
TESTPORT=60055

build:
	docker build -t $(IMAGE_NAME) -f Dockerfile .

bash:
	docker run -it \
           -v ~/github/TrenaMultiScore/misc/erythropoiesis/marjorieDemos/shiny/data:/home/shiny/data \
           -p 9010:3838 \
           --entrypoint /bin/bash \
           tms_brand_105_06apr2020:latest


shinyRun:
	docker run -it \
           -v ~/github/TrenaMultiScore/misc/erythropoiesis/marjorieDemos/shiny/data:/home/shiny/data \
           -p 9010:3838 \
           tms_brand_105_06apr2020:latest

	docker run -ti --rm -p $(TESTPORT):3838 rocker/shiny bash

dockerRun:
	docker run -p $(TESTPORT):3838 $(IMAGE_NAME)

save:
	time docker save $(IMAGE_NAME) >  ~/tmp/docker_$(IMAGE_NAME).image
	ls -lat ~/tmp/docker_$(IMAGE_NAME).image

scp:
	time scp ~/tmp/docker_$(IMAGE_NAME).image pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/docker_$(IMAGE_NAME).image
#	ssh pshannon@trena docker load < /local/tmp/docker_$(IMAGE_NAME).image


