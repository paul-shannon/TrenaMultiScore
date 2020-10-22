ID = 22oct2020
IMAGE_NAME = tms_navigator_$(ID)
TESTPORT=60055

build:
	docker build -t $(IMAGE_NAME) -f Dockerfile .

bash:
	docker run -it \
           -v ~/github/TrenaMultiScore/misc/erythropoiesis/marjorieDemos/shiny/data:/home/shiny/data \
           -p 9010:3838 \
           --entrypoint /bin/bash \
           $(IMAGE_NAME)



