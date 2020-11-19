ID = 18nov2020
IMAGE_NAME = tms_navigator_$(ID)
TESTPORT=60055

build:
	docker build -t $(IMAGE_NAME) -f Dockerfile .

bash:
	docker run -it \
           -p 9010:3838 \
           --entrypoint /bin/bash \
           $(IMAGE_NAME)

test:
	sudo docker run -p 6000:3838 tms_navigator_18nov2020 R -e "shiny::runApp('tmsNavigator.R', port=3838, host='127.0.0.1')"


