FROM rocker/shiny
WORKDIR /home/shiny
RUN mkdir /home/shiny/data
COPY installPrerequisites.R .
COPY exploreBindingSites.R .
RUN R -f ./installPrerequisites.R .
CMD ["R", "-e", "shiny::runApp('/home/shiny/exploreBindingSites.R', port=3838, host='0.0.0.0')"]





