FROM factor_app_09oct2020
WORKDIR /home/shiny
COPY installPrerequisites.R .
COPY tmsShinyNavigator.R .
COPY tbl.afp.RData .
#RUN R -f ./installPrerequisites.R .
CMD ["R", "-e", "shiny::runApp('/home/shiny/tmsShinyNavigator.R', port=3838, host='0.0.0.0')"]





