FROM bioshinyapp_18nov2020
WORKDIR /home/app
COPY *.R ./
COPY tbl.3.0.250000.500000-itraq.RData ./
 
 CMD ["R", "-e", "shiny::runApp('/home/shiny/tmsNavigator.R', port=3838, host='0.0.0.0')"]





