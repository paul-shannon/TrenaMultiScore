NAME=klf1

default:
	@echo knit
	@echo publish
	@echo viewRemote

knit:
	R -e "rmarkdown::render('$(NAME).Rmd')"
	open $(NAME).html

publish:
	scp klf1.html  paulshannnon@pshannon.net:public_html/pshannon.net/tms/klf1/
	scp style.css  paulshannnon@pshannon.net:public_html/pshannon.net/tms/klf1/
	scp *.Rmd paulshannnon@pshannon.net:public_html/pshannon.net/tms/klf1/

viewRemote:
	open https://pshannon.net/tms/klf1/klf1.html

