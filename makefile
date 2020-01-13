all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes trenaMultiScore)

install:
	(cd ..; R CMD INSTALL --no-test-load trenaMultiScore)

check:
	(cd ..; R CMD check `ls -t trenaMultiScore) | head -1`)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

