default:
	@echo roxy install build test check all [roxy instll test]

all:  roxy install

roxy:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes TrenaMultiScore)

install:
	(cd ..; R CMD INSTALL --no-test-load TrenaMultiScore)

check:
	(cd ..; R CMD check `ls -t TrenaMultiScore) | head -1`)

test:
	(cd inst/unitTests; R -f test_TrenaMultiScore.R)


