VERSION ?= latest

build:
	docker build -t fastreer .

tag:
	docker tag fastreer gkanogiannis/fastreer:$(VERSION)

push:
	docker push gkanogiannis/fastreer:$(VERSION)

test:
	docker run --rm fastreer --check
