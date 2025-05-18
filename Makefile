FASTREER_VERSION ?= $(shell git describe --tags --abbrev=0)

PHONY: docker-build docker-tag docker-push docker-test \
        pypi-build pypi-upload

# ----------------------
# Docker-related targets
# ----------------------
docker-build:
	docker build -t fastreer .

docker-tag:
	docker tag fastreer gkanogiannis/fastreer:$(FASTREER_VERSION)

docker-push:
	docker push gkanogiannis/fastreer:$(FASTREER_VERSION)

docker-test:
	docker run --rm fastreer --check

# ----------------------
# PyPI-related targets
# ----------------------
PYPI_ROOT = fastreer-pypi
PACKAGE_DIR = ${PYPI_ROOT}/src/fastreer
SRC_JAR_DIR = inst/java
DEST_JAR_DIR = $(PACKAGE_DIR)/$(SRC_JAR_DIR)

# Clean build artifacts
pypi-clean:
	rm -rf ${PYPI_ROOT}/build ${PYPI_ROOT}/dist ${PYPI_ROOT}/**/*.egg-info
	rm -rf $(PACKAGE_DIR) $(PYPI_ROOT)/LICENSE.md $(PYPI_ROOT)/README.md

pypi-build: pypi-clean
	@echo "📦 Copying .jar files from $(SRC_JAR_DIR) to $(DEST_JAR_DIR)"
	@mkdir -p $(DEST_JAR_DIR)
	cp fastreeR.py $(PACKAGE_DIR)/cli.py
	cp $(SRC_JAR_DIR)/*.jar $(DEST_JAR_DIR)
	cp LICENSE.md README.md $(PYPI_ROOT)
	@echo "🔧 Building Python package for version $(FASTREER_VERSION)"
	@python -m build ${PYPI_ROOT}

# Local install
pypi-local-install: pypi-build
	pip install $(PYPI_ROOT)/dist/fastreer-*.whl --force-reinstall --no-cache-dir

# Upload to PyPI via twine (optional)
pypi-upload: pypi-build
	twine upload $(PYPI_ROOT)/dist/*

# ----------------------
# Galaxy-related targets
# ----------------------
GALAXY_ROOT = fastreer-galaxy

galaxy-update-versions:
	bash $(GALAXY_ROOT)/update_galaxy_versions.sh $(FASTREER_VERSION)
