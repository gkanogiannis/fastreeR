FASTREER_VERSION ?= $(shell git describe --tags --abbrev=0)

.PHONY: docker-build docker-tag docker-push docker-test \
        pypi-build pypi-upload check-version

# ----------------------
# Version consistency
# ----------------------
# Fails if FASTREER_VERSION in fastreeR.py differs from Version in DESCRIPTION.
# (setup.cfg is not checked: the pypi-publish workflow sets its version from the tag.)
check-version:
	@desc_ver=$$(sed -n 's/^Version:[[:space:]]*//p' DESCRIPTION); \
	py_ver=$$(sed -n 's/^FASTREER_VERSION[[:space:]]*=[[:space:]]*"\(.*\)".*/\1/p' fastreeR.py); \
	if [ -z "$$desc_ver" ]; then \
		echo "Could not read Version from DESCRIPTION." >&2; exit 1; \
	fi; \
	if [ "$$desc_ver" != "$$py_ver" ]; then \
		echo "Version mismatch: DESCRIPTION=$$desc_ver, fastreeR.py FASTREER_VERSION=$$py_ver" >&2; \
		echo "Update FASTREER_VERSION in fastreeR.py to match DESCRIPTION." >&2; \
		exit 1; \
	fi; \
	echo "Versions match: $$desc_ver"

# ----------------------
# Docker-related targets
# ----------------------
docker-build: check-version
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
PYPI_ROOT = deploy/fastreer-pypi
PACKAGE_DIR = ${PYPI_ROOT}/src/fastreer
SRC_JAR_DIR = inst/java
DEST_JAR_DIR = $(PACKAGE_DIR)/$(SRC_JAR_DIR)

# Clean build artifacts
pypi-clean:
	rm -rf ${PYPI_ROOT}/build ${PYPI_ROOT}/dist ${PYPI_ROOT}/**/*.egg-info
	rm -rf $(PACKAGE_DIR) $(PYPI_ROOT)/LICENSE.md $(PYPI_ROOT)/README.md

pypi-build: check-version pypi-clean
	@echo "Copying .jar files from $(SRC_JAR_DIR) to $(DEST_JAR_DIR)"
	@mkdir -p $(DEST_JAR_DIR)
	cp fastreeR.py $(PACKAGE_DIR)/cli.py
	cp $(SRC_JAR_DIR)/*.jar $(DEST_JAR_DIR)
	@mkdir -p $(PYPI_ROOT)
	cp LICENSE.md README.md $(PYPI_ROOT)
	@echo "Building Python package for version $(FASTREER_VERSION)"
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
GALAXY_ROOT = deploy/fastreer-galaxy
TOOLSHED_KEY ?= somekey

galaxy-update-versions:
	bash $(GALAXY_ROOT)/galaxy_update_versions.sh $(GALAXY_ROOT) $(FASTREER_VERSION)

galaxy-test-tools:
	-@bash $(GALAXY_ROOT)/galaxy_test_tools.sh $(GALAXY_ROOT)

galaxy-update-testshed:
	-@bash $(GALAXY_ROOT)/galaxy_update_shed.sh $(GALAXY_ROOT) testtoolshed $(TOOLSHED_KEY)

galaxy-update-toolshed:
	-@bash $(GALAXY_ROOT)/galaxy_update_shed.sh $(GALAXY_ROOT) toolshed $(TOOLSHED_KEY)
