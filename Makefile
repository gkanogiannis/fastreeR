VERSION := $(shell git describe --tags --abbrev=0)

PHONY: docker-build docker-tag docker-push docker-test \
        pypi-sync-version pypi-build pypi-upload

# ----------------------
# Docker-related targets
# ----------------------
docker-build:
	docker build -t fastreer .

docker-tag:
	docker tag fastreer gkanogiannis/fastreer:$(VERSION)

docker-push:
	docker push gkanogiannis/fastreer:$(VERSION)

docker-test:
	docker run --rm fastreer --check

# ----------------------
# PyPI-related targets
# ----------------------
PACKAGE_DIR = fastreer
SRC_JAR_DIR = inst/java
DEST_JAR_DIR = $(PACKAGE_DIR)/inst/java

# Update setup.cfg version dynamically
pypi-sync-version:
	@echo "🛠️  Syncing version $(VERSION) into setup.cfg"
	@sed -i.bak "s/^version = .*/version = $(VERSION)/" setup.cfg && rm setup.cfg.bak

# Clean build artifacts
pypi-clean:
	rm -rf build dist *.egg-info
	rm -rf $(PACKAGE_DIR)/cli.py
	rm -rf $(DEST_JAR_DIR)/*.jar

# Build PyPI wheel
pypi-build: pypi-clean pypi-sync-version
	@echo "📦 Copying .jar files from $(SRC_JAR_DIR) to $(DEST_JAR_DIR)"
	@mkdir -p $(DEST_JAR_DIR)
	cp fastreeR.py $(PACKAGE_DIR)/cli.py
	cp $(SRC_JAR_DIR)/*.jar $(DEST_JAR_DIR)
	@echo "🔧 Building Python package for version $(VERSION)"
	@python -m build

# Local install
pypi-local-install: pypi-build
	pip install dist/fastreer-*.whl --force-reinstall --no-cache-dir

# Upload to PyPI via twine (optional)
pypi-upload:
	twine upload dist/*
