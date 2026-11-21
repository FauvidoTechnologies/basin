# Releases use this repo's own PyPI token from ./.env (git-ignored), never ~/.pypirc:
# twine gets an empty config file, an explicit upload URL, and only the TWINE_* values
# set here, so other projects' credentials cannot be picked up. See .env.example.
ENV_FILE := $(CURDIR)/.env
VERSION := $(shell sed -n 's/^__version__ = "\(.*\)"/\1/p' basin/version.py)
# read KEY from .env, dropping quotes and spaces
token = $$(sed -n 's/^$(1) *= *//p' $(ENV_FILE) | tr -d "\"' ")
upload = env -u TWINE_REPOSITORY -u TWINE_REPOSITORY_URL -u TWINE_CERT \
	TWINE_USERNAME=__token__ TWINE_PASSWORD="$(call token,$(1))" \
	twine upload --non-interactive --config-file /dev/null --repository-url $(2) dist/*

.PHONY: test example clean build check-version check check-token test-release new-release

test:
	poetry run pytest tests/ -v

example:
	poetry run python example.py

clean:
	rm -rf dist build *.egg-info

check-version:
	@grep -q '^version = "$(VERSION)"' pyproject.toml || \
		{ echo "version mismatch: basin/version.py=$(VERSION), pyproject.toml differs"; exit 1; }

build: clean check-version
	poetry build

check: build
	twine check dist/*

check-token:
	@test -n "$(call token,$(TOKEN_VAR))" || \
		{ echo "no $(TOKEN_VAR) in $(ENV_FILE): see .env.example"; exit 1; }

test-release: TOKEN_VAR := testpypi_api_token
test-release: check
	@$(MAKE) --no-print-directory check-token TOKEN_VAR=$(TOKEN_VAR)
	@$(call upload,testpypi_api_token,https://test.pypi.org/legacy/)

new-release: TOKEN_VAR := pypy_api_token
new-release: check
	@$(MAKE) --no-print-directory check-token TOKEN_VAR=$(TOKEN_VAR)
	@$(call upload,pypy_api_token,https://upload.pypi.org/legacy/)
