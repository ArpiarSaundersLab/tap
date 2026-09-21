# Releasing TAP

This document gives a simple, repeatable release workflow for publishing the package to PyPI.

## 1. Prepare the release

Before you cut a release:

- Confirm the version is correct in [pyproject.toml](pyproject.toml).
- Confirm the changelog notes are ready.
- Confirm the README and install instructions match the package name and supported Python version.
- Confirm the package builds cleanly from a fresh environment.
- Confirm the tests pass for the release branch.

## 2. Update the version

Edit the package version in [pyproject.toml](pyproject.toml):

```toml
[project]
version = "0.0.2"
```

Use semantic versioning:

- `0.0.1` = initial release
- `0.1.0` = feature release
- `1.0.0` = stable release

## 3. Create a release branch or tag

On the release branch:

```bash
git checkout -b release-0.0.2
git add .
git commit -m "Prepare release 0.0.2"
git tag v0.0.2
```

## 4. Build the distribution artifacts

From the project root:

```bash
python3 -m pip install --upgrade pip build twine
python3 -m build
```

This creates:

- `dist/*.tar.gz` (source distribution)
- `dist/*.whl` (built wheel)

## 5. Validate the package locally

Before uploading, test the built wheel in a fresh environment:

```bash
python3 -m venv /tmp/tap-release-test
source /tmp/tap-release-test/bin/activate
pip install --upgrade pip
pip install dist/*.whl
python -c "import tap; print(tap.__file__)"
```

Optionally run a minimal smoke test:

```bash
python - <<'PY'
import tap
print('tap imported successfully')
print(getattr(tap, '__version__', 'no-version'))
PY
```

## 6. Upload to TestPyPI (recommended)

```bash
python3 -m twine upload --repository testpypi dist/*
```

Then test installation from TestPyPI:

```bash
python3 -m venv /tmp/tap-testpypi
source /tmp/tap-testpypi/bin/activate
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ tap
python -c "import tap; print('installed from TestPyPI')"
```

## 7. Upload to PyPI

Once the TestPyPI install works:

```bash
python3 -m twine upload dist/*
```

## 8. Push the release to GitHub

```bash
git push origin release-0.0.2
git push origin v0.0.2
```

Then create a GitHub release from the tag and add notes.

## 9. Final checks

- Confirm the package appears on PyPI.
- Confirm README renders correctly on PyPI.
- Confirm the install command is correct.
- Confirm the version is correct.
- Confirm the GitHub release is published.

## Recommended release checklist

- [ ] Version bumped in [pyproject.toml](pyproject.toml)
- [ ] README updated
- [ ] Tests run successfully
- [ ] Build succeeds
- [ ] Wheel installs from a clean environment
- [ ] Uploaded to TestPyPI
- [ ] Installed successfully from TestPyPI
- [ ] Uploaded to PyPI
- [ ] Git tag created
- [ ] GitHub release published
- [ ] Documentation links checked
